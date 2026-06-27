/**
 * tagsites.js — client-side logic for the Results page
 *
 * Native canvas multi-panel plot: line tracks (row 1) + feature annotations (row 2)
 * + sequence strip (row 3), all sharing one x-axis with drag-to-zoom, pan, scroll
 * zoom, hover tooltips, and a clickable legend.
 *
 * Display hierarchy (sequence strip) by cell width:
 *   ≥ 9px  → bold ClustalX-colored letter (white when pending/committed)
 *   < 9px  → filled ClustalX color box (pending/committed override with amber/green)
 *
 * Click/dblclick behavior:
 *   click anywhere in data area  → toggle pending residue  (residue_click)
 *   dblclick in sequence strip   → commit residue          (residue_dblclick)
 *   dblclick in plot panels      → reset zoom to full range
 *
 * 3Dmol uses a timer-based single/double-click detector (no native dblclick).
 */

(function () {
  "use strict";

  /* ── Module state ──────────────────────────────────────────────────────────── */

  var viewer = null;
  var seqArray = [];
  var inputNames = {};
  var pendingSet = new Set();
  var committedSet = new Set();
  var perResidueColors = [];
  var currentBg = "transparent";

  // Plot data (sent from Python via tagsites_set_plot)
  var lineTracks    = [];   // [{name, color, values[]}]  values[i] = score at pos i+1, null = NaN
  var rangeFeatures = [];   // [{source, start, stop, desc, color, yRow}]
  var hiddenTracks  = new Set();
  var plotTitle     = "";

  // 3D-viewer click timer
  var structLastPos  = null;
  var structLastTime = 0;
  var DBL_CLICK_MS   = 350;

  // Zoom/pan state — null means full range
  var currentRange = null;
  var dragState    = null;  // {startX, lastX, mode, startR0, startR1}
  var wasDrag      = false;
  var dragMode     = "zoom";  // "zoom" or "pan" — toggled by toolbar buttons

  /* ── Layout constants (CSS px, drawn in DPR-scaled ctx) ───────────────────── */

  var LEFT_GUTTER  = 62;   // y-axis labels
  var RIGHT_GUTTER = 14;
  var TOP_GUTTER   = 28;   // title
  var SEQ_H        = 58;   // fixed height for sequence strip
  var PANEL_GAP    = 8;    // gap between panels
  // line panel gets 65% of the remaining height; feature panel gets the rest
  var LINE_FRAC    = 0.65;

  /* ── ClustalX residue text colors ──────────────────────────────────────────── */

  const CLUSTALX_COLORS = {
    A: "#33a02c", I: "#33a02c", L: "#33a02c", M: "#33a02c", V: "#33a02c",
    F: "#e07b00", Y: "#e07b00", W: "#e07b00",
    H: "#1f4fcc", K: "#1f4fcc", R: "#1f4fcc",
    D: "#cc2222", E: "#cc2222",
    S: "#9900cc", T: "#9900cc",
    N: "#008888", Q: "#008888",
    G: "#7a6000", P: "#7a6000", C: "#7a6000",
  };

  /* ── Helpers ───────────────────────────────────────────────────────────────── */

  function shinySet(inputId, value) {
    if (window.Shiny && window.Shiny.setInputValue) {
      Shiny.setInputValue(inputId, value, { priority: "event" });
    }
  }

  // Full x-range: [0.5, n+0.5] where n = sequence length.
  function fullRange() {
    var n = seqArray.length || 1;
    return [0.5, n + 0.5];
  }

  // Currently visible x-range (null → full range).
  function getXRange() {
    return currentRange || fullRange();
  }

  // Clamp [r0,r1] so r0 ≥ 0.5 and r1 ≤ n+0.5, preserving the span (pan stops at edge).
  // Also enforces a minimum span of 1 residue.
  function clampRange(r0, r1) {
    var fr = fullRange();
    var MIN = fr[0], MAX = fr[1];
    var span = Math.max(r1 - r0, 1.0);
    span = Math.min(span, MAX - MIN);
    if (r0 < MIN) { r0 = MIN; r1 = MIN + span; }
    if (r1 > MAX) { r1 = MAX; r0 = MAX - span; }
    if (r0 < MIN) r0 = MIN;
    return [r0, r1];
  }

  // Return the canvas element and its CSS dimensions + data-area width.
  function getCanvasInfo() {
    var canvas = document.getElementById("ts_plot_div");
    if (!canvas) return null;
    var cssW = canvas.clientWidth;
    var cssH = canvas.clientHeight;
    if (cssW <= 0 || cssH <= 0) return null;
    var dataW = Math.max(cssW - LEFT_GUTTER - RIGHT_GUTTER, 1);
    return { canvas: canvas, cssW: cssW, cssH: cssH, dataW: dataW };
  }

  // Data position (1-based) → canvas CSS x.
  function posToX(pos, inf) {
    var r = getXRange();
    return LEFT_GUTTER + (pos - r[0]) / (r[1] - r[0]) * inf.dataW;
  }

  // Canvas CSS x → nearest data position, clamped to [1,n].
  function xToPos(x, inf) {
    var n = seqArray.length || 1;
    var r = getXRange();
    var p = Math.round(r[0] + (x - LEFT_GUTTER) / inf.dataW * (r[1] - r[0]));
    return Math.max(1, Math.min(n, p));
  }

  // True if x is inside the data area (right of left gutter, left of right margin).
  function inDataArea(x, inf) {
    return x >= LEFT_GUTTER && x <= inf.cssW - RIGHT_GUTTER;
  }

  // Smallest "nice" interval ≥ raw from a fixed set.
  function niceInterval(raw) {
    var steps = [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000];
    for (var i = 0; i < steps.length; i++) {
      if (steps[i] >= raw) return steps[i];
    }
    return steps[steps.length - 1];
  }

  // Generate ~n nice y-tick values spanning [yMin, yMax].
  function niceYTicks(yMin, yMax, n) {
    var span = yMax - yMin;
    var rawStep = span / n;
    var mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
    var mults = [1, 2, 2.5, 5, 10];
    var step = mag * mults[mults.length - 1];
    for (var i = 0; i < mults.length; i++) {
      if (mag * mults[i] >= rawStep) { step = mag * mults[i]; break; }
    }
    var ticks = [];
    var start = Math.ceil(yMin / step) * step;
    for (var v = start; v <= yMax + 1e-9; v += step) {
      ticks.push(parseFloat(v.toPrecision(6)));
    }
    return ticks;
  }

  // Format a y-tick value compactly.
  function fmtVal(v) {
    if (!isFinite(v)) return "";
    if (Math.abs(v) >= 100) return v.toFixed(0);
    if (Math.abs(v) >= 1)   return v.toFixed(1);
    return v.toFixed(2);
  }

  // Compute panel boundary positions (all in CSS px, top-down).
  function getPanelLayout(cssH) {
    // content area sits below title and above fixed SEQ_H, with gaps between panels
    var contentH = Math.max(cssH - TOP_GUTTER - SEQ_H - PANEL_GAP * 2, 40);
    var lineH = Math.round(contentH * LINE_FRAC);
    var featH = contentH - lineH;
    return {
      lineTop: TOP_GUTTER,
      lineH:   lineH,
      featTop: TOP_GUTTER + lineH + PANEL_GAP,
      featH:   featH,
      seqTop:  TOP_GUTTER + lineH + PANEL_GAP + featH + PANEL_GAP,
      seqH:    SEQ_H,
    };
  }

  /* ── Core render ───────────────────────────────────────────────────────────── */

  function render() {
    var inf = getCanvasInfo();
    if (!inf) return;
    var canvas = inf.canvas;

    // resize backing store to match CSS size × DPR (only when size changes)
    var dpr = window.devicePixelRatio || 1;
    var bsW = Math.round(inf.cssW * dpr);
    var bsH = Math.round(inf.cssH * dpr);
    if (canvas.width !== bsW || canvas.height !== bsH) {
      canvas.width  = bsW;
      canvas.height = bsH;
    }

    var ctx = canvas.getContext("2d");
    ctx.save();
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, inf.cssW, inf.cssH);

    var layout = getPanelLayout(inf.cssH);

    drawTitle(ctx, inf, layout);
    drawLinePanel(ctx, inf, layout);
    drawFeaturePanel(ctx, inf, layout);
    drawSeqStrip(ctx, inf, layout);

    // drag-to-zoom rubber band overlaid on top
    if (dragState && wasDrag && dragState.mode === "zoom") {
      var rb0 = Math.min(dragState.startX, dragState.lastX);
      var rb1 = Math.max(dragState.startX, dragState.lastX);
      var rbTop  = layout.lineTop;
      var rbBot  = layout.seqTop + layout.seqH;
      ctx.fillStyle  = "rgba(100,149,237,0.13)";
      ctx.fillRect(rb0, rbTop, rb1 - rb0, rbBot - rbTop);
      ctx.strokeStyle = "rgba(100,149,237,0.65)";
      ctx.lineWidth   = 1;
      ctx.strokeRect(rb0, rbTop, rb1 - rb0, rbBot - rbTop);
    }

    ctx.restore();
  }

  /* ── Title ─────────────────────────────────────────────────────────────────── */

  function drawTitle(ctx, inf, layout) {
    if (!plotTitle) return;
    ctx.fillStyle    = "#333";
    ctx.font         = "bold 13px sans-serif";
    ctx.textAlign    = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(plotTitle, inf.cssW / 2, TOP_GUTTER / 2 + 1);
  }

  /* ── Line panel (row 1: continuous per-position scores) ──────────────────── */

  function drawLinePanel(ctx, inf, layout) {
    var top = layout.lineTop, h = layout.lineH;

    // panel background + border
    ctx.fillStyle   = "#fff";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_GUTTER, top, inf.dataW, h);

    // y-axis label "Score" (vertical, in gutter)
    ctx.save();
    ctx.translate(11, top + h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.font         = "11px sans-serif";
    ctx.textAlign    = "center";
    ctx.textBaseline = "middle";
    ctx.fillStyle    = "#666";
    ctx.fillText("Score", 0, 0);
    ctx.restore();

    if (lineTracks.length === 0) return;

    // compute shared y range across all visible tracks
    var yMin = Infinity, yMax = -Infinity;
    lineTracks.forEach(function (tr) {
      if (hiddenTracks.has(tr.name)) return;
      tr.values.forEach(function (v) {
        if (v !== null && isFinite(v)) {
          if (v < yMin) yMin = v;
          if (v > yMax) yMax = v;
        }
      });
    });
    if (!isFinite(yMin)) { yMin = 0; yMax = 1; }
    if (yMin === yMax)   { yMin -= 0.5; yMax += 0.5; }
    var ySpan = yMax - yMin;

    // y-axis grid lines and tick labels
    var yTicks = niceYTicks(yMin, yMax, 4);
    ctx.textAlign    = "right";
    ctx.textBaseline = "middle";
    ctx.font         = "10px sans-serif";
    ctx.fillStyle    = "#666";
    yTicks.forEach(function (yv) {
      var py = top + h - (yv - yMin) / ySpan * h;
      if (py < top || py > top + h) return;
      // faint grid
      ctx.strokeStyle = "#f0f0f0";
      ctx.lineWidth   = 1;
      ctx.beginPath();
      ctx.moveTo(LEFT_GUTTER, py);
      ctx.lineTo(LEFT_GUTTER + inf.dataW, py);
      ctx.stroke();
      // tick label
      ctx.fillText(fmtVal(yv), LEFT_GUTTER - 3, py);
    });

    // clip line traces to panel area
    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.clip();

    var r = getXRange();
    lineTracks.forEach(function (tr) {
      if (hiddenTracks.has(tr.name)) return;
      ctx.strokeStyle = tr.color;
      ctx.lineWidth   = 1.5;
      ctx.lineJoin    = "round";
      ctx.beginPath();
      var started = false;
      for (var i = 0; i < tr.values.length; i++) {
        var pos = i + 1;
        // skip positions far outside visible range (with 1-residue margin)
        if (pos < Math.floor(r[0]) - 1 || pos > Math.ceil(r[1]) + 1) { started = false; continue; }
        var v = tr.values[i];
        if (v === null || !isFinite(v)) { started = false; continue; }
        var px = posToX(pos, inf);
        var py = top + h - (v - yMin) / ySpan * h;
        if (!started) { ctx.moveTo(px, py); started = true; }
        else          { ctx.lineTo(px, py); }
      }
      ctx.stroke();
    });

    ctx.restore();  // remove clip
  }

  /* ── Feature panel (row 2: range annotations) ──────────────────────────────── */

  // Row positions match Python's y_positions = {"Phobius":2, "Pfam":4, "modification":6}
  var FEAT_ROWS = ["Phobius", "Pfam", "modification"];
  var FEAT_YVALS = { Phobius: 2, Pfam: 4, modification: 6 };
  var FEAT_YMIN  = 0, FEAT_YMAX = 8;  // axis range matching Python

  function drawFeaturePanel(ctx, inf, layout) {
    var top = layout.featTop, h = layout.featH;
    var ySpan = FEAT_YMAX - FEAT_YMIN;

    // panel background + border
    ctx.fillStyle   = "#fafafa";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_GUTTER, top, inf.dataW, h);

    // row labels and subtle separators
    ctx.textAlign    = "right";
    ctx.textBaseline = "middle";
    ctx.font         = "10px sans-serif";
    ctx.fillStyle    = "#777";
    FEAT_ROWS.forEach(function (name) {
      var yv  = FEAT_YVALS[name];
      var py  = top + h - (yv - FEAT_YMIN) / ySpan * h;
      // row separator
      ctx.strokeStyle = "#ebebeb";
      ctx.lineWidth   = 0.5;
      ctx.beginPath();
      ctx.moveTo(LEFT_GUTTER, py);
      ctx.lineTo(LEFT_GUTTER + inf.dataW, py);
      ctx.stroke();
      // label
      ctx.fillText(name, LEFT_GUTTER - 3, py);
    });

    if (rangeFeatures.length === 0) return;

    var r = getXRange();
    var rowPx   = (2 / ySpan) * h;          // pixels per 2-unit row spacing
    var bandH   = Math.max(6, rowPx * 0.55); // ~55% of row height, min 6px

    // clip feature rects to data area
    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.clip();

    rangeFeatures.forEach(function (feat) {
      var yv = FEAT_YVALS[feat.source];
      if (yv === undefined) return;
      // skip entirely off-screen features
      if (feat.stop < r[0] || feat.start > r[1]) return;

      var x0 = posToX(feat.start - 0.5, inf);
      var x1 = posToX(feat.stop  + 0.5, inf);
      if (x1 <= x0 + 0.5) return;

      var py = top + h - (yv - FEAT_YMIN) / ySpan * h;

      ctx.globalAlpha = 0.78;
      ctx.fillStyle   = feat.color;
      ctx.fillRect(x0, py - bandH / 2, x1 - x0, bandH);
      ctx.globalAlpha = 1;
      ctx.strokeStyle = feat.color;
      ctx.lineWidth   = 0.5;
      ctx.strokeRect(x0, py - bandH / 2, x1 - x0, bandH);

      // inline label when rect is wide enough
      var w = x1 - x0;
      if (w > 28) {
        var label  = feat.desc.length > 14 ? feat.desc.slice(0, 13) + "…" : feat.desc;
        var fsize  = Math.min(10, Math.max(7, w * 0.12));
        ctx.font         = fsize + "px sans-serif";
        ctx.fillStyle    = "#fff";
        ctx.textAlign    = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(label, (x0 + x1) / 2, py);
      }
    });

    ctx.restore();
  }

  /* ── Sequence strip (row 3) ────────────────────────────────────────────────── */

  function drawSeqStrip(ctx, inf, layout) {
    var top      = layout.seqTop, h = layout.seqH;
    var TICK_H   = 16;
    var LETTER_H = h - TICK_H;
    var LETTER_THRESHOLD = 9;

    // strip background + border + separator between letters and ticks
    ctx.fillStyle   = "#f8f8f8";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.beginPath();
    ctx.moveTo(LEFT_GUTTER, top + LETTER_H);
    ctx.lineTo(LEFT_GUTTER + inf.dataW, top + LETTER_H);
    ctx.stroke();

    if (seqArray.length === 0) return;

    var r     = getXRange();
    var cellW = inf.dataW / (r[1] - r[0]);
    var lo    = Math.max(1, Math.floor(r[0]));
    var hi    = Math.min(seqArray.length, Math.ceil(r[1]));

    for (var pos = lo; pos <= hi; pos++) {
      var aa   = seqArray[pos - 1];
      var cx   = posToX(pos, inf);
      var x0   = posToX(pos - 0.5, inf);
      var x1   = posToX(pos + 0.5, inf);
      var w    = x1 - x0;
      if (w < 0.5) continue;

      var isPending   = pendingSet.has(pos);
      var isCommitted = committedSet.has(pos);
      var baseColor   = CLUSTALX_COLORS[aa] || "#555555";

      if (cellW >= LETTER_THRESHOLD) {
        // letter mode
        var bgColor = isPending ? "#ffc107" : isCommitted ? "#28a745" : "#ffffff";
        ctx.fillStyle = bgColor;
        ctx.fillRect(x0, top, w, LETTER_H);
        ctx.strokeStyle = "#e0e0e0";
        ctx.lineWidth   = 0.5;
        ctx.strokeRect(x0, top, w, LETTER_H);
        var letterColor = (isPending || isCommitted) ? "#ffffff" : baseColor;
        var fontSize    = Math.min(14, Math.max(8, w * 0.75));
        ctx.fillStyle    = letterColor;
        ctx.font         = "bold " + fontSize + "px 'Courier New', monospace";
        ctx.textAlign    = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(aa, cx, top + LETTER_H / 2);
      } else {
        // box mode
        var boxColor = isPending ? "#ffc107" : isCommitted ? "#28a745" : baseColor;
        ctx.fillStyle = boxColor;
        ctx.fillRect(x0, top + 1, Math.max(w - 0.5, 0.5), LETTER_H - 2);
      }
    }

    // auto-spaced position ticks (labels ≥ 45px apart)
    var step = niceInterval(Math.ceil(45 / Math.max(cellW, 0.1)));
    ctx.fillStyle    = "#777";
    ctx.font         = "10px sans-serif";
    ctx.textAlign    = "center";
    ctx.textBaseline = "middle";
    var firstTick = Math.ceil(r[0] / step) * step;
    for (var t = firstTick; t <= r[1]; t += step) {
      if (t < 1 || t > seqArray.length) continue;
      var tx = posToX(t, inf);
      ctx.strokeStyle = "#bbb";
      ctx.lineWidth   = 1;
      ctx.beginPath();
      ctx.moveTo(tx, top + LETTER_H);
      ctx.lineTo(tx, top + LETTER_H + 4);
      ctx.stroke();
      ctx.fillText(String(t), tx, top + LETTER_H + TICK_H / 2 + 2);
    }
  }

  /* ── Legend ────────────────────────────────────────────────────────────────── */

  function buildLegend() {
    var wrap = document.getElementById("ts-plot-wrap");
    if (!wrap) return;
    // remove existing
    var existing = document.getElementById("ts-legend");
    if (existing) existing.parentNode.removeChild(existing);
    if (lineTracks.length === 0) return;

    var div = document.createElement("div");
    div.id = "ts-legend";

    lineTracks.forEach(function (tr) {
      var item = document.createElement("div");
      item.className = "ts-legend-item" + (hiddenTracks.has(tr.name) ? " ts-hidden" : "");
      item.setAttribute("data-track", tr.name);

      var swatch = document.createElement("span");
      swatch.className    = "ts-legend-swatch";
      swatch.style.background = tr.color;

      var label = document.createElement("span");
      label.textContent = tr.name;

      item.appendChild(swatch);
      item.appendChild(label);
      item.addEventListener("click", function () {
        var name = this.getAttribute("data-track");
        if (hiddenTracks.has(name)) hiddenTracks.delete(name);
        else hiddenTracks.add(name);
        this.classList.toggle("ts-hidden", hiddenTracks.has(name));
        render();
      });
      div.appendChild(item);
    });

    wrap.appendChild(div);
  }

  /* ── Tooltip ───────────────────────────────────────────────────────────────── */

  function ensureTooltip() {
    if (document.getElementById("ts-tooltip")) return;
    var wrap = document.getElementById("ts-plot-wrap");
    if (!wrap) return;
    var tip = document.createElement("div");
    tip.id = "ts-tooltip";
    wrap.appendChild(tip);
  }

  function showTooltip(x, y, lines) {
    var tip = document.getElementById("ts-tooltip");
    if (!tip) return;
    tip.textContent = lines.join("\n");
    tip.style.display = "block";
    // flip left if near right edge
    var wrap = document.getElementById("ts-plot-wrap");
    var wrapW = wrap ? wrap.offsetWidth : 600;
    var tipW  = tip.offsetWidth || 140;
    var left  = x + 14;
    if (left + tipW > wrapW - 8) left = x - tipW - 10;
    tip.style.left = left + "px";
    tip.style.top  = (y + 4) + "px";
  }

  function hideTooltip() {
    var tip = document.getElementById("ts-tooltip");
    if (tip) tip.style.display = "none";
  }

  function updateTooltip(e) {
    if (seqArray.length === 0 && lineTracks.length === 0) { hideTooltip(); return; }
    var canvas = document.getElementById("ts_plot_div");
    if (!canvas) return;
    var rect = canvas.getBoundingClientRect();
    var cx   = e.clientX - rect.left;
    var cy   = e.clientY - rect.top;
    var inf  = getCanvasInfo();
    if (!inf || !inDataArea(cx, inf)) { hideTooltip(); return; }

    var pos    = xToPos(cx, inf);
    var layout = getPanelLayout(inf.cssH);
    var lines  = ["Position: " + pos];

    // show track values when hovering the line panel
    if (cy >= layout.lineTop && cy <= layout.lineTop + layout.lineH) {
      lineTracks.forEach(function (tr) {
        if (hiddenTracks.has(tr.name)) return;
        var v = tr.values[pos - 1];
        if (v !== null && isFinite(v)) lines.push(tr.name + ": " + v.toFixed(3));
      });
    }

    // show feature info when hovering the feature panel
    if (cy >= layout.featTop && cy <= layout.featTop + layout.featH) {
      rangeFeatures.forEach(function (feat) {
        if (pos >= feat.start && pos <= feat.stop) {
          lines.push(feat.source + ": " + feat.desc + " (" + feat.start + "–" + feat.stop + ")");
        }
      });
    }

    showTooltip(cx, cy, lines);
  }

  /* ── Canvas interactions: zoom, pan, click, dblclick ──────────────────────── */

  function initCanvasInteractions(canvas) {

    // ── Drag-to-zoom / shift-drag-to-pan ───────────────────────────────────────

    canvas.addEventListener("mousedown", function (e) {
      if (e.button !== 0 && e.button !== 1) return;
      wasDrag = false;
      var rect = canvas.getBoundingClientRect();
      var cx   = e.clientX - rect.left;
      var r    = getXRange();
      dragState = {
        startX:  cx,
        lastX:   cx,
        mode:    (e.shiftKey || e.button === 1) ? "pan" : dragMode,
        startR0: r[0],
        startR1: r[1],
      };
      e.preventDefault();
    });

    canvas.addEventListener("mousemove", function (e) {
      var rect = canvas.getBoundingClientRect();
      var cx   = e.clientX - rect.left;

      if (dragState) {
        if (Math.abs(cx - dragState.startX) > 4) wasDrag = true;
        dragState.lastX = cx;

        if (wasDrag && dragState.mode === "pan") {
          var inf  = getCanvasInfo();
          if (inf) {
            var span   = dragState.startR1 - dragState.startR0;
            var dataDx = (dragState.startX - cx) / inf.dataW * span;
            currentRange = clampRange(dragState.startR0 + dataDx, dragState.startR1 + dataDx);
          }
        }
        if (wasDrag) render();
      }

      updateTooltip(e);
    });

    canvas.addEventListener("mouseup", function (e) {
      if (dragState && wasDrag && dragState.mode === "zoom") {
        var rect = canvas.getBoundingClientRect();
        var cx   = e.clientX - rect.left;
        var x0   = Math.min(dragState.startX, cx);
        var x1   = Math.max(dragState.startX, cx);
        var inf  = getCanvasInfo();
        if (inf && x1 - x0 > 4) {
          var r    = getXRange();
          var r0new = r[0] + (x0 - LEFT_GUTTER) / inf.dataW * (r[1] - r[0]);
          var r1new = r[0] + (x1 - LEFT_GUTTER) / inf.dataW * (r[1] - r[0]);
          currentRange = clampRange(r0new, r1new);
          render();
        }
      }
      dragState = null;
    });

    canvas.addEventListener("mouseleave", function () {
      // don't cancel dragState here — user may re-enter; just hide tooltip
      hideTooltip();
    });

    // scroll-wheel zoom centered on cursor position
    canvas.addEventListener("wheel", function (e) {
      e.preventDefault();
      var rect   = canvas.getBoundingClientRect();
      var cx     = e.clientX - rect.left;
      var inf    = getCanvasInfo();
      if (!inf || !inDataArea(cx, inf)) return;
      var r      = getXRange();
      var factor = e.deltaY > 0 ? 1.18 : 1 / 1.18;
      var pivot  = r[0] + (cx - LEFT_GUTTER) / inf.dataW * (r[1] - r[0]);
      currentRange = clampRange(pivot - (pivot - r[0]) * factor, pivot + (r[1] - pivot) * factor);
      render();
    }, { passive: false });

    // ── Click / dblclick (residue selection) ────────────────────────────────────
    // Two rapid clicks cancel each other (toggle-on then toggle-off) before dblclick
    // fires — net state is unchanged, then dblclick commits or resets zoom cleanly.

    canvas.addEventListener("click", function (e) {
      if (wasDrag) { wasDrag = false; return; }
      var rect = canvas.getBoundingClientRect();
      var cx   = e.clientX - rect.left;
      var inf  = getCanvasInfo();
      if (!inf || !inDataArea(cx, inf)) return;
      shinySet(inputNames.residue_click, xToPos(cx, inf));
    });

    canvas.addEventListener("dblclick", function (e) {
      if (wasDrag) { wasDrag = false; return; }
      var rect   = canvas.getBoundingClientRect();
      var cx     = e.clientX - rect.left;
      var cy     = e.clientY - rect.top;
      var inf    = getCanvasInfo();
      var layout = getPanelLayout(inf ? inf.cssH : 520);

      // dblclick in the sequence strip → commit/remove residue
      // dblclick elsewhere (plot panels, gutter) → reset zoom to full range
      if (cy >= layout.seqTop && inf && inDataArea(cx, inf)) {
        shinySet(inputNames.residue_dblclick, xToPos(cx, inf));
      } else {
        currentRange = null;
        render();
      }
    });
  }

  window.addEventListener("resize", function () { render(); });

  /* ── 3Dmol viewer ──────────────────────────────────────────────────────────── */

  function initViewer(pdbStr) {
    var container = document.getElementById("ts-viewer-container");
    if (!container) return;

    viewer = $3Dmol.createViewer(container, {
      backgroundColor: "white", backgroundAlpha: 0,
      hoverDuration: 100,
      antialias: true,
    });
    viewer.addModel(pdbStr, "pdb");
    viewer.setStyle({}, { cartoon: { color: "spectrum" } });

    // single/double-click via timer (3Dmol has no native dblclick)
    viewer.setClickable({}, true, function (atom) {
      if (!atom || atom.resi == null) return;
      var pos = parseInt(atom.resi);
      var now = Date.now();
      if (structLastPos === pos && now - structLastTime < DBL_CLICK_MS) {
        structLastPos = null;
        shinySet(inputNames.residue_dblclick, pos);
      } else {
        structLastPos  = pos;
        structLastTime = now;
        shinySet(inputNames.residue_click, pos);
      }
    });

    viewer.setHoverable({}, true,
      function (atom) {
        if (atom._label) return;
        atom._label = viewer.addLabel(
          atom.resn + " " + atom.resi + (atom.chain ? " " + atom.chain : ""),
          { position: atom, backgroundColor: "rgba(0,0,0,0.75)",
            fontColor: "white", fontSize: 12, borderThickness: 0 }
        );
        viewer.render();
      },
      function (atom) {
        if (atom._label) {
          viewer.removeLabel(atom._label);
          delete atom._label;
          viewer.render();
        }
      }
    );

    viewer.zoomTo();
    viewer.render();
  }

  function setViewerBackground(bg) {
    currentBg = bg;
    document.querySelectorAll(".ts-bg-btn").forEach(function (btn) {
      var isBg = btn.getAttribute("onclick") || "";
      btn.classList.toggle("ts-bg-btn-active", isBg.indexOf("'" + bg + "'") !== -1);
    });
    if (!viewer) return;
    if (bg === "transparent") {
      viewer.setBackgroundColor(0xffffff, 0);
    } else {
      viewer.setBackgroundColor(bg);
    }
    viewer.render();
  }

  function applyViewerColors() {
    if (!viewer) return;

    if (perResidueColors.length === 0) {
      viewer.setStyle({}, { cartoon: { color: "spectrum" } });
    } else {
      viewer.setStyle({}, { cartoon: {} });
      for (var i = 0; i < perResidueColors.length; i++) {
        viewer.setStyle({ resi: i + 1 }, { cartoon: { color: perResidueColors[i] || "#888888" } });
      }
    }

    pendingSet.forEach(function (pos) {
      viewer.setStyle({ resi: pos }, { cartoon: { color: "#ffc107" } });
      viewer.addStyle({ resi: pos }, { stick: { colorscheme: "amino", radius: 0.2 } });
    });
    committedSet.forEach(function (pos) {
      viewer.setStyle({ resi: pos }, { cartoon: { color: "#28a745" } });
      viewer.addStyle({ resi: pos }, { stick: { colorscheme: "amino", radius: 0.2 } });
    });

    viewer.render();
  }

  /* ── Exposed globals ─────────────────────────────────────────────────────────── */

  window.tsSetMode = function (mode) {
    dragMode = mode;
    document.querySelectorAll(".ts-mode-btn").forEach(function (btn) {
      var isTarget = (btn.getAttribute("onclick") || "").indexOf("'" + mode + "'") !== -1;
      btn.classList.toggle("ts-mode-active", isTarget);
    });
  };

  window.tsResetZoom = function () {
    currentRange = null;
    render();
  };

  window.tsSetColorBy = function (btn, value, inputId) {
    document.querySelectorAll(".ts-colorby-btn").forEach(function (b) {
      b.classList.remove("ts-colorby-active");
    });
    btn.classList.add("ts-colorby-active");
    shinySet(inputId, value);
  };

  window.tsSetBackground = function (bg) { setViewerBackground(bg); };

  window.tsDownloadStructure = function () {
    if (!viewer) return;
    var container = document.getElementById("ts-viewer-container");
    var origW = container.offsetWidth;
    var origH = container.offsetHeight;
    var SCALE = 3;

    container.style.width  = (origW * SCALE) + "px";
    container.style.height = (origH * SCALE) + "px";
    viewer.resize();
    viewer.render();

    setTimeout(function () {
      var uri = viewer.getCanvas().toDataURL("image/png");
      container.style.width  = origW + "px";
      container.style.height = origH + "px";
      viewer.resize();
      viewer.render();
      var a = document.createElement("a");
      a.href = uri;
      a.download = "structure.png";
      a.click();
    }, 120);
  };

  window.tsRemoveSite = function (pos, inputId) {
    Shiny.setInputValue(inputId, pos, { priority: "event" });
  };

  /* ── Custom message handlers (Python → JS) ─────────────────────────────────── */

  // Main initialization: receives title, lineTracks, rangeFeatures, seq, inputs
  Shiny.addCustomMessageHandler("tagsites_set_plot", function (msg) {
    inputNames    = msg.inputs || {};
    plotTitle     = msg.title  || "";
    lineTracks    = msg.lineTracks    || [];
    rangeFeatures = msg.rangeFeatures || [];
    seqArray      = (msg.seq || "").split("");
    hiddenTracks  = new Set();
    pendingSet.clear();
    committedSet.clear();
    perResidueColors = [];
    currentRange     = null;
    dragMode         = "zoom";
    // reset toolbar active state to zoom
    document.querySelectorAll(".ts-mode-btn").forEach(function (btn) {
      btn.classList.toggle("ts-mode-active",
        (btn.getAttribute("onclick") || "").indexOf("'zoom'") !== -1);
    });

    // attach interactions exactly once per canvas element
    var canvas = document.getElementById("ts_plot_div");
    if (canvas && !canvas._tsInteractive) {
      canvas._tsInteractive = true;
      initCanvasInteractions(canvas);
    }
    ensureTooltip();
    buildLegend();
    render();
  });

  Shiny.addCustomMessageHandler("tagsites_init_struct", function (msg) {
    inputNames = msg.inputs;
    var pane = document.getElementById("ts-struct-pane");
    if (pane) pane.style.display = "flex";
    initViewer(msg.pdb);
  });

  Shiny.addCustomMessageHandler("tagsites_set_states", function (msg) {
    pendingSet   = new Set(msg.pending);
    committedSet = new Set(msg.committed);
    render();
    applyViewerColors();
  });

  Shiny.addCustomMessageHandler("tagsites_set_colors", function (msg) {
    perResidueColors = msg.colors || [];
    applyViewerColors();
  });

  Shiny.addCustomMessageHandler("tagsites_set_bg", function (msg) {
    setViewerBackground(msg.bg);
  });
})();
