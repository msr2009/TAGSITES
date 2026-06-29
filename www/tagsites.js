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
  var plotYMax      = 1.1;  // fixed y ceiling for the line panel
  var legendHitBoxes = [];  // [{name, x0, y0, x1, y1}] — rebuilt on each render

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
  var TOP_GUTTER   = 6;    // top margin
  var SEQ_H        = 58;   // fixed height for sequence strip
  var LEGEND_H     = 22;   // fixed-height legend band between line and feature panels
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
    // fixed overhead: title + legend band + seq strip + 3 gaps
    var contentH  = Math.max(cssH - TOP_GUTTER - LEGEND_H - SEQ_H - PANEL_GAP * 3, 40);
    var lineH     = Math.round(contentH * LINE_FRAC);
    var featH     = contentH - lineH;
    var lineTop   = TOP_GUTTER;
    var legendTop = lineTop + lineH + PANEL_GAP;
    var featTop   = legendTop + LEGEND_H + PANEL_GAP;
    var seqTop    = featTop + featH + PANEL_GAP;
    return {
      lineTop:   lineTop,
      lineH:     lineH,
      legendTop: legendTop,
      legendH:   LEGEND_H,
      featTop:   featTop,
      featH:     featH,
      seqTop:    seqTop,
      seqH:      SEQ_H,
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

    drawLinePanel(ctx, inf, layout);
    drawLegendBand(ctx, inf, layout);
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

    // run name in top-left corner of the panel
    if (plotTitle) {
      ctx.save();
      ctx.font         = "bold 11px sans-serif";
      ctx.fillStyle    = "#555";
      ctx.textAlign    = "left";
      ctx.textBaseline = "top";
      ctx.fillText(plotTitle, LEFT_GUTTER + 5, top + 4);
      ctx.restore();
    }

    if (lineTracks.length === 0) return;

    // fixed y ceiling (plotYMax); auto-range only the minimum
    var yMin = Infinity;
    var yMax = plotYMax;
    lineTracks.forEach(function (tr) {
      if (hiddenTracks.has(tr.name)) return;
      tr.values.forEach(function (v) {
        if (v !== null && isFinite(v) && v < yMin) yMin = v;
      });
    });
    if (!isFinite(yMin)) yMin = 0;
    if (yMin >= yMax)    yMin = yMax - 1;
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

  var FEAT_ROWS = ["isoforms", "Phobius", "Pfam", "modification"];

  function drawFeaturePanel(ctx, inf, layout) {
    var top = layout.featTop, h = layout.featH;

    // panel background + border
    ctx.fillStyle   = "#fafafa";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_GUTTER, top, inf.dataW, h);

    // only show rows that actually have data
    var activeRows = FEAT_ROWS.filter(function (name) {
      return rangeFeatures.some(function (f) { return f.source === name; });
    });
    if (activeRows.length === 0) return;

    // assign evenly-spaced y positions within the panel
    var n = activeRows.length;
    var dynYPos = {};
    activeRows.forEach(function (name, idx) { dynYPos[name] = idx + 1; });
    var dynYMin = 0, dynYMax = n + 1;
    var ySpan   = dynYMax - dynYMin;

    // row labels and subtle separators
    ctx.textAlign    = "right";
    ctx.textBaseline = "middle";
    ctx.font         = "10px sans-serif";
    ctx.fillStyle    = "#777";
    activeRows.forEach(function (name) {
      var py = top + h - (dynYPos[name] - dynYMin) / ySpan * h;
      ctx.strokeStyle = "#ebebeb";
      ctx.lineWidth   = 0.5;
      ctx.beginPath();
      ctx.moveTo(LEFT_GUTTER, py);
      ctx.lineTo(LEFT_GUTTER + inf.dataW, py);
      ctx.stroke();
      ctx.fillText(name.charAt(0).toUpperCase() + name.slice(1), LEFT_GUTTER - 3, py);
    });

    var r     = getXRange();
    var bandH = Math.max(6, (h / ySpan) * 0.55);

    // clip feature rects to data area
    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.clip();

    rangeFeatures.forEach(function (feat) {
      var yv = dynYPos[feat.source];
      if (yv === undefined) return;
      if (feat.stop < r[0] || feat.start > r[1]) return;

      var x0 = posToX(feat.start - 0.5, inf);
      var x1 = posToX(feat.stop  + 0.5, inf);
      if (x1 <= x0 + 0.5) return;

      var py = top + h - (yv - dynYMin) / ySpan * h;

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
        var label = feat.desc.length > 14 ? feat.desc.slice(0, 13) + "…" : feat.desc;
        var fsize = Math.min(10, Math.max(7, w * 0.12));
        ctx.font         = fsize + "px sans-serif";
        ctx.fillStyle    = "#fff";
        ctx.textAlign    = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(label, (x0 + x1) / 2, py);
      }
    });

    ctx.restore();
  }

  /* ── Legend band (between line panel and feature panel) ─────────────────────── */

  function drawLegendBand(ctx, inf, layout) {
    var top = layout.legendTop;
    var h   = LEGEND_H;

    ctx.fillStyle = "#ffffff";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);

    legendHitBoxes = [];
    if (lineTracks.length === 0) return;

    var SWATCH_W = 14, SWATCH_H = 3, GAP = 4, ITEM_PAD = 10, LEFT_PAD = 8;
    var cy = top + h / 2;
    var cx = LEFT_GUTTER + LEFT_PAD;
    var rightEdge = LEFT_GUTTER + inf.dataW - LEFT_PAD;

    ctx.font         = "10px sans-serif";
    ctx.textBaseline = "middle";
    ctx.textAlign    = "left";

    lineTracks.forEach(function (tr) {
      var textW = ctx.measureText(tr.name).width;
      var itemW = SWATCH_W + GAP + textW + ITEM_PAD;
      if (cx + itemW > rightEdge) return;

      var hidden = hiddenTracks.has(tr.name);
      ctx.globalAlpha = hidden ? 0.3 : 1.0;

      // swatch as a short colored stroke
      ctx.strokeStyle = tr.color;
      ctx.lineWidth   = SWATCH_H;
      ctx.beginPath();
      ctx.moveTo(cx, cy);
      ctx.lineTo(cx + SWATCH_W, cy);
      ctx.stroke();

      // track name
      ctx.fillStyle = "#333";
      ctx.fillText(tr.name, cx + SWATCH_W + GAP, cy);

      ctx.globalAlpha = 1;

      legendHitBoxes.push({ name: tr.name, x0: cx, y0: top, x1: cx + itemW, y1: top + h });
      cx += itemW;
    });
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

    // clip residue drawing to the data area so partially-visible edge residues don't bleed into gutters
    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, top, inf.dataW, LETTER_H);
    ctx.clip();

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

    ctx.restore();  // end residue clip

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

  /* ── Legend (now drawn on canvas — buildLegend kept as no-op for safety) ─────── */

  function buildLegend() { /* legend is drawn in drawLegendBand(); nothing to do here */ }

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

    // no tooltip over the legend band
    if (cy >= layout.legendTop && cy <= layout.legendTop + LEGEND_H) { hideTooltip(); return; }

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

    // ── Drag tracked on document so it survives the mouse leaving the canvas ────
    // Private handler refs so we can remove them exactly.
    var _onDragMove = null;
    var _onDragUp   = null;

    canvas.addEventListener("mousedown", function (e) {
      if (e.button !== 0 && e.button !== 1) return;
      var rect = canvas.getBoundingClientRect();
      var cx   = e.clientX - rect.left;
      var r    = getXRange();
      wasDrag   = false;
      dragState = {
        startX:  cx,
        lastX:   cx,
        mode:    (e.shiftKey || e.button === 1) ? "pan" : dragMode,
        startR0: r[0],
        startR1: r[1],
      };
      e.preventDefault();

      _onDragMove = function (ev) {
        var rect2 = canvas.getBoundingClientRect();
        var mx    = ev.clientX - rect2.left;
        if (Math.abs(mx - dragState.startX) > 4) wasDrag = true;
        dragState.lastX = mx;
        if (wasDrag && dragState.mode === "pan") {
          var inf  = getCanvasInfo();
          if (inf) {
            var span   = dragState.startR1 - dragState.startR0;
            var dataDx = (dragState.startX - mx) / inf.dataW * span;
            currentRange = clampRange(dragState.startR0 + dataDx, dragState.startR1 + dataDx);
          }
        }
        if (wasDrag) render();
      };

      _onDragUp = function (ev) {
        document.removeEventListener("mousemove", _onDragMove);
        document.removeEventListener("mouseup",   _onDragUp);
        if (wasDrag && dragState && dragState.mode === "zoom") {
          var rect2 = canvas.getBoundingClientRect();
          var mx    = ev.clientX - rect2.left;
          var x0    = Math.min(dragState.startX, mx);
          var x1    = Math.max(dragState.startX, mx);
          var inf   = getCanvasInfo();
          if (inf && x1 - x0 > 4) {
            var r2   = getXRange();
            currentRange = clampRange(
              r2[0] + (x0 - LEFT_GUTTER) / inf.dataW * (r2[1] - r2[0]),
              r2[0] + (x1 - LEFT_GUTTER) / inf.dataW * (r2[1] - r2[0])
            );
          }
        }
        dragState = null;
        render();
      };

      document.addEventListener("mousemove", _onDragMove);
      document.addEventListener("mouseup",   _onDragUp);
    });

    // ── Tooltip + pointer cursor; drag is handled by document listeners above ───

    canvas.addEventListener("mousemove", function (e) {
      if (dragState) return;  // drag active — skip tooltip
      var rect   = canvas.getBoundingClientRect();
      var cy     = e.clientY - rect.top;
      var inf    = getCanvasInfo();
      var layout = getPanelLayout(inf ? inf.cssH : 520);
      var inLegend = inf && cy >= layout.legendTop && cy <= layout.legendTop + LEGEND_H;
      canvas.style.cursor = inLegend ? "pointer" : "crosshair";
      updateTooltip(e);
    });

    canvas.addEventListener("mouseleave", function () {
      if (!dragState) {
        hideTooltip();
        canvas.style.cursor = "crosshair";
      }
    });

    // ── Click: legend toggle or residue selection ────────────────────────────────

    canvas.addEventListener("click", function (e) {
      if (wasDrag) { wasDrag = false; return; }
      var rect   = canvas.getBoundingClientRect();
      var cx     = e.clientX - rect.left;
      var cy     = e.clientY - rect.top;
      var inf    = getCanvasInfo();
      var layout = getPanelLayout(inf ? inf.cssH : 520);

      // legend click → toggle track visibility
      if (cy >= layout.legendTop && cy <= layout.legendTop + LEGEND_H) {
        for (var i = 0; i < legendHitBoxes.length; i++) {
          var box = legendHitBoxes[i];
          if (cx >= box.x0 && cx <= box.x1) {
            if (hiddenTracks.has(box.name)) hiddenTracks.delete(box.name);
            else hiddenTracks.add(box.name);
            render();
            return;
          }
        }
        return;
      }

      if (!inf || !inDataArea(cx, inf)) return;
      shinySet(inputNames.residue_click, xToPos(cx, inf));
    });

    // ── Dblclick: seq strip → commit residue, elsewhere → reset zoom ────────────

    canvas.addEventListener("dblclick", function (e) {
      if (wasDrag) { wasDrag = false; return; }
      var rect   = canvas.getBoundingClientRect();
      var cx     = e.clientX - rect.left;
      var cy     = e.clientY - rect.top;
      var inf    = getCanvasInfo();
      var layout = getPanelLayout(inf ? inf.cssH : 520);

      // legend dblclick → no-op (single click already handles toggle)
      if (cy >= layout.legendTop && cy <= layout.legendTop + LEGEND_H) return;

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

    // Clear any existing viewer so stale canvases don't stack on re-init
    if (viewer) {
      viewer.clear();
      viewer = null;
    }
    container.innerHTML = "";

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

    viewer.resize();
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

  window.tsZoomIn = function () {
    var r    = getXRange();
    var mid  = (r[0] + r[1]) / 2;
    var half = (r[1] - r[0]) / 2 / 1.5;
    currentRange = clampRange(mid - half, mid + half);
    render();
  };

  window.tsZoomOut = function () {
    var r    = getXRange();
    var mid  = (r[0] + r[1]) / 2;
    var half = (r[1] - r[0]) / 2 * 1.5;
    currentRange = clampRange(mid - half, mid + half);
    render();
  };

  // Switch drag mode; clicking the active button reverts to default (zoom/set-region).
  window.tsSetMode = function (mode) {
    dragMode = (dragMode === mode && mode !== "zoom") ? "zoom" : mode;
    document.querySelectorAll(".ts-mode-btn").forEach(function (btn) {
      var onclick = btn.getAttribute("onclick") || "";
      var isMode  = onclick.indexOf("tsSetMode") !== -1;
      btn.classList.toggle("ts-mode-active", isMode && onclick.indexOf("'" + dragMode + "'") !== -1);
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
    plotYMax      = (typeof msg.yMax === "number") ? msg.yMax : 1.1;
    lineTracks    = msg.lineTracks    || [];
    rangeFeatures = msg.rangeFeatures || [];
    seqArray      = (msg.seq || "").split("");
    hiddenTracks  = new Set();
    pendingSet.clear();
    committedSet.clear();
    perResidueColors = [];
    currentRange     = null;
    dragMode         = "zoom";
    // reset toolbar — set-region (zoom) is the default active mode
    document.querySelectorAll(".ts-mode-btn").forEach(function (btn) {
      var onclick = btn.getAttribute("onclick") || "";
      btn.classList.toggle("ts-mode-active",
        onclick.indexOf("tsSetMode") !== -1 && onclick.indexOf("'zoom'") !== -1);
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
    // when no data colors are set (None), show the default N→C rainbow legend
    var legend = msg.legend || (perResidueColors.length === 0 ? {type: "rainbow"} : null);
    renderStructLegend(legend);
  });

  // ── Structure color legend ─────────────────────────────────────────────────────

  // viridis stops for the gradient bar (matches residue_colors_gradient in Python)
  var _VIRIDIS_CSS = "linear-gradient(to right," +
    "#440154,#48306e,#3e4989,#31688e,#26838f,#1f9e89,#35b779,#6ece58,#b5de2b,#fde725)";

  // rainbow N→C matches 3Dmol's default chain-spectrum coloring
  var _RAINBOW_CSS = "linear-gradient(to right,#0000ff,#00ffff,#00ff00,#ffff00,#ff0000)";

  function renderStructLegend(legend) {
    var el = document.getElementById("ts-struct-legend");
    if (!el) return;
    if (!legend) { el.innerHTML = ""; return; }

    if (legend.type === "gradient") {
      el.innerHTML =
        '<div class="ts-sleg-label">' + legend.label + '</div>' +
        '<div class="ts-sleg-bar-row">' +
          '<span class="ts-sleg-tick">' + legend.vmin + '</span>' +
          '<div class="ts-sleg-bar" style="background:' + _VIRIDIS_CSS + '"></div>' +
          '<span class="ts-sleg-tick">' + legend.vmax + '</span>' +
        '</div>';
    } else if (legend.type === "rainbow") {
      el.innerHTML =
        '<div class="ts-sleg-label">N → C</div>' +
        '<div class="ts-sleg-bar-row">' +
          '<span class="ts-sleg-tick">N</span>' +
          '<div class="ts-sleg-bar" style="background:' + _RAINBOW_CSS + '"></div>' +
          '<span class="ts-sleg-tick">C</span>' +
        '</div>';
    } else if (legend.type === "categorical") {
      var items = (legend.items || []).map(function (it) {
        return '<div class="ts-sleg-item">' +
          '<span class="ts-sleg-swatch" style="background:' + it.color + '"></span>' +
          '<span class="ts-sleg-name">' + it.label + '</span>' +
        '</div>';
      }).join("");
      el.innerHTML = '<div class="ts-sleg-cat">' + items + '</div>';
    }
  }

  Shiny.addCustomMessageHandler("tagsites_set_bg", function (msg) {
    setViewerBackground(msg.bg);
  });

  // Update task log textareas in-place so scroll position and user-resize are preserved.
  Shiny.addCustomMessageHandler("tagsites_update_logs", function (msg) {
    (msg.updates || []).forEach(function (u) {
      var el = document.getElementById(u.id);
      if (!el) return;
      // only auto-scroll if the user hasn't scrolled up
      var atBottom = el.scrollHeight - el.scrollTop <= el.clientHeight + 4;
      el.value = u.log;
      if (atBottom) el.scrollTop = el.scrollHeight;
    });
  });
})();
