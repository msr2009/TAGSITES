/**
 * tagsites.js — client-side logic for the Results page
 *
 * Native canvas multi-panel plot, top to bottom: line tracks (continuous scores),
 * feature annotations, isoforms (one row per isoform, hidden when <= 1), score
 * heatmap, and sequence strip. All panes share one x-axis with drag-to-zoom, pan,
 * scroll zoom, hover tooltips, and a clickable legend. The line/annotations/isoform
 * panes are content-driven with minimum heights; the score row and sequence strip
 * are fixed height.
 *
 * Display hierarchy (sequence strip) by cell width:
 *   ≥ 9px  → bold ClustalX-colored letter (white when committed)
 *   < 9px  → filled ClustalX color box (committed override with green)
 *
 * Click/dblclick behavior:
 *   click anywhere in data area  → add residue to committed sites (residue_click)
 *   dblclick in plot panels      → reset zoom to full range
 *   click on 3D structure        → add residue to committed sites (struct_click)
 */

(function () {
  "use strict";

  /* ── Module state ──────────────────────────────────────────────────────────── */

  var viewer = null;
  var seqArray = [];
  var inputNames = {};
  var committedSet = new Set();
  var perResidueColors = [];
  var currentBg = "transparent";

  // Plot data (sent from Python via tagsites_set_plot)
  var lineTracks    = [];   // [{name, color, values[]}]  values[i] = score at pos i+1, null = NaN
  var rangeFeatures = [];   // [{source, start, stop, desc, color, yRow}]
  var hiddenTracks  = new Set();
  var plotTitle     = "";
  var plotYMax      = 1.0;  // fixed y ceiling for the line panel
  var legendHitBoxes = [];  // [{name, x0, y0, x1, y1}] — rebuilt on each render

  // Isoform pane (see utils/results.py build_plot_payload / _build_isoform_pane)
  var isoIsoforms     = [];  // [{label, accession, isQuery, present, skipped, inserts, domains}]

  // isoform label plus its UniProt accession, when known and not already the label itself
  function isoLabelWithAcc(iso) {
    return iso.label + (iso.accession && iso.accession !== iso.label ? " (" + iso.accession + ")" : "");
  }
  var isoExonBoundaries = [];  // [pos, ...] — reserved for future genewise integration; empty today
  var isoHitBoxes      = [];  // [{isoIdx, x, y, r}] — insert-marker hit targets, rebuilt on each render

  // Tag-site score heatmap row (see utils/scoring.py, issue #32)
  var scoreTrack = [];   // scoreTrack[i] = int score at pos i+1, null = no data
  var scoreFlags = [];   // scoreFlags[i] = [failed-criterion labels] at pos i+1
  var scoreMasked = [];       // scoreMasked[i] = true if pos i+1 is masked outright (issue #38)
  var scoreMaskReasons = [];  // scoreMaskReasons[i] = [mask-criterion labels] at pos i+1
  var scoreMax   = 1;    // normalization ceiling for the white→green fill
  var MASKED_HEX = "#d9b3b3";  // must match modules.results_server._MASKED_HEX
  var suggestedSites = [];  // positions picked by "Add suggested" — marked with * above the heatmap

  // Zoom/pan state — null means full range
  var currentRange = null;
  var dragState    = null;  // {startX, lastX, mode, startR0, startR1}
  var wasDrag      = false;
  var dragMode     = "zoom";  // "zoom" or "pan" — toggled by toolbar buttons

  /* ── Layout constants (CSS px, drawn in DPR-scaled ctx) ───────────────────── */

  var LEFT_GUTTER  = 62;   // y-axis labels
  var RIGHT_GUTTER = 14;
  var TOP_GUTTER   = 6;    // top margin
  var SCORE_FILL_H = 24;   // score heatmap bar height (excludes marker margin)
  var SEQ_TICK_H   = 16;   // axis-number row above the sequence letters
  var SEQ_H        = SCORE_FILL_H + SEQ_TICK_H;   // letters (= score bar height) + axis numbers
  var SCORE_MARK_H = 8;    // top margin reserved for suggested-site triangle markers
  var SCORE_H      = SCORE_FILL_H + SCORE_MARK_H;   // total row height (fill + marker margin)
  var LEGEND_H     = 22;   // fixed-height legend band between line and feature panels
  var PANEL_GAP    = 8;    // gap between panels

  // The line, annotations, and isoform panes are content-driven with minimum heights;
  // only the legend/score/sequence rows above have fixed heights.
  var MIN_LINE_H   = 80;   // line panel minimum (continuous tracks)
  var LINE_TRACK_H = 40;   // notional per-track height, scales the line panel with track count
  var MIN_FEAT_H   = 80;   // annotations panel minimum (range features)
  var FEAT_ROW_H   = 34;   // per-active-row height inside the annotations panel
  var ISO_ROW_H    = 26;   // per-isoform row height inside the isoform pane
  var ISO_PAD      = 10;   // top+bottom padding inside the isoform pane
  var PANE_LABEL_H = 18;   // space reserved above each dynamic pane for its title

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

  // Wrap a residue position in a fresh array so re-clicking the SAME position
  // still invalidates server-side: py-shiny's reactive Value skips updates
  // when the new value `is` the old one, and CPython caches small ints, so a
  // bare repeated int silently fails to re-trigger the click handler.
  function shinySetPos(inputId, pos) {
    shinySet(inputId, [pos]);
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

  // Format a y-tick value compactly.
  function fmtVal(v) {
    if (!isFinite(v)) return "";
    if (Math.abs(v) >= 100) return v.toFixed(0);
    if (Math.abs(v) >= 1)   return v.toFixed(1);
    return v.toFixed(2);
  }

  // Range-feature rows that currently have data (order matches the annotations panel).
  function activeFeatRows() {
    return FEAT_ROWS.filter(function (name) {
      return rangeFeatures.some(function (f) { return f.source === name; });
    });
  }

  // Compute panel boundary + height positions (all in CSS px, top-down).
  //
  // The line, annotations, and isoform panes are content-driven with minimum heights;
  // the legend band, score row, and sequence strip stay fixed. If the canvas is taller
  // than the sum of natural heights, the line panel absorbs the extra space so the
  // canvas is never under-filled.
  function getPanelLayout(cssH) {
    var lineH = Math.max(MIN_LINE_H, lineTracks.length * LINE_TRACK_H);
    var featH = Math.max(MIN_FEAT_H, activeFeatRows().length * FEAT_ROW_H);
    var isoH  = isoIsoforms.length > 1 ? isoIsoforms.length * ISO_ROW_H + ISO_PAD : 0;

    // 5 gaps: line→legend, legend→feat, feat→iso (if present), iso/feat→score, score→seq
    var nGaps = isoH > 0 ? 5 : 4;
    // each dynamic pane (line/feat/iso) reserves PANE_LABEL_H above it for its title
    var nLabels = 2 + (isoH > 0 ? 1 : 0);
    var naturalTotal = TOP_GUTTER + PANE_LABEL_H * nLabels + lineH + LEGEND_H + featH + isoH
      + SCORE_H + SEQ_H + PANEL_GAP * nGaps;
    var slack = Math.max(cssH - naturalTotal, 0);
    lineH += slack;  // give any leftover canvas height to the line panel

    var lineLabelTop = TOP_GUTTER;
    var lineTop   = lineLabelTop + PANE_LABEL_H;
    var legendTop = lineTop + lineH + PANEL_GAP;
    var featLabelTop = legendTop + LEGEND_H + PANEL_GAP;
    var featTop   = featLabelTop + PANE_LABEL_H;
    var isoLabelTop = featTop + featH + PANEL_GAP;
    var isoTop    = isoH > 0 ? isoLabelTop + PANE_LABEL_H : isoLabelTop;
    var scoreTop  = isoTop + isoH + (isoH > 0 ? PANEL_GAP : 0);
    var seqTop    = scoreTop + SCORE_H + PANEL_GAP;
    return {
      lineLabelTop: lineLabelTop,
      lineTop:   lineTop,
      lineH:     lineH,
      legendTop: legendTop,
      legendH:   LEGEND_H,
      featLabelTop: featLabelTop,
      featTop:   featTop,
      featH:     featH,
      isoLabelTop: isoLabelTop,
      isoTop:    isoTop,
      isoH:      isoH,
      scoreTop:  scoreTop,
      scoreH:    SCORE_H,
      seqTop:    seqTop,
      seqH:      SEQ_H,
      naturalTotal: naturalTotal,
    };
  }

  /* ── Core render ───────────────────────────────────────────────────────────── */

  function render() {
    // set the canvas's CSS height to fit its content-driven panes before measuring —
    // getPanelLayout(0) reports the natural total height with zero slack.
    var canvasEl = document.getElementById("ts_plot_div");
    if (canvasEl) {
      var naturalH = getPanelLayout(0).naturalTotal;
      canvasEl.style.height = naturalH + "px";
    }

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
    drawIsoformPane(ctx, inf, layout);
    drawScoreHeatmap(ctx, inf, layout);
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

  /* ── Line panel (row 1: continuous per-position scores) ──────────────────── */

  function drawLinePanel(ctx, inf, layout) {
    var top = layout.lineTop, h = layout.lineH;

    // panel background + border
    ctx.fillStyle   = "#fff";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);

    // y-axis label "Value" (vertical, in gutter)
    ctx.save();
    ctx.translate(11, top + h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.font         = "11px sans-serif";
    ctx.textAlign    = "center";
    ctx.textBaseline = "middle";
    ctx.fillStyle    = "#666";
    ctx.fillText("Value", 0, 0);
    ctx.restore();

    // pane title, above the panel
    ctx.save();
    ctx.font         = "bold 13px sans-serif";
    ctx.fillStyle    = "#555";
    ctx.textAlign    = "left";
    ctx.textBaseline = "top";
    ctx.fillText("Conservation and structure", LEFT_GUTTER, layout.lineLabelTop);
    ctx.restore();

    if (lineTracks.length === 0) return;

    // y-axis is pinned to 0 at the bottom (scores are normalized to [0, 1])
    var yMin = 0;
    var yMax = plotYMax;
    if (yMin >= yMax) yMax = yMin + 1;
    var ySpan = yMax - yMin;

    // fixed tick set at quarter intervals, filtered to what's on-screen
    var yTicks = [0, 0.25, 0.5, 0.75, 1.0].filter(function (v) { return v <= yMax + 1e-9; });
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

  var FEAT_ROWS = ["Phobius", "Pfam", "modification", "hydrophobic_patch", "UniProt", "UniProt_site"];
  var FEAT_ROW_LABELS = {hydrophobic_patch: "Hydro. patch", UniProt: "UniProt", UniProt_site: "UniProt site"};

  function drawFeaturePanel(ctx, inf, layout) {
    var top = layout.featTop, h = layout.featH;

    // panel background + border
    ctx.fillStyle   = "#fafafa";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_GUTTER, top, inf.dataW, h);

    // pane title, above the panel
    ctx.save();
    ctx.font         = "bold 13px sans-serif";
    ctx.fillStyle    = "#555";
    ctx.textAlign    = "left";
    ctx.textBaseline = "top";
    ctx.fillText("Annotations", LEFT_GUTTER, layout.featLabelTop);
    ctx.restore();

    // only show rows that actually have data
    var activeRows = activeFeatRows();
    if (activeRows.length === 0) return;

    // assign evenly-spaced y positions within the panel
    var n = activeRows.length;
    var dynYPos = {};
    activeRows.forEach(function (name, idx) { dynYPos[name] = idx + 1; });
    // small half-row margin above/below the outermost rows instead of a full row's worth
    var dynYMin = 0.5, dynYMax = n + 0.5;
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
      var label = FEAT_ROW_LABELS[name] || (name.charAt(0).toUpperCase() + name.slice(1));
      ctx.fillText(label, LEFT_GUTTER - 3, py);
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

  /* ── Isoform pane (one row per isoform, longest first) ───────────────────────── */

  function drawIsoformPane(ctx, inf, layout) {
    isoHitBoxes = [];
    if (layout.isoH <= 0 || isoIsoforms.length === 0) return;

    var top = layout.isoTop, h = layout.isoH;

    ctx.fillStyle   = "#fafafa";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_GUTTER, top, inf.dataW, h);

    // pane title, above the panel
    ctx.save();
    ctx.font         = "bold 13px sans-serif";
    ctx.fillStyle    = "#555";
    ctx.textAlign    = "left";
    ctx.textBaseline = "top";
    ctx.fillText("Isoforms", LEFT_GUTTER, layout.isoLabelTop);
    ctx.restore();

    var n = isoIsoforms.length;
    var innerH = h - ISO_PAD;
    var rowH   = innerH / n;
    var r      = getXRange();

    // grey exon-boundary lines, drawn behind the isoform rows (reserved for future
    // genewise integration — exonBoundaries is empty until that work lands)
    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.clip();
    ctx.strokeStyle = "#ddd";
    ctx.lineWidth   = 1;
    isoExonBoundaries.forEach(function (pos) {
      if (pos < r[0] || pos > r[1]) return;
      var x = posToX(pos, inf);
      ctx.beginPath();
      ctx.moveTo(x, top + ISO_PAD / 2);
      ctx.lineTo(x, top + h - ISO_PAD / 2);
      ctx.stroke();
    });
    ctx.restore();

    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.clip();

    isoIsoforms.forEach(function (iso, idx) {
      var rowTop = top + ISO_PAD / 2 + idx * rowH;
      var rowMid = rowTop + rowH / 2;
      var barH   = Math.max(4, rowH * 0.55);

      // row label in the left gutter
      ctx.fillStyle    = iso.isQuery ? "#333" : "#888";
      ctx.font         = (iso.isQuery ? "bold " : "") + "9px sans-serif";
      ctx.textAlign    = "right";
      ctx.textBaseline = "middle";
      var label = iso.label.length > 12 ? iso.label.slice(0, 11) + "…" : iso.label;
      ctx.fillText(label, LEFT_GUTTER - 3, rowMid);

      // centerline under internal gaps only — i.e. "skipped" spans that fall strictly
      // between two present spans, not before the first or after the last (no line
      // is drawn for missing N-/C-terminal sequence, only for an internal deletion)
      var presentSpans = iso.present || [];
      if (presentSpans.length) {
        var firstPos = presentSpans[0][0];
        var lastPos  = presentSpans[presentSpans.length - 1][1];
        ctx.strokeStyle = "#ccc";
        ctx.lineWidth   = 1;
        (iso.skipped || []).forEach(function (span) {
          var start = span[0], stop = span[1];
          if (start <= firstPos || stop >= lastPos) return;  // N-/C-terminal, skip
          if (stop < r[0] || start > r[1]) return;
          var x0 = posToX(start - 0.5, inf);
          var x1 = posToX(stop + 0.5, inf);
          ctx.beginPath();
          ctx.moveTo(x0, rowMid);
          ctx.lineTo(x1, rowMid);
          ctx.stroke();
        });
      }

      // present spans: grey by default, colored only where they overlap a projected
      // Pfam domain (sub-segmented so the colored region matches the domain exactly)
      var GREY = "#bbbbbb";
      presentSpans.forEach(function (span) {
        var start = span[0], stop = span[1];
        if (stop < r[0] || start > r[1]) return;

        var domains = (iso.domains || [])
          .filter(function (d) { return d.stop >= start && d.start <= stop; })
          .map(function (d) { return {start: Math.max(d.start, start), stop: Math.min(d.stop, stop),
                                      color: d.color}; })
          .sort(function (a, b) { return a.start - b.start; });

        var segs = [];
        var cursor = start;
        domains.forEach(function (d) {
          if (d.start > cursor) segs.push({start: cursor, stop: d.start - 1, color: GREY});
          segs.push({start: d.start, stop: d.stop, color: d.color});
          cursor = d.stop + 1;
        });
        if (cursor <= stop) segs.push({start: cursor, stop: stop, color: GREY});

        segs.forEach(function (seg) {
          var x0 = posToX(seg.start - 0.5, inf);
          var x1 = posToX(seg.stop + 0.5, inf);
          if (x1 <= x0 + 0.5) return;
          ctx.globalAlpha = iso.isQuery ? 1.0 : 0.45;
          ctx.fillStyle    = seg.color;
          ctx.fillRect(x0, rowTop + (rowH - barH) / 2, x1 - x0, barH);
          ctx.globalAlpha = 1;
        });
      });

      // insertion markers: downward-pointing caret above the bar, with a matching
      // line on the bar itself marking the exact insertion site, hoverable
      var barTop    = rowTop + (rowH - barH) / 2;
      var barBottom = barTop + barH;
      (iso.inserts || []).forEach(function (ins) {
        if (ins.after < r[0] || ins.after > r[1]) return;
        var cx = posToX(ins.after + 0.5, inf);
        var triW = 9, triH = barH * 0.7 * 1.5;
        var triBottom = barTop + barH * 0.7, triTop = triBottom - triH;
        ctx.fillStyle = "#d84315";
        ctx.beginPath();
        ctx.moveTo(cx - triW / 2, triTop);
        ctx.lineTo(cx + triW / 2, triTop);
        ctx.lineTo(cx, triBottom);
        ctx.closePath();
        ctx.fill();

        ctx.strokeStyle = "#d84315";
        ctx.lineWidth   = 1.5;
        ctx.beginPath();
        ctx.moveTo(cx, barTop);
        ctx.lineTo(cx, barBottom);
        ctx.stroke();

        isoHitBoxes.push({isoIdx: idx, insert: ins, x0: cx - 4, x1: cx + 4,
                          y0: triTop, y1: barBottom});
      });
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

  /* ── Tag-site score heatmap (white→green, sits directly above the seq strip) ── */

  function drawScoreHeatmap(ctx, inf, layout) {
    var markTop = layout.scoreTop, markH = SCORE_MARK_H;
    var top = markTop + markH, h = layout.scoreH - markH;

    ctx.fillStyle   = "#f8f8f8";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_GUTTER, top, inf.dataW, h);

    // row label in the left gutter, matching the sequence strip's style
    ctx.fillStyle    = "#777";
    ctx.font         = "10px sans-serif";
    ctx.textAlign    = "right";
    ctx.textBaseline = "middle";
    ctx.fillText("Score", LEFT_GUTTER - 6, top + h / 2);

    if (scoreTrack.length === 0) return;

    var r     = getXRange();
    var lo    = Math.max(1, Math.floor(r[0]));
    var hi    = Math.min(scoreTrack.length, Math.ceil(r[1]));

    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.clip();

    for (var pos = lo; pos <= hi; pos++) {
      var s = scoreTrack[pos - 1];
      if (s === null || s === undefined) continue;
      var x0 = posToX(pos - 0.5, inf);
      var x1 = posToX(pos + 0.5, inf);
      var w  = x1 - x0;
      if (w < 0.5) continue;
      ctx.fillStyle = scoreMasked[pos - 1] ? MASKED_HEX : scoreCmapHex(s / scoreMax);
      ctx.fillRect(x0, top + 1, Math.max(w - 0.5, 0.5), h - 2);
    }

    ctx.restore();

    // downward-facing green triangle above each suggested site, in the reserved top margin
    if (suggestedSites.length) {
      ctx.save();
      ctx.beginPath();
      ctx.rect(LEFT_GUTTER, markTop, inf.dataW, markH);
      ctx.clip();
      ctx.fillStyle = "#28a745";
      for (var i = 0; i < suggestedSites.length; i++) {
        var sp = suggestedSites[i];
        if (sp < lo || sp > hi) continue;
        var cx = posToX(sp, inf);
        var triW = 8, triH = markH - 1;
        ctx.beginPath();
        ctx.moveTo(cx - triW / 2, markTop);
        ctx.lineTo(cx + triW / 2, markTop);
        ctx.lineTo(cx, markTop + triH);
        ctx.closePath();
        ctx.fill();
      }
      ctx.restore();
    }
  }

  // viridis colormap stops — must match utils.results._VIRIDIS so the heatmap
  // and structure coloring (utils.scoring.score_cmap_hex) agree
  var _SCORE_VIRIDIS = [
    [0.000, 68,  1,  84], [0.111, 72,  40, 120], [0.222, 62,  74, 137],
    [0.333, 49, 104, 142], [0.444, 38, 130, 142], [0.556, 31, 158, 137],
    [0.667, 53, 183, 121], [0.778, 110, 206,  88], [0.889, 181, 222,  43],
    [1.000, 253, 231,  37],
  ];

  function scoreCmapHex(frac) {
    frac = Math.max(0, Math.min(1, frac));
    var stops = _SCORE_VIRIDIS;
    for (var i = 0; i < stops.length - 1; i++) {
      var t0 = stops[i][0], t1 = stops[i + 1][0];
      if (frac <= t1) {
        var f = t1 > t0 ? (frac - t0) / (t1 - t0) : 0;
        var r = Math.round(stops[i][1] + f * (stops[i + 1][1] - stops[i][1]));
        var g = Math.round(stops[i][2] + f * (stops[i + 1][2] - stops[i][2]));
        var b = Math.round(stops[i][3] + f * (stops[i + 1][3] - stops[i][3]));
        return hexFromRgb(r, g, b);
      }
    }
    var last = stops[stops.length - 1];
    return hexFromRgb(last[1], last[2], last[3]);
  }

  function hexFromRgb(r, g, b) {
    function hx(v) { var s = v.toString(16); return s.length < 2 ? "0" + s : s; }
    return "#" + hx(r) + hx(g) + hx(b);
  }

  /* ── Sequence strip (row 3) ────────────────────────────────────────────────── */

  function drawSeqStrip(ctx, inf, layout) {
    var top      = layout.seqTop, h = layout.seqH;
    var TICK_H   = SEQ_TICK_H;
    var LETTER_H = h - TICK_H;
    var LETTER_THRESHOLD = 9;

    // strip background + border + separator between letters and ticks
    ctx.fillStyle   = "#f8f8f8";
    ctx.fillRect(LEFT_GUTTER, top, inf.dataW, h);
    ctx.strokeStyle = "#d8d8d8";
    ctx.lineWidth   = 1;
    ctx.beginPath();
    ctx.moveTo(LEFT_GUTTER, top + LETTER_H);
    ctx.lineTo(LEFT_GUTTER + inf.dataW, top + LETTER_H);
    ctx.stroke();

    // row label in the left gutter, matching the score row's style
    ctx.fillStyle    = "#777";
    ctx.font         = "10px sans-serif";
    ctx.textAlign    = "right";
    ctx.textBaseline = "middle";
    ctx.fillText("Sequence", LEFT_GUTTER - 6, top + LETTER_H / 2);

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

      var isCommitted = committedSet.has(pos);
      var baseColor   = CLUSTALX_COLORS[aa] || "#555555";

      if (cellW >= LETTER_THRESHOLD) {
        // letter mode
        var bgColor = isCommitted ? "#28a745" : "#ffffff";
        ctx.fillStyle = bgColor;
        ctx.fillRect(x0, top, w, LETTER_H);
        if (!isCommitted) {
          ctx.strokeStyle = "#e0e0e0";
          ctx.lineWidth   = 0.5;
          ctx.strokeRect(x0, top, w, LETTER_H);
        }
        var letterColor = isCommitted ? "#ffffff" : baseColor;
        var fontSize    = Math.min(14, Math.max(8, w * 0.75));
        ctx.fillStyle    = letterColor;
        ctx.font         = "bold " + fontSize + "px 'Courier New', monospace";
        ctx.textAlign    = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(aa, cx, top + LETTER_H / 2);
      } else {
        // box mode
        var boxColor = isCommitted ? "#28a745" : baseColor;
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

    // show isoform info when hovering the isoform pane
    if (layout.isoH > 0 && cy >= layout.isoTop && cy <= layout.isoTop + layout.isoH) {
      // insertion markers take priority — they're a point feature, not a span
      var hitInsert = isoHitBoxes.find(function (box) {
        return cx >= box.x0 && cx <= box.x1 && cy >= box.y0 && cy <= box.y1;
      });
      if (hitInsert) {
        var isoLabel = isoLabelWithAcc(isoIsoforms[hitInsert.isoIdx]);
        lines = ["+" + hitInsert.insert.length + " aa insert after " + hitInsert.insert.after +
                 " — " + isoLabel];
      } else {
        var n = isoIsoforms.length;
        var innerH = layout.isoH - ISO_PAD;
        var rowH   = innerH / n;
        var rowIdx = Math.floor((cy - layout.isoTop - ISO_PAD / 2) / rowH);
        var iso    = isoIsoforms[rowIdx];
        if (iso) {
          var isoLbl = isoLabelWithAcc(iso);
          var inPresent = (iso.present || []).some(function (s) { return pos >= s[0] && pos <= s[1]; });
          if (inPresent) {
            var domain = (iso.domains || []).find(function (d) { return pos >= d.start && pos <= d.stop; });
            lines = [isoLbl + (domain ? " — Pfam: " + domain.desc + " (" + domain.start +
                     "–" + domain.stop + ")" : " — present")];
          } else {
            var skip = (iso.skipped || []).find(function (s) { return pos >= s[0] && pos <= s[1]; });
            if (skip) lines = [isoLbl + " — absent (" + skip[0] + "–" + skip[1] + ")"];
          }
        }
      }
    }

    // show tag-site score + failed-criterion flags when hovering the score row
    if (cy >= layout.scoreTop && cy <= layout.scoreTop + layout.scoreH) {
      var s = scoreTrack[pos - 1];
      if (s !== null && s !== undefined) {
        var aa = seqArray[pos - 1] || "";
        if (scoreMasked[pos - 1]) {
          var reasons = scoreMaskReasons[pos - 1] || [];
          lines = [aa + pos + ": masked" + (reasons.length ? " [" + reasons.join(", ") + "]" : "")];
        } else {
          var flags = scoreFlags[pos - 1] || [];
          lines = [aa + pos + ": " + s + (flags.length ? " [" + flags.join(", ") + "]" : "")];
        }
      }
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

      // isoform pane click → no residue selection (rows aren't tag sites)
      if (layout.isoH > 0 && cy >= layout.isoTop && cy <= layout.isoTop + layout.isoH) return;

      if (!inf || !inDataArea(cx, inf)) return;
      shinySetPos(inputNames.residue_click, xToPos(cx, inf));
    });

    // ── Dblclick: reset zoom to full range (legend dblclick is a no-op) ─────────

    canvas.addEventListener("dblclick", function (e) {
      if (wasDrag) { wasDrag = false; return; }
      var rect   = canvas.getBoundingClientRect();
      var cy     = e.clientY - rect.top;
      var inf    = getCanvasInfo();
      var layout = getPanelLayout(inf ? inf.cssH : 520);

      if (cy >= layout.legendTop && cy <= layout.legendTop + LEGEND_H) return;

      currentRange = null;
      render();
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

    viewer.setClickable({}, true, function (atom) {
      if (!atom || atom.resi == null) return;
      shinySetPos(inputNames.struct_click, parseInt(atom.resi));
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

    committedSet.forEach(function (pos) {
      viewer.setStyle({ resi: pos }, { cartoon: { color: "#28a745" } });
      viewer.addStyle({ resi: pos }, { stick: { colorscheme: "amino", radius: 0.2 } });
      viewer.addStyle({ resi: pos }, { sphere: { colorscheme: "amino", scale: 0.4 } });
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
    var runTitleEl = document.getElementById("ts-run-title");
    if (runTitleEl) runTitleEl.textContent = plotTitle;
    plotYMax      = 1.0;
    lineTracks    = msg.lineTracks    || [];
    // sasa and hydrophobicity lines clutter the default view — start them hidden
    var defaultHidden = lineTracks
      .map(function (tr) { return tr.name; })
      .filter(function (n) { return /_sasa$/.test(n) || /_hydro$/.test(n); });
    rangeFeatures = msg.rangeFeatures || [];
    seqArray      = (msg.seq || "").split("");
    scoreTrack    = msg.scoreTrack    || [];
    scoreFlags    = msg.scoreFlags    || [];
    scoreMasked   = msg.scoreMasked   || [];
    scoreMaskReasons = msg.scoreMaskReasons || [];
    scoreMax      = (typeof msg.scoreMax === "number" && msg.scoreMax > 0) ? msg.scoreMax : 1;
    suggestedSites = msg.suggestedSites || [];
    var isoformPane   = msg.isoformPane || {};
    isoIsoforms       = isoformPane.isoforms || [];
    isoExonBoundaries = isoformPane.exonBoundaries || [];
    hiddenTracks  = new Set(defaultHidden);
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

  Shiny.addCustomMessageHandler("tagsites_clear_struct", function (msg) {
    if (viewer) {
      viewer.clear();
      viewer = null;
    }
    var container = document.getElementById("ts-viewer-container");
    if (container) container.innerHTML = "";
    var pane = document.getElementById("ts-struct-pane");
    if (pane) pane.style.display = "none";
  });

  Shiny.addCustomMessageHandler("tagsites_set_states", function (msg) {
    committedSet = new Set(msg.committed);
    render();
    applyViewerColors();
  });

  // Incremental rescore: score row + isoform pane only. Unlike tagsites_set_plot (the
  // load handler), this preserves zoom, committed sites, hidden tracks, and drag mode —
  // needed so toggling an Advanced/Rescoring control doesn't wipe the user's view.
  Shiny.addCustomMessageHandler("tagsites_set_scores", function (msg) {
    scoreTrack       = msg.scoreTrack       || [];
    scoreFlags       = msg.scoreFlags       || [];
    scoreMasked      = msg.scoreMasked      || [];
    scoreMaskReasons = msg.scoreMaskReasons || [];
    scoreMax         = (typeof msg.scoreMax === "number" && msg.scoreMax > 0) ? msg.scoreMax : 1;
    suggestedSites   = msg.suggestedSites   || [];
    if (msg.isoformPane) {
      isoIsoforms       = msg.isoformPane.isoforms       || [];
      isoExonBoundaries = msg.isoformPane.exonBoundaries || [];
    }
    render();
  });

  Shiny.addCustomMessageHandler("tagsites_set_colors", function (msg) {
    perResidueColors = msg.colors || [];
    applyViewerColors();
    // when no data colors are set (None), show the default N→C rainbow legend
    var legend = msg.legend || (perResidueColors.length === 0 ? {type: "rainbow"} : null);
    renderStructLegend(legend);
  });

  // ── Structure color legend ─────────────────────────────────────────────────────

  // gradient CSS strings — must match the Python colormaps in utils/results.py
  var _CMAP_CSS = {
    viridis: "linear-gradient(to right," +
      "#440154,#48306e,#3e4989,#31688e,#26838f,#1f9e89,#35b779,#6ece58,#b5de2b,#fde725)",
    plasma:  "linear-gradient(to right,#0d0887,#7e03a8,#cb4679,#f89540,#f0f921)",
    cool:    "linear-gradient(to right,#00ffff,#ff00ff)",
    // pLDDT gradient: stops at band boundaries (0%, 50%, 70%, 100%)
    plddt:   "linear-gradient(to right,#ff7d45 0%,#ffdb13 50%,#65cbf3 70%,#0053d6 100%)",
    // blue-white-red diverging: polar (buried/hydrophilic) → hydrophobic surface exposure
    bwr:     "linear-gradient(to right,#2166ac,#f7f7f7,#b2182b)",
  };

  // rainbow N→C matches 3Dmol's default chain-spectrum coloring
  var _RAINBOW_CSS = "linear-gradient(to right,#0000ff,#00ffff,#00ff00,#ffff00,#ff0000)";

  function renderStructLegend(legend) {
    var el = document.getElementById("ts-struct-legend");
    if (!el) return;
    if (!legend) { el.innerHTML = ""; return; }

    if (legend.type === "gradient") {
      var gradCSS = _CMAP_CSS[legend.cmap] || _CMAP_CSS.viridis;
      el.innerHTML =
        '<div class="ts-sleg-label">' + legend.label + '</div>' +
        '<div class="ts-sleg-bar-row">' +
          '<span class="ts-sleg-tick">' + legend.vmin + '</span>' +
          '<div class="ts-sleg-bar" style="background:' + gradCSS + '"></div>' +
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
  // Show a spinner next to the file input when a .zip bundle is selected.
  document.addEventListener("change", function (e) {
    var inp = e.target;
    if (inp.type !== "file" || !inp.files || !inp.files.length) return;
    var cardBody = inp.closest(".card-body");
    if (!cardBody) return;
    var spinner = cardBody.querySelector(".ts-bundle-spinner");
    if (!spinner) return;
    spinner.style.display = inp.files[0].name.toLowerCase().endsWith(".zip") ? "flex" : "none";
  });

  Shiny.addCustomMessageHandler("tagsites_bundle_done", function (msg) {
    document.querySelectorAll(".ts-bundle-spinner").forEach(function (el) {
      el.style.display = "none";
    });
  });

  // Disable a "Fetch"/"Search" button and show its spinner while the (blocking,
  // synchronous) server-side call is in flight, so repeat clicks can't queue up.
  // doneMessage is sent by the server once the call finishes (success or failure).
  function wireBlockingButton(btnClass, wrapClass, spinnerClass, doneMessage) {
    document.addEventListener("click", function (e) {
      var btn = e.target.closest(btnClass);
      if (!btn) return;
      var container = btn.closest(wrapClass);
      btn.disabled = true;
      var spinner = container && container.querySelector(spinnerClass);
      if (spinner) spinner.style.display = "flex";
    });

    Shiny.addCustomMessageHandler(doneMessage, function (msg) {
      document.querySelectorAll(btnClass).forEach(function (btn) {
        btn.disabled = false;
      });
      document.querySelectorAll(spinnerClass).forEach(function (el) {
        el.style.display = "none";
      });
    });
  }

  wireBlockingButton(".ts-genomic-fetch-btn", ".ts-genomic-fetch-wrap",
                      ".ts-genomic-spinner", "tagsites_genomic_fetch_done");
  wireBlockingButton(".ts-uniprot-search-btn", ".ts-uniprot-search-wrap",
                      ".ts-uniprot-spinner", "tagsites_uniprot_search_done");

  Shiny.addCustomMessageHandler("tagsites_trigger_download", function (msg) {
    var link = document.getElementById("progress-download_results");
    if (link) link.click();
  });

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
