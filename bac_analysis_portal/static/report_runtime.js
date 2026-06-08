function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function getTaskMethod(task) {
  return String(task?.params?.method || task?.method || "").trim();
}

function hasNonEmptyObject(value) {
  return Boolean(
    value
    && typeof value === "object"
    && !Array.isArray(value)
    && Object.keys(value).length
  );
}

function hasRenderableCgviewPayload(payload) {
  const cgview = payload?.cgview;
  if (!cgview || typeof cgview !== "object") return false;
  const features = Array.isArray(cgview.features) ? cgview.features : [];
  const tracks = Array.isArray(cgview.tracks) ? cgview.tracks : [];
  const contigs = Array.isArray(cgview.sequence?.contigs) ? cgview.sequence.contigs : [];
  const legendItems = Array.isArray(cgview.legend?.items) ? cgview.legend.items : [];
  return features.length > 0 && tracks.length > 0 && contigs.length > 0 && legendItems.length > 0;
}

function normalizeCgviewLookupKey(value) {
  return String(value || "")
    .trim()
    .toLowerCase()
    .replace(/\.gbk$/i, "")
    .replace(/\s+/g, "")
    .replace(/[|/\\]+/g, "_");
}

function parseCgviewCoordinate(value) {
  const numeric = Number(String(value ?? "").trim());
  return Number.isFinite(numeric) && numeric > 0 ? Math.round(numeric) : null;
}

function parseCgviewStrand(value) {
  const normalized = String(value || "").trim().toLowerCase();
  if (!normalized) return 1;
  if (normalized === "-" || normalized === "-1" || normalized.includes("minus") || normalized.includes("reverse")) {
    return -1;
  }
  return 1;
}

function buildCgviewContigAliasMap(payload) {
  const contigs = Array.isArray(payload?.cgview?.sequence?.contigs) ? payload.cgview.sequence.contigs : [];
  const aliasMap = new Map();
  contigs.forEach((contig) => {
    if (!contig || typeof contig !== "object") return;
    const primary = String(contig.name || contig.id || contig.seqID || contig.accession || "").trim();
    if (!primary) return;
    [
      contig.name,
      contig.id,
      contig.seqID,
      contig.accession,
      contig.label,
      primary,
    ].forEach((candidate) => {
      const key = normalizeCgviewLookupKey(candidate);
      if (key) aliasMap.set(key, primary);
    });
  });
  return aliasMap;
}

function getCgviewColumnValue(row, columns, candidates) {
  if (!Array.isArray(row) || !Array.isArray(columns)) return "";
  return getFirstMatchingCell(row, columns, candidates);
}

function appendCgviewOverlayFeatures(payload, currentMap, sections) {
  const cgview = payload?.cgview;
  if (!cgview || typeof cgview !== "object") return payload;
  const aliasMap = buildCgviewContigAliasMap(payload);
  if (!aliasMap.size) return payload;

  const mapKey = normalizeCgviewLookupKey(currentMap?.key);
  const isMainMap = currentMap?.role === "main" || mapKey === "main";
  const currentMapAliases = new Set([
    mapKey,
    normalizeCgviewLookupKey(currentMap?.asset_name),
    normalizeCgviewLookupKey(currentMap?.label),
  ].filter(Boolean));
  const overlayFeatures = [];
  const dedupe = new Set();

  const maybePushFeature = ({
    row,
    columns,
    source,
    legend,
    labelPrefix = "",
    nameCandidates,
    contigCandidates,
    startCandidates,
    endCandidates,
    strandCandidates,
    noteCandidates,
    favorite = false,
  }) => {
    const rawContig = getCgviewColumnValue(row, columns, contigCandidates);
    const resolvedContig = aliasMap.get(normalizeCgviewLookupKey(rawContig)) || "";
    if (!resolvedContig) return;
    if (!isMainMap && !currentMapAliases.has(normalizeCgviewLookupKey(rawContig)) && !currentMapAliases.has(normalizeCgviewLookupKey(resolvedContig))) {
      return;
    }
    const start = parseCgviewCoordinate(getCgviewColumnValue(row, columns, startCandidates));
    const stop = parseCgviewCoordinate(getCgviewColumnValue(row, columns, endCandidates));
    if (start === null || stop === null) return;
    const mapStart = Math.max(1, Math.min(start, stop));
    const mapStop = Math.max(start, stop);
    const geneName = getCgviewColumnValue(row, columns, nameCandidates) || legend;
    const detail = getCgviewColumnValue(row, columns, noteCandidates);
    const displayName = labelPrefix ? `${labelPrefix} ${geneName}` : geneName;
    const dedupeKey = [
      source,
      resolvedContig,
      mapStart,
      mapStop,
      displayName,
    ].join("|");
    if (dedupe.has(dedupeKey)) return;
    dedupe.add(dedupeKey);
    overlayFeatures.push({
      name: displayName,
      type: "misc_feature",
      source,
      contig: resolvedContig,
      legend,
      start: mapStart,
      stop: mapStop,
      strand: parseCgviewStrand(getCgviewColumnValue(row, columns, strandCandidates)),
      tags: [source, legend],
      score: 1,
      favorite,
      visible: true,
      ...(detail ? { tags: [source, legend, detail] } : {}),
    });
  };

  const virulenceColumns = sections?.resistance_virulence?.virulence_elements?.columns || [];
  const virulenceRows = sections?.resistance_virulence?.virulence_elements?.rows || [];
  virulenceRows.forEach((row) => {
    maybePushFeature({
      row,
      columns: virulenceColumns,
      source: "portal-virulence",
      legend: "毒力基因",
      labelPrefix: "VF",
      nameCandidates: ["基因名称", "VF名称"],
      contigCandidates: ["Contig名称", "所在序列", "sequence_id"],
      startCandidates: ["起始碱基", "基因起始", "gene_start", "start"],
      endCandidates: ["终止碱基", "基因终止", "gene_end", "end"],
      strandCandidates: ["正负链", "strand"],
      noteCandidates: ["VF分类", "产物", "product"],
      favorite: true,
    });
  });

  const resistanceColumns = sections?.resistance_virulence?.resistance_elements?.columns || [];
  const resistanceRows = sections?.resistance_virulence?.resistance_elements?.rows || [];
  resistanceRows.forEach((row) => {
    maybePushFeature({
      row,
      columns: resistanceColumns,
      source: "portal-resistance",
      legend: "耐药基因",
      labelPrefix: "ARG",
      nameCandidates: ["基因名称"],
      contigCandidates: ["Contig名称", "所在序列", "sequence_id"],
      startCandidates: ["起始碱基", "基因起始", "gene_start", "start"],
      endCandidates: ["终止碱基", "基因终止", "gene_end", "end"],
      strandCandidates: ["正负链", "strand"],
      noteCandidates: ["耐药药物", "产物", "product"],
      favorite: true,
    });
  });

  const mgeColumns = sections?.mge_monitoring?.elements?.columns || [];
  const mgeRows = sections?.mge_monitoring?.elements?.rows || [];
  mgeRows.forEach((row) => {
    maybePushFeature({
      row,
      columns: mgeColumns,
      source: "portal-mge",
      legend: "移动元件",
      nameCandidates: ["元件类型", "mge_type", "注释", "annotation"],
      contigCandidates: ["所在序列", "sequence_id", "Contig名称"],
      startCandidates: ["起始", "gene_start", "start"],
      endCandidates: ["终止", "gene_end", "end"],
      strandCandidates: ["strand"],
      noteCandidates: ["注释", "识别方法", "annotation", "method"],
      favorite: true,
    });
  });

  if (!overlayFeatures.length) return payload;
  cgview.features = [...(Array.isArray(cgview.features) ? cgview.features : []), ...overlayFeatures];
  return payload;
}

function getCgviewOverlayTrackSources(viewer) {
  if (!viewer?.tracks) return [];
  const toList = (value) => {
    if (Array.isArray(value)) return value;
    if (value && typeof value !== "string" && typeof value.length === "number") return Array.from(value);
    return [value];
  };
  return viewer.tracks()
    .filter((track) => track?.dataMethod === "source")
    .map((track) => toList(track?.dataKeys))
    .flat()
    .map((value) => String(value || ""))
    .filter((value) => ["portal-resistance", "portal-virulence", "portal-mge"].includes(value));
}

function applyCgviewOverlayVisibility(viewer, activeSources) {
  if (!viewer?.tracks) return;
  const trackSources = new Set(["portal-resistance", "portal-virulence", "portal-mge"]);
  const toList = (value) => {
    if (Array.isArray(value)) return value;
    if (value && typeof value !== "string" && typeof value.length === "number") return Array.from(value);
    return [value];
  };
  viewer.tracks()
    .filter((track) => track?.dataMethod === "source")
    .each((_, track) => {
      const keys = toList(track?.dataKeys);
      const overlayKeys = keys.map((value) => String(value || "")).filter((value) => trackSources.has(value));
      if (!overlayKeys.length) return;
      viewer.updateTracks(track, { visible: overlayKeys.some((value) => activeSources.has(value)) });
    });
  viewer.drawFull();
}

function normalizeCgviewPayload(payload, mapLabel) {
  const cgview = payload?.cgview;
  if (!cgview || typeof cgview !== "object") return payload;

  const legendPalette = {
    CDS: { swatchColor: "rgba(30, 98, 106, 1)", decoration: "arrow" },
    "毒力基因": { swatchColor: "rgba(187, 68, 90, 1)", decoration: "arc" },
    "耐药基因": { swatchColor: "rgba(209, 126, 56, 1)", decoration: "arc" },
    "移动元件": { swatchColor: "rgba(95, 92, 153, 1)", decoration: "arc" },
    "GC Content": { swatchColor: "rgba(32, 58, 92, 1)", decoration: "arc" },
    "GC Skew+": { swatchColor: "rgba(20, 140, 72, 1)", decoration: "arc" },
    "GC Skew-": { swatchColor: "rgba(166, 32, 130, 1)", decoration: "arc" },
    tRNA: { swatchColor: "rgba(199, 129, 72, 1)", decoration: "arc" },
    rRNA: { swatchColor: "rgba(164, 83, 83, 1)", decoration: "arc" },
    ncRNA: { swatchColor: "rgba(116, 92, 156, 1)", decoration: "arc" },
    tmRNA: { swatchColor: "rgba(116, 92, 156, 1)", decoration: "arc" },
    repeat_region: { swatchColor: "rgba(108, 119, 132, 1)", decoration: "arc" },
    misc_feature: { swatchColor: "rgba(126, 138, 122, 1)", decoration: "arc" },
    mobile_element: { swatchColor: "rgba(143, 101, 66, 1)", decoration: "arc" },
  };
  const features = Array.isArray(cgview.features) ? cgview.features : [];
  const seenLegendNames = new Set();
  const legendItems = [];

  features.forEach((feature) => {
    if (!feature || typeof feature !== "object") return;
    const rawType = String(feature.legend || feature.type || "misc_feature").trim();
    const legendName = rawType || "misc_feature";
    feature.legend = legendName;
    const paletteEntry = legendPalette[legendName] || legendPalette.misc_feature;
    if (!seenLegendNames.has(legendName)) {
      seenLegendNames.add(legendName);
      legendItems.push({
        name: legendName,
        swatchColor: paletteEntry.swatchColor,
        decoration: paletteEntry.decoration,
      });
    }
  });

  ["GC Content", "GC Skew+", "GC Skew-"].forEach((legendName) => {
    if (seenLegendNames.has(legendName)) return;
    const paletteEntry = legendPalette[legendName];
    if (!paletteEntry) return;
    seenLegendNames.add(legendName);
    legendItems.push({
      name: legendName,
      swatchColor: paletteEntry.swatchColor,
      decoration: paletteEntry.decoration,
    });
  });

  cgview.name = mapLabel || cgview.name || "CGView 环形图";
  cgview.settings = {
    ...(cgview.settings || {}),
    format: "circular",
    backgroundColor: "rgba(247, 244, 238, 1)",
    showShading: false,
    arrowHeadLength: 0.45,
    initialMapThicknessProportion: 0.06,
    maxMapThicknessProportion: 0.115,
  };
  cgview.backbone = {
    ...(cgview.backbone || {}),
    color: "rgba(57, 76, 97, 0.9)",
    colorAlternate: "rgba(176, 186, 198, 0.95)",
    thickness: 12,
    decoration: "arrow",
  };
  cgview.ruler = {
    ...(cgview.ruler || {}),
    color: "rgba(86, 98, 110, 0.82)",
    font: "SansSerif, plain, 11",
    tickCount: 14,
    tickWidth: 1,
    tickLength: 5,
    rulerPadding: 14,
    spacing: 2,
  };
  cgview.annotation = {
    ...(cgview.annotation || {}),
    color: "rgba(63, 73, 84, 0.78)",
    font: "SansSerif, plain, 11",
    labelPlacement: "default",
    onlyDrawFavorites: false,
    labelLineLength: 16,
    priorityMax: 40,
  };
  cgview.dividers = {
    ...(cgview.dividers || {}),
    visible: true,
    color: "rgba(121, 133, 145, 0.18)",
    thickness: 1,
    spacing: 2,
  };
  cgview.legend = {
    ...(cgview.legend || {}),
    position: "top-right",
    on: "canvas",
    backgroundColor: "rgba(252, 250, 246, 0.96)",
    defaultFont: "SansSerif, plain, 12",
    defaultFontColor: "rgba(52, 62, 73, 1)",
    textAlignment: "left",
    defaultMinArcLength: 1,
    visible: true,
    items: legendItems,
  };
  const baseTracks = (Array.isArray(cgview.tracks) ? cgview.tracks : []).map((track, index) => {
    const dataKeys = Array.isArray(track?.dataKeys) ? track.dataKeys : [track?.dataKeys];
    const keyText = dataKeys.map((value) => String(value || "")).join(" ");
    if (track?.dataMethod === "sequence" && keyText.includes("gc-content")) {
      return {
        ...track,
        name: "GC Content",
        position: "inside",
        thicknessRatio: 0.4,
        dataOptions: {
          ...(track?.dataOptions || {}),
          window: 2400,
          step: 120,
          deviation: "average",
        },
      };
    }
    if (track?.dataMethod === "sequence" && keyText.includes("gc-skew")) {
      return {
        ...track,
        name: "GC Skew",
        position: "inside",
        thicknessRatio: 0.34,
        dataOptions: {
          ...(track?.dataOptions || {}),
          window: 2200,
          step: 110,
          deviation: "average",
        },
      };
    }
    return {
      ...track,
      separateFeaturesBy: track?.separateFeaturesBy || "strand",
      position: track?.position || "both",
      thicknessRatio: index === 0 ? 0.34 : (track?.thicknessRatio || 1),
    };
  });
  const overlayTrackConfigs = [
    { name: "毒力标记", source: "portal-virulence", position: "outside", thicknessRatio: 0.08 },
    { name: "耐药标记", source: "portal-resistance", position: "outside", thicknessRatio: 0.08 },
    { name: "移动元件标记", source: "portal-mge", position: "outside", thicknessRatio: 0.07 },
  ];
  overlayTrackConfigs.forEach((track) => {
    if (!features.some((feature) => feature?.source === track.source)) return;
    if (baseTracks.some((item) => item?.dataMethod === "source" && Array.isArray(item?.dataKeys) && item.dataKeys.includes(track.source))) {
      return;
    }
    baseTracks.push({
      name: track.name,
      dataType: "feature",
      dataMethod: "source",
      dataKeys: [track.source],
      separateFeaturesBy: "strand",
      position: track.position,
      thicknessRatio: track.thicknessRatio,
      favorite: true,
      visible: true,
    });
  });
  const hasGcContentTrack = baseTracks.some((track) => track?.dataMethod === "sequence" && String(track?.dataKeys || "").includes("gc-content"));
  const hasGcSkewTrack = baseTracks.some((track) => track?.dataMethod === "sequence" && String(track?.dataKeys || "").includes("gc-skew"));
  if (!hasGcContentTrack) {
    baseTracks.push({
      name: "GC Content",
      dataType: "plot",
      dataMethod: "sequence",
      dataKeys: ["gc-content"],
      dataOptions: {
        window: 2400,
        step: 120,
        deviation: "average",
      },
      position: "inside",
      thicknessRatio: 0.4,
      favorite: false,
      visible: true,
    });
  }
  if (!hasGcSkewTrack) {
    baseTracks.push({
      name: "GC Skew",
      dataType: "plot",
      dataMethod: "sequence",
      dataKeys: ["gc-skew"],
      dataOptions: {
        window: 2200,
        step: 110,
        deviation: "average",
      },
      position: "inside",
      thicknessRatio: 0.34,
      favorite: false,
      visible: true,
    });
  }
  cgview.tracks = baseTracks;

  return payload;
}

function isMetaReport(data) {
  if (getTaskMethod(data?.task) === "meta") return true;
  const binning = data?.sections?.binning_results || {};
  return Boolean(
    hasNonEmptyObject(binning?.quality?.summary)
    || hasNonEmptyObject(binning?.taxonomy?.summary)
    || hasNonEmptyObject(binning?.viral_assembly?.summary)
    || (Array.isArray(binning?.quality?.table?.rows) && binning.quality.table.rows.length)
    || (Array.isArray(binning?.taxonomy?.table?.rows) && binning.taxonomy.table.rows.length)
    || (Array.isArray(binning?.viral_assembly?.table?.rows) && binning.viral_assembly.table.rows.length)
  );
}

function formatDateTime(value) {
  if (!value) return "-";
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? String(value) : date.toLocaleString();
}

function formatBases(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "-";
  if (numeric >= 1e9) return `${(numeric / 1e9).toFixed(2)} Gb`;
  if (numeric >= 1e6) return `${(numeric / 1e6).toFixed(2)} Mb`;
  if (numeric >= 1e3) return `${(numeric / 1e3).toFixed(2)} Kb`;
  return `${numeric} bp`;
}

function formatRate(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "-";
  if (numeric <= 1) return `${(numeric * 100).toFixed(2)}%`;
  return `${numeric.toFixed(2)}%`;
}

function formatChartValue(value, formatter) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "-";
  if (formatter === "percent") return formatRate(numeric);
  if (formatter === "bases") return formatBases(numeric);
  return numeric.toFixed(numeric >= 100 ? 0 : 2);
}

function parsePercentDisplay(value) {
  const numeric = Number(String(value ?? "").replace("%", "").trim());
  return Number.isFinite(numeric) ? numeric : null;
}

function parseOptionalNumber(value) {
  if (value === null || value === undefined || value === "") return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
}

function assessAssemblyQuality({ totalLength, contigCount, completeness, contamination }) {
  const hasMetrics = totalLength !== null || contigCount !== null || completeness !== null || contamination !== null;
  if (!hasMetrics) {
    return {
      hasMetrics: false,
      state: "neutral",
      title: "组装结果未生成",
      note: "组装结果未生成",
      body: "当前样本未生成可判读的组装与 CheckM 指标，可能因组装步骤失败或流程提前终止。",
    };
  }

  const reasons = [];
  if (totalLength !== null && totalLength < 500000) reasons.push("总长度不足 500 kb");
  if (contigCount !== null && contigCount > 200) reasons.push("Contig 数过多");
  if (completeness !== null && completeness <= 80) reasons.push("完整性偏低");
  if (contamination !== null && contamination >= 10) reasons.push("污染率偏高");

  const success = (totalLength !== null && totalLength >= 500000)
    && (contigCount !== null && contigCount <= 200)
    && (completeness !== null && completeness > 80)
    && (contamination !== null && contamination < 10);

  if (success) {
    return {
      hasMetrics: true,
      state: "success",
      title: "组装结果整体稳健",
      note: "组装质量良好",
      body: "",
    };
  }

  const highRisk = (totalLength !== null && totalLength < 500000)
    || (contigCount !== null && contigCount > 200)
    || (completeness !== null && completeness <= 50)
    || (contamination !== null && contamination >= 20);

  return {
    hasMetrics: true,
    state: highRisk ? "danger" : "warning",
    title: highRisk ? "组装完整性存在明显风险" : "组装结果需结合质量指标复核",
    note: reasons.length ? reasons.join(" / ") : "建议复核组装结果",
    body: "",
  };
}

function getQualityCardState(q20, q30) {
  if (q20 === null || q30 === null) return "neutral";
  if (q20 !== null && q30 !== null && q20 < 90 && q30 < 80) return "danger";
  return "success";
}

function getQMetricState(label, value) {
  if (value === null) return "neutral";
  if (label === "Q20") return value >= 90 ? "success" : "danger";
  if (label === "Q30") return value >= 80 ? "success" : "danger";
  return "neutral";
}

function getCompletenessState(value) {
  if (value === null) return "neutral";
  if (value > 90) return "success";
  if (value >= 50) return "warning";
  return "danger";
}

function getContaminationState(value) {
  if (value === null) return "neutral";
  if (value < 5) return "success";
  if (value <= 20) return "warning";
  return "danger";
}

function clampPercent(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return 0;
  return Math.max(0, Math.min(100, numeric));
}

function buildChartInsight(text) {
  if (!text) return "";
  return `
    <div class="chart-insight" role="note" aria-label="图表判读">
      <span class="chart-insight-label">图表判读</span>
      <p>${escapeHtml(text)}</p>
    </div>
  `;
}

function summarizeRawQc(section, readLabel) {
  const q20 = Number(section?.before_summary?.q20_rate);
  const q30 = Number(section?.before_summary?.q30_rate);
  const gc = Number(section?.before_summary?.gc_content);
  const meanLength = Number(section?.before_summary?.mean_length);
  const qualityPieces = [];
  if (Number.isFinite(q20) && Number.isFinite(q30)) {
    if (q20 >= 0.9 && q30 >= 0.8) {
      qualityPieces.push(`${readLabel} 原始测序质量整体达标`);
    } else if (q20 < 0.85 || q30 < 0.75) {
      qualityPieces.push(`${readLabel} 原始测序质量偏低，建议重点关注过滤后结果`);
    } else {
      qualityPieces.push(`${readLabel} 原始测序质量接近阈值，建议结合后续质控结果综合判读`);
    }
  }
  if (Number.isFinite(gc)) {
    if (gc >= 0.35 && gc <= 0.65) {
      qualityPieces.push(`GC 比例处于常见范围`);
    } else {
      qualityPieces.push(`GC 比例偏离常见范围`);
    }
  }
  if (Number.isFinite(meanLength)) {
    qualityPieces.push(`平均读长约为 ${meanLength.toFixed(0)} bp`);
  }
  return qualityPieces.join("，") || "";
}

function summarizeInsertSize(fastp) {
  const peak = Number(fastp?.insert_size?.peak);
  const unknown = Number(fastp?.insert_size?.unknown);
  if (!Number.isFinite(peak) && !Number.isFinite(unknown)) return "";
  const summary = [];
  if (Number.isFinite(peak)) {
    summary.push(`插入片段主峰位于 ${peak.toFixed(0)} bp`);
  }
  if (Number.isFinite(unknown)) {
    if (unknown > 500000) {
      summary.push(`未识别配对 reads 较多`);
    } else {
      summary.push(`未识别配对 reads 处于可接受范围`);
    }
  }
  return summary.join("，");
}

function hasFastpReadContent(readSection) {
  if (!readSection || typeof readSection !== "object") return false;
  const beforeSummary = readSection.before_summary || {};
  if (Number.isFinite(Number(beforeSummary.total_reads)) && Number(beforeSummary.total_reads) > 0) {
    return true;
  }
  const qualityCurves = readSection.quality_curves?.mean;
  if (Array.isArray(qualityCurves) && qualityCurves.length) {
    return true;
  }
  const contentCurves = readSection.content_curves || {};
  return ["A", "T", "C", "G", "GC"].some((key) => Array.isArray(contentCurves?.[key]) && contentCurves[key].length);
}

function isPairedEndFastp(fastp, rawQc) {
  const sequencing = String(fastp?.sequencing || "").trim().toLowerCase();
  if (sequencing) {
    if (sequencing.includes("paired")) return true;
    if (sequencing.includes("single")) return false;
    if (sequencing === "pe") return true;
    if (sequencing === "se") return false;
  }
  return hasFastpReadContent(rawQc?.paired_end?.right);
}

function summarizeBaseDistribution(baseCurves, readLabel) {
  const gcValues = Array.isArray(baseCurves?.GC) ? baseCurves.GC : [];
  if (!gcValues.length) return "";
  const gcPercent = gcValues.map((value) => Number(value) * 100).filter(Number.isFinite);
  if (!gcPercent.length) return "";
  const min = Math.min(...gcPercent);
  const max = Math.max(...gcPercent);
  const mean = gcPercent.reduce((sum, value) => sum + value, 0) / gcPercent.length;
  const spread = max - min;
  if (spread <= 8) {
    return `${readLabel} 碱基组成整体平稳，GC 比例均值约为 ${mean.toFixed(1)}%，沿读长波动较小。`;
  }
  return `${readLabel} 碱基组成沿读长存在一定波动，GC 比例均值约为 ${mean.toFixed(1)}%，建议结合接头与过滤结果综合判断。`;
}

function summarizeCoverage(section) {
  const mean = Number(section?.mean_depth);
  const max = Number(section?.max_depth);
  const contigCount = Number(section?.contig_count);
  const notes = [];
  if (Number.isFinite(mean)) {
    if (mean >= 100) {
      notes.push(`整体覆盖深度充足`);
    } else if (mean >= 30) {
      notes.push(`整体覆盖深度基本可用`);
    } else {
      notes.push(`整体覆盖深度偏低`);
    }
  }
  if (Number.isFinite(mean) && Number.isFinite(max) && mean > 0) {
    const ratio = max / mean;
    if (ratio >= 8) {
      notes.push(`局部存在明显高深度峰值`);
    } else {
      notes.push(`覆盖深度分布相对平稳`);
    }
  }
  if (Number.isFinite(contigCount)) {
    notes.push(`共覆盖 ${contigCount.toFixed(0)} 条组装片段`);
  }
  return notes.join("，");
}

function summarizeGeneLengthDistribution(section) {
  const labels = Array.isArray(section?.x_values) ? section.x_values : [];
  const points = Array.isArray(section?.points) ? section.points.map((value) => Number(value) || 0) : [];
  if (!labels.length || !points.length) return "";
  const maxValue = Math.max(...points);
  const maxIndex = points.indexOf(maxValue);
  const total = points.reduce((sum, value) => sum + value, 0);
  const topThree = [...points].sort((a, b) => b - a).slice(0, 3).reduce((sum, value) => sum + value, 0);
  const concentration = total > 0 ? (topThree / total) * 100 : 0;
  const dominantLabel = labels[maxIndex] || "主要长度区间";
  if (concentration >= 55) {
    return `基因长度主要集中于 ${dominantLabel} 等少数区间，长度分布呈明显集中趋势。`;
  }
  return `基因长度以 ${dominantLabel} 区间最为常见，但整体仍呈较分散分布。`;
}

function summarizeContigDepth(section) {
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) return "";
  const grouped = { "基因组": [], "质粒": [] };
  points.forEach((point) => {
    const type = point?.type === "质粒" ? "质粒" : "基因组";
    const depth = Number(point?.depth);
    if (Number.isFinite(depth)) grouped[type].push(depth);
  });
  const avgGenome = grouped["基因组"].length ? grouped["基因组"].reduce((a, b) => a + b, 0) / grouped["基因组"].length : null;
  const avgPlasmid = grouped["质粒"].length ? grouped["质粒"].reduce((a, b) => a + b, 0) / grouped["质粒"].length : null;
  if (avgGenome === null && avgPlasmid === null) return "";
  if (avgGenome !== null && avgPlasmid !== null) {
    if (avgPlasmid > avgGenome * 1.3) {
      return `质粒序列平均深度高于基因组背景，提示部分质粒可能存在较高拷贝特征。`;
    }
    if (avgGenome > avgPlasmid * 1.3) {
      return `基因组序列整体深度高于质粒，提示质粒丰度相对有限。`;
    }
    return `基因组与质粒序列平均深度接近，整体丰度差异不明显。`;
  }
  return avgPlasmid !== null ? `当前结果以质粒序列深度信息为主。` : `当前结果以基因组序列深度信息为主。`;
}

function summarizeRelationship(section) {
  const leftNodes = Array.isArray(section?.nodes_left) ? section.nodes_left : [];
  const rightNodes = Array.isArray(section?.nodes_right) ? section.nodes_right : [];
  const links = Array.isArray(section?.links) ? section.links : [];
  if (!leftNodes.length || !rightNodes.length || !links.length) return "";
  const dominantLeft = leftNodes[0];
  const dominantRight = rightNodes[0];
  const dominantLink = [...links].sort((a, b) => (Number(b.value) || 0) - (Number(a.value) || 0))[0];
  const leftName = dominantLeft?.name || section?.left_label || "主要分类";
  const rightName = dominantRight?.name || section?.right_label || "主要基因";
  if (dominantLink?.source && dominantLink?.target) {
    return `${leftName} 为当前最主要的分类来源，${dominantLink.source} 与 ${dominantLink.target} 的关联最为集中。`;
  }
  return `${leftName} 与 ${rightName} 为当前结果中的主要关系主体。`;
}

function getMetricByKey(metrics, key) {
  return Array.isArray(metrics) ? metrics.find((metric) => metric?.key === key) : null;
}

function buildSummaryCard({ icon, label, title, body, meta, state = "neutral", target = "" }) {
  return `
    <article class="summary-brief-card${target ? " is-jumpable" : ""}" data-state="${escapeHtml(state)}"${target ? ` data-summary-target="${escapeHtml(target)}" tabindex="0" role="button"` : ""}>
      <div class="summary-brief-head">
        <span class="summary-brief-icon" aria-hidden="true">${escapeHtml(icon || "•")}</span>
        <span class="summary-brief-label">${escapeHtml(label)}</span>
      </div>
      <h3>${escapeHtml(title || "-")}</h3>
      <p>${escapeHtml(body || "-")}</p>
      ${meta ? `<span class="summary-brief-meta">${escapeHtml(meta)}</span>` : ""}
    </article>
  `;
}

function bindSummaryCardJumps(container) {
  if (!(container instanceof HTMLElement)) return;
  container.querySelectorAll("[data-summary-target]").forEach((card) => {
    const jump = () => {
      const selector = card.dataset.summaryTarget || "";
      const target = selector ? document.querySelector(selector) : null;
      if (!target) return;
      ensureMobileSectionVisible(target);
      target.scrollIntoView({ behavior: "smooth", block: "start" });
      if (selector.startsWith("#")) {
        history.replaceState(null, "", selector);
      }
    };
    card.addEventListener("click", jump);
    card.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        jump();
      }
    });
  });
}

function escapeRegExp(value) {
  return String(value || "").replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function clearReportSearchHighlights() {
  document.querySelectorAll("mark.report-search-hit").forEach((mark) => {
    const parent = mark.parentNode;
    if (!parent) return;
    parent.replaceChild(document.createTextNode(mark.textContent || ""), mark);
    parent.normalize();
  });
  reportSearchMatches = [];
  reportSearchIndex = -1;
}

function updateReportSearchCounter() {
  const counter = document.getElementById("report-search-count");
  if (!counter) return;
  counter.textContent = reportSearchMatches.length
    ? `${reportSearchIndex + 1} / ${reportSearchMatches.length}`
    : "0 / 0";
}

function setActiveReportSearchMatch(index) {
  if (!reportSearchMatches.length) {
    reportSearchIndex = -1;
    updateReportSearchCounter();
    return;
  }
  reportSearchIndex = (index + reportSearchMatches.length) % reportSearchMatches.length;
  reportSearchMatches.forEach((node, nodeIndex) => {
    node.classList.toggle("is-current", nodeIndex === reportSearchIndex);
  });
  const active = reportSearchMatches[reportSearchIndex];
  active.scrollIntoView({ behavior: "smooth", block: "center" });
  updateReportSearchCounter();
}

function highlightTextNodeForReportSearch(textNode, matcher) {
  const text = textNode.nodeValue || "";
  matcher.lastIndex = 0;
  if (!matcher.test(text)) return [];
  matcher.lastIndex = 0;
  const fragment = document.createDocumentFragment();
  const nodes = [];
  let cursor = 0;
  let match = matcher.exec(text);
  while (match) {
    const start = match.index;
    const value = match[0];
    if (start > cursor) {
      fragment.appendChild(document.createTextNode(text.slice(cursor, start)));
    }
    const mark = document.createElement("mark");
    mark.className = "report-search-hit";
    mark.textContent = value;
    fragment.appendChild(mark);
    nodes.push(mark);
    cursor = start + value.length;
    match = matcher.exec(text);
  }
  if (cursor < text.length) {
    fragment.appendChild(document.createTextNode(text.slice(cursor)));
  }
  textNode.parentNode?.replaceChild(fragment, textNode);
  return nodes;
}

function collectReportSearchTextNodes(root) {
  const nodes = [];
  const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
    acceptNode(node) {
      const text = node.nodeValue || "";
      const parent = node.parentElement;
      if (!text.trim() || !parent) return NodeFilter.FILTER_REJECT;
      if (parent.closest("script, style, mark, .report-search-toolbar")) return NodeFilter.FILTER_REJECT;
      return NodeFilter.FILTER_ACCEPT;
    },
  });
  let node = walker.nextNode();
  while (node) {
    nodes.push(node);
    node = walker.nextNode();
  }
  return nodes;
}

function applyReportSearch(query) {
  clearReportSearchHighlights();
  const keyword = String(query || "").trim();
  const root = document.querySelector(".report-content");
  if (!keyword || !root) {
    updateReportSearchCounter();
    return;
  }
  const matcher = new RegExp(escapeRegExp(keyword), "gi");
  const textNodes = collectReportSearchTextNodes(root);
  reportSearchMatches = textNodes.flatMap((node) => highlightTextNodeForReportSearch(node, matcher));
  setActiveReportSearchMatch(0);
}

function initializeReportSearch() {
  const input = document.getElementById("report-search-input");
  const previous = document.getElementById("report-search-prev");
  const next = document.getElementById("report-search-next");
  const clear = document.getElementById("report-search-clear");
  if (!(input instanceof HTMLInputElement)) return;
  input.addEventListener("input", () => applyReportSearch(input.value));
  input.addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
      event.preventDefault();
      setActiveReportSearchMatch(reportSearchIndex + (event.shiftKey ? -1 : 1));
    }
    if (event.key === "Escape") {
      input.value = "";
      clearReportSearchHighlights();
      updateReportSearchCounter();
    }
  });
  previous?.addEventListener("click", () => setActiveReportSearchMatch(reportSearchIndex - 1));
  next?.addEventListener("click", () => setActiveReportSearchMatch(reportSearchIndex + 1));
  clear?.addEventListener("click", () => {
    input.value = "";
    clearReportSearchHighlights();
    updateReportSearchCounter();
    input.focus();
  });
  updateReportSearchCounter();
}

function countPrioritySerotypeHits(section) {
  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  if (!columns.length || !rows.length) return 0;
  const serotypeIndex = columns.indexOf("血清型");
  const virulenceIndex = columns.indexOf("毒力基因");
  if (serotypeIndex < 0 || virulenceIndex < 0) return 0;
  return rows.filter((row) => {
    const serotype = String(row?.[serotypeIndex] ?? "").trim();
    const virulence = String(row?.[virulenceIndex] ?? "").trim();
    const validSerotype = serotype && serotype !== "-";
    const validVirulence = virulence && virulence !== "-";
    return validSerotype && validVirulence;
  }).length;
}

let currentReportScenario = "research";
let reportSearchMatches = [];
let reportSearchIndex = -1;
const MOBILE_COMPACT_REPORT_KINDS = new Set([
  "sars-cov-2",
  "hmpv",
  "denv",
  "zikav",
  "chikv",
  "hpiv",
  "hadv",
  "norovirus",
  "enterovirus",
  "hepatovirus",
  "bandavirus",
  "orthohantavirus",
  "astroviridae",
  "rhinovirus",
  "seasonal_hcov",
  "rotavirus",
  "rsv",
  "monkeypox",
  "influenza",
]);

function getFirstMatchingCell(row, columns, candidates) {
  if (!Array.isArray(columns) || !Array.isArray(row)) return "";
  for (const name of candidates) {
    const index = columns.indexOf(name);
    if (index >= 0) {
      const value = String(row[index] ?? "").trim();
      if (value) return value;
    }
  }
  return "";
}

function dedupeList(items, limit = 5) {
  const seen = new Set();
  const output = [];
  for (const raw of items || []) {
    const value = String(raw || "").trim();
    if (!value || value === "-" || seen.has(value)) continue;
    seen.add(value);
    output.push(value);
    if (output.length >= limit) break;
  }
  return output;
}

function collectColumnValues(rows, columns, candidates, limit = 200) {
  if (!Array.isArray(rows) || !Array.isArray(columns) || !rows.length || !columns.length) return [];
  const collected = [];
  rows.forEach((row) => {
    candidates.forEach((name) => {
      const value = getFirstMatchingCell(row, columns, [name]);
      if (value) collected.push(value);
    });
  });
  return dedupeList(collected, limit);
}

function classifyClinicalPathogen(speciesName) {
  const species = String(speciesName || "").toLowerCase();
  const commonPathogens = [
    "klebsiella pneumoniae",
    "escherichia coli",
    "staphylococcus aureus",
    "acinetobacter baumannii",
    "pseudomonas aeruginosa",
    "enterococcus faecium",
    "enterococcus faecalis",
    "streptococcus pneumoniae",
    "salmonella enterica",
    "neisseria meningitidis",
    "mycobacterium tuberculosis",
  ];
  const infectionMap = [
    { pattern: ["klebsiella pneumoniae", "escherichia coli"], text: "常见于血流感染、泌尿系统感染、腹腔感染及医院获得性肺炎等场景" },
    { pattern: ["acinetobacter baumannii", "pseudomonas aeruginosa"], text: "常见于呼吸机相关肺炎、血流感染、创面感染及医院获得性感染" },
    { pattern: ["staphylococcus aureus"], text: "常见于血流感染、皮肤软组织感染、肺炎及骨关节感染" },
    { pattern: ["streptococcus pneumoniae"], text: "常见于社区获得性肺炎、菌血症、脑膜炎和鼻窦炎" },
    { pattern: ["neisseria meningitidis"], text: "常见于流行性脑脊髓膜炎和侵袭性血流感染" },
    { pattern: ["salmonella enterica"], text: "常见于肠道感染、菌血症及部分侵袭性感染" },
    { pattern: ["mycobacterium tuberculosis"], text: "常见于肺结核，也可引起播散性感染与肺外结核" },
  ];
  const common = commonPathogens.some((item) => species.includes(item));
  const infectionHint = infectionMap.find((item) => item.pattern.some((keyword) => species.includes(keyword)))?.text
    || "可能的感染类型需结合送检部位、临床表现与基础疾病综合判断";
  return {
    common,
    significance: common ? "属于临床常见致病菌或机会致病菌，具有明确感染相关性" : "临床意义需结合样本来源、分离背景与患者症状进一步判断",
    infectionHint,
  };
}

function extractDominantClinicalSpecies(rawValue) {
  const text = String(rawValue || "").trim();
  if (!text) return "--";
  const cleaned = text
    .replace(/\s+noSpe\b.*$/i, "")
    .replace(/\s+unclassified\b.*$/i, "")
    .trim();
  const match = cleaned.match(/^(.+?)\s*\([^)]*\)/);
  if (match?.[1]) {
    return match[1].trim();
  }
  return cleaned;
}

function classifyPublicHealthPathogen(speciesName) {
  const species = String(speciesName || "").toLowerCase();
  const reportablePatterns = [
    "neisseria meningitidis",
    "vibrio cholerae",
    "salmonella typhi",
    "salmonella paratyphi",
    "shigella",
    "yersinia pestis",
    "bacillus anthracis",
    "brucella",
    "corynebacterium diphtheriae",
    "mycobacterium tuberculosis",
  ];
  const priorityPatterns = [
    "neisseria meningitidis",
    "klebsiella pneumoniae",
    "acinetobacter baumannii",
    "pseudomonas aeruginosa",
    "escherichia coli",
    "salmonella enterica",
    "staphylococcus aureus",
    "mycobacterium tuberculosis",
  ];
  const reportable = reportablePatterns.some((item) => species.includes(item));
  const priority = priorityPatterns.some((item) => species.includes(item));
  let note = "当前更适合作为一般监测对象，是否构成重点公共卫生事件需结合流行病学背景判定。";
  if (reportable) {
    note = "属于法定报告或重点关注病原体，建议结合个案信息及时开展规范化上报与流行病学调查。";
  } else if (priority) {
    note = "属于重点监测或医院感染监测中常见关注病原体，具有持续公共卫生监测意义。";
  }
  return {
    reportable,
    priority,
    note,
  };
}

function normalizeMarkerToken(value) {
  return String(value || "")
    .toLowerCase()
    .replace(/[^a-z0-9\u4e00-\u9fa5]+/g, "");
}

function splitMarkerCandidates(name) {
  return String(name || "")
    .split(/[\/;,，、]|(?:\s+\|\s+)|(?:\s+or\s+)/i)
    .map((item) => item.trim())
    .filter(Boolean);
}

function tokenizeGeneLikeValues(values) {
  const tokens = new Set();
  (Array.isArray(values) ? values : []).forEach((raw) => {
    const text = String(raw || "").trim();
    if (!text || text === "-") return;
    const pieces = text
      .split(/[\/;,，、\s()|\[\]]+/)
      .map((item) => normalizeMarkerToken(item))
      .filter(Boolean);
    pieces.forEach((item) => tokens.add(item));
    const whole = normalizeMarkerToken(text);
    if (whole) tokens.add(whole);
  });
  return tokens;
}

function annotateMarkerHits(markers, detectedGenes) {
  const detected = dedupeList(Array.isArray(detectedGenes) ? detectedGenes : [], 300);
  const detectedTokenSet = tokenizeGeneLikeValues(detected);
  return (Array.isArray(markers) ? markers : []).map((item) => {
    const candidates = splitMarkerCandidates(item?.name);
    const matchedHits = [];
    candidates.forEach((candidate) => {
      const normalizedCandidate = normalizeMarkerToken(candidate);
      if (!normalizedCandidate) return;
      if (!detectedTokenSet.has(normalizedCandidate)) return;
      detected.forEach((gene) => {
        const geneTokens = tokenizeGeneLikeValues([gene]);
        if (geneTokens.has(normalizedCandidate)) {
          matchedHits.push(gene);
        }
      });
    });
    return {
      ...item,
      status: matchedHits.length ? "detected" : "not_detected",
      matched_hits: dedupeList(matchedHits, 20),
    };
  });
}

function matchSerotypeProfiles(profiles, observedSerotype) {
  const text = String(observedSerotype || "").trim().toLowerCase();
  if (!text) return [];
  return (Array.isArray(profiles) ? profiles : []).filter((item) => {
    const pattern = String(item?.pattern || "").trim().toLowerCase();
    return pattern && text.includes(pattern);
  });
}

function extractPrioritySerotypeKnowledge(section) {
  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  if (!columns.length || !rows.length) return [];
  const speciesIndex = columns.indexOf("物种");
  const serotypeIndex = columns.indexOf("血清型");
  const matchedSerotypeIndex = columns.indexOf("知识库命中血清型");
  const panelIndex = columns.indexOf("知识库血清型面板");
  const virulenceIndex = columns.indexOf("血清型-毒力关联");
  const resistanceIndex = columns.indexOf("血清型-耐药关联");
  const regionalIndex = columns.indexOf("血清型-地域分布");
  const interpretationIndex = columns.indexOf("血清型知识库提示");
  const items = [];
  for (const row of rows) {
    if (!Array.isArray(row)) continue;
    const matchedSerotype = String(matchedSerotypeIndex >= 0 ? (row[matchedSerotypeIndex] ?? "") : "").trim();
    const rawSerotype = String(serotypeIndex >= 0 ? (row[serotypeIndex] ?? "") : "").trim();
    if (!matchedSerotype && !rawSerotype) continue;
    items.push({
      species: String(speciesIndex >= 0 ? (row[speciesIndex] ?? "") : "").trim(),
      serotype: matchedSerotype || rawSerotype,
      panel: String(panelIndex >= 0 ? (row[panelIndex] ?? "") : "").trim(),
      virulence: String(virulenceIndex >= 0 ? (row[virulenceIndex] ?? "") : "").trim(),
      resistance: String(resistanceIndex >= 0 ? (row[resistanceIndex] ?? "") : "").trim(),
      regional: String(regionalIndex >= 0 ? (row[regionalIndex] ?? "") : "").trim(),
      interpretation: String(interpretationIndex >= 0 ? (row[interpretationIndex] ?? "") : "").trim(),
    });
  }
  return items;
}

function extractClinicalInterpretation(data) {
  const knowledgeMeta = data?.sections?.knowledge_interpretation?.clinical;
  const metrics = Array.isArray(data?.overview_metrics) ? data.overview_metrics : [];
  const sections = data?.sections || {};
  const speciesMetric = getMetricByKey(metrics, "species_estimation");
  const speciesName = extractDominantClinicalSpecies(speciesMetric?.items?.[0]?.display || speciesMetric?.items?.[1]?.display || "--");
  const classification = classifyClinicalPathogen(speciesName);
  const resistanceSection = sections?.resistance_virulence?.resistance_elements || {};
  const virulenceSection = sections?.resistance_virulence?.virulence_elements || {};
  const resistanceColumns = Array.isArray(resistanceSection.columns) ? resistanceSection.columns : [];
  const virulenceColumns = Array.isArray(virulenceSection.columns) ? virulenceSection.columns : [];
  const resistanceRows = Array.isArray(resistanceSection.rows) ? resistanceSection.rows : [];
  const virulenceRows = Array.isArray(virulenceSection.rows) ? virulenceSection.rows : [];
  const resistanceGenes = dedupeList(resistanceRows.map((row) => getFirstMatchingCell(row, resistanceColumns, ["基因名称", "耐药基因", "Gene", "Best_Hit_ARO"])), 6);
  const resistanceClasses = dedupeList(resistanceRows.map((row) => getFirstMatchingCell(row, resistanceColumns, ["耐药药物", "Drug Class", "Drug", "Resistance Mechanism"])), 6);
  const virulenceGenes = dedupeList(virulenceRows.map((row) => getFirstMatchingCell(row, virulenceColumns, ["毒力基因", "基因名称", "Gene"])), 6);
  const prioritySerotypeKnowledge = extractPrioritySerotypeKnowledge(sections?.priority_serotype || {});
  const mlstKnowledge = sections?.mlst?.knowledge_summary || {};
  const resistanceCount = resistanceRows.length;
  const virulenceCount = virulenceRows.length;
  const priorityCount = countPrioritySerotypeHits(sections?.priority_serotype || {});
  const highResistance = resistanceCount >= 10 || resistanceClasses.length >= 3;
  const highVirulence = virulenceCount >= 10 || priorityCount > 0;
  const riskLevel = highResistance && highVirulence ? "高风险" : (highResistance || highVirulence || resistanceCount >= 3 || virulenceCount >= 3) ? "中风险" : "低风险";
  const qMetric = getMetricByKey(metrics, "q_metrics");
  const q20 = qMetric?.items?.[0]?.display || "--";
  const q30 = qMetric?.items?.[1]?.display || "--";
  const totalBases = getMetricByKey(metrics, "total_bases")?.display || "--";
  const assemblyMetric = getMetricByKey(metrics, "assembly_profile");
  const contigCount = assemblyMetric?.contig_count ?? "--";
  const totalLength = assemblyMetric?.display || assemblyMetric?.items?.[0]?.display || "--";
  const treatmentHint = resistanceClasses.length
    ? `提示可能对${resistanceClasses.join("、")}相关药物类别存在耐药风险，建议结合药敏试验结果综合评估经验性与目标性用药。`
    : "当前未检出明确关键耐药基因，仍建议结合药敏试验与感染部位进行用药决策。";
  const virulenceHint = highVirulence
    ? `检出 ${virulenceCount} 个毒力相关条目，提示该菌株可能具备较强致病潜能。`
    : virulenceCount > 0
      ? `检出 ${virulenceCount} 个毒力相关条目，提示存在一定致病潜能。`
      : "当前未检出明确高风险毒力条目。";
  const fallback = {
    speciesName,
    significance: classification.significance,
    commonLabel: classification.common ? "属于常见临床致病菌" : "并非常见临床致病菌",
    infectionHint: classification.infectionHint,
    resistanceCount,
    virulenceCount,
    resistanceGenes,
    resistanceClasses,
    virulenceGenes,
    riskLevel,
    treatmentHint,
    virulenceHint,
    q20,
    q30,
    totalBases,
    contigCount,
    totalLength,
    taskName: String(data?.task?.name || data?.task?.id || "-"),
    sampleType: String(data?.task?.sample_type || data?.task?.inputtype || "-"),
    submittingUnit: String(data?.task?.owner_group || data?.task?.owner || "-"),
    prioritySerotypeKnowledge,
    mlstKnowledge,
  };
  if (knowledgeMeta?.status === "ready") {
    return {
      ...fallback,
      ...knowledgeMeta,
      resistanceGenes: Array.isArray(knowledgeMeta?.resistanceGenes) ? knowledgeMeta.resistanceGenes : fallback.resistanceGenes,
      resistanceClasses: Array.isArray(knowledgeMeta?.resistanceClasses) ? knowledgeMeta.resistanceClasses : fallback.resistanceClasses,
      virulenceGenes: Array.isArray(knowledgeMeta?.virulenceGenes) ? knowledgeMeta.virulenceGenes : fallback.virulenceGenes,
      evidence: Array.isArray(knowledgeMeta?.evidence) ? knowledgeMeta.evidence : [],
      recommendations: Array.isArray(knowledgeMeta?.recommendations) ? knowledgeMeta.recommendations : [],
      prioritySerotypeKnowledge,
      mlstKnowledge,
    };
  }
  return fallback;
}

function extractCdcInterpretation(data) {
  const knowledgeMeta = data?.sections?.knowledge_interpretation?.cdc;
  const metrics = Array.isArray(data?.overview_metrics) ? data.overview_metrics : [];
  const sections = data?.sections || {};
  const publicHealthSupport = sections?.public_health_support || {};
  const prioritySerotypeKnowledge = extractPrioritySerotypeKnowledge(sections?.priority_serotype || {});
  const mlstKnowledge = sections?.mlst?.knowledge_summary || {};
  const speciesMetric = getMetricByKey(metrics, "species_estimation");
  const speciesName = extractDominantClinicalSpecies(speciesMetric?.items?.[0]?.display || speciesMetric?.items?.[1]?.display || "--");
  const taxonomyRows = Array.isArray(sections?.species_identification?.species?.rows)
    ? sections.species_identification.species.rows
    : [];
  const topRow = taxonomyRows.length
    ? [...taxonomyRows].sort((a, b) => Number(b?.["比例数值"] || 0) - Number(a?.["比例数值"] || 0))[0]
    : null;
  const genus = String(topRow?.["属"] || "").trim() || "--";
  const family = String(topRow?.["科"] || "").trim() || "--";
  const mlst = sections?.mlst?.result || {};
  const serotype = sections?.serotype?.result || {};
  const resistanceRows = sections?.resistance_virulence?.resistance_elements?.rows || [];
  const resistanceColumns = sections?.resistance_virulence?.resistance_elements?.columns || [];
  const virulenceRows = sections?.resistance_virulence?.virulence_elements?.rows || [];
  const virulenceColumns = sections?.resistance_virulence?.virulence_elements?.columns || [];
  const resistanceGenes = collectColumnValues(resistanceRows, resistanceColumns, ["基因名称", "耐药基因", "Gene", "Best_Hit_ARO"], 300);
  const virulenceGenes = collectColumnValues(virulenceRows, virulenceColumns, ["毒力基因", "基因名称", "Gene"], 300);
  const resistanceClasses = dedupeList(resistanceRows.map((row) => getFirstMatchingCell(row, resistanceColumns, ["耐药药物", "Drug Class", "Drug", "Resistance Mechanism"])), 8);
  const pathogenProfile = classifyPublicHealthPathogen(speciesName);
  const supportMatched = publicHealthSupport?.status === "matched";
  const primarySource = String(publicHealthSupport?.primary_source || "").trim();
  const supportPriority = String(publicHealthSupport?.who_priority_group || "").trim();
  const domesticPriority = String(publicHealthSupport?.china_priority_group || "").trim();
  const surveillanceDomain = String(publicHealthSupport?.surveillance_domain || "").trim();
  const supportReportingHint = String(publicHealthSupport?.reporting_hint || "").trim();
  const supportPhenotypes = Array.isArray(publicHealthSupport?.supported_resistance_phenotypes)
    ? publicHealthSupport.supported_resistance_phenotypes
    : [];
  const monitoringTags = Array.isArray(publicHealthSupport?.monitoring_tags)
    ? publicHealthSupport.monitoring_tags
    : [];
  const supportEvidence = Array.isArray(publicHealthSupport?.support_evidence)
    ? publicHealthSupport.support_evidence
    : [];
  const keySerotypes = Array.isArray(publicHealthSupport?.key_serotypes)
    ? publicHealthSupport.key_serotypes
    : [];
  const keyResistanceGenes = Array.isArray(publicHealthSupport?.key_resistance_genes)
    ? publicHealthSupport.key_resistance_genes
    : [];
  const keyVirulenceGenes = Array.isArray(publicHealthSupport?.key_virulence_genes)
    ? publicHealthSupport.key_virulence_genes
    : [];
  const serotypeProfiles = Array.isArray(publicHealthSupport?.serotype_profiles)
    ? publicHealthSupport.serotype_profiles
    : [];
  const clusterHint = String(publicHealthSupport?.cluster_hint || "").trim();
  const typingHint = String(publicHealthSupport?.typing_hint || "").trim();
  const resistanceFocus = String(publicHealthSupport?.resistance_focus || "").trim();
  const geneBlob = resistanceGenes.join(" ").toLowerCase();
  const classBlob = resistanceClasses.join(" ").toLowerCase();
  const hasEsbl = /(ctx-m|shv|tem)/.test(geneBlob) || /cephalosporin|β-内酰胺|beta-lactam/.test(classBlob);
  const hasCarbapenemase = /(kpc|ndm|oxa-48|oxa-23|vim|imp)/.test(geneBlob) || /carbapenem/.test(classBlob);
  const highVirulence = virulenceRows.length >= 10 || /(tox|capsule|invasion|adhesion)/.test(virulenceGenes.join(" ").toLowerCase());
  const typingFragments = [];
  if (mlst?.sequence_type && mlst.sequence_type !== "-") {
    typingFragments.push(`MLST：${mlst.sequence_type}`);
  }
  if (serotype?.predicted_serotype && serotype.predicted_serotype !== "-") {
    typingFragments.push(`血清型：${serotype.predicted_serotype}`);
  }
  const hasTyping = typingFragments.length > 0;
  const observedSerotype = String(serotype?.predicted_serotype || "").trim();
  const matchedSerotypeProfiles = matchSerotypeProfiles(serotypeProfiles, observedSerotype);
  const transmissionRisk = hasCarbapenemase || (pathogenProfile.reportable && highVirulence) || supportPriority === "critical"
    ? "高"
    : (hasEsbl || pathogenProfile.priority || virulenceRows.length > 0 || supportPriority === "high")
      ? "中"
      : "低";
  const outbreakPotential = highVirulence || pathogenProfile.reportable
    ? "提示存在较高聚集性事件或暴发风险，建议关注时空分布与接触史。"
    : "当前未见强烈暴发指征，但仍建议结合地区监测背景持续观察。";
  const historyRelation = hasTyping
    ? `当前样本已获得 ${typingFragments.join("，")}。${typingHint || "如本地历史菌株数据库存在同型别记录，建议进一步开展同源性比较与传播链追踪。"}`
    : (typingHint || "当前未形成足够的分型/相似性证据，暂不能直接判断与本地历史菌株的相关性。");
  const annotatedResistanceMarkers = annotateMarkerHits(keyResistanceGenes, resistanceGenes);
  const annotatedVirulenceMarkers = annotateMarkerHits(keyVirulenceGenes, virulenceGenes);
  const fallback = {
    speciesName,
    genus,
    family,
    reportableLabel: supportReportingHint || (pathogenProfile.reportable ? "法定报告或重点报告病原体" : pathogenProfile.priority ? "重点监测病原体" : "一般监测病原体"),
    pathogenNote: String(publicHealthSupport?.public_health_significance || "").trim() || pathogenProfile.note,
    resistanceCount: resistanceRows.length,
    virulenceCount: virulenceRows.length,
    resistanceGenes,
    virulenceGenes,
    resistanceClasses,
    hasEsbl,
    hasCarbapenemase,
    outbreakPotential,
    historyRelation,
    transmissionRisk,
    typingSummary: typingFragments.length ? typingFragments.join("；") : "未获得稳定分型/相似性结果",
    primarySource,
    surveillanceDomain: surveillanceDomain || "--",
    domesticPriorityGroup: domesticPriority || "--",
    monitoringTags: monitoringTags.join("、"),
    domesticResistanceFocus: resistanceFocus,
    clusterHint: clusterHint || outbreakPotential,
    supportEvidence,
    keySerotypes,
    observedSerotype,
    matchedSerotypeProfiles,
    keyResistanceGenes: annotatedResistanceMarkers,
    keyVirulenceGenes: annotatedVirulenceMarkers,
    whoPriorityGroup: supportPriority || "--",
    whoPhenotypes: supportPhenotypes.length ? supportPhenotypes.join("；") : "",
    monitoringAdvice: transmissionRisk === "高"
      ? "建议尽快开展病例核查、接触者排查及同源性复核，并纳入重点监测。"
      : transmissionRisk === "中"
        ? "建议纳入持续监测，结合时空分布与既往记录评估是否存在传播聚集。"
        : "建议保留本次监测结果并进行常规趋势跟踪。",
    prioritySerotypeKnowledge,
    mlstKnowledge,
  };
  if (knowledgeMeta?.status === "ready") {
    return {
      ...fallback,
      ...knowledgeMeta,
      resistanceGenes: Array.isArray(knowledgeMeta?.resistanceGenes) ? knowledgeMeta.resistanceGenes : fallback.resistanceGenes,
      virulenceGenes: Array.isArray(knowledgeMeta?.virulenceGenes) ? knowledgeMeta.virulenceGenes : fallback.virulenceGenes,
      resistanceClasses: Array.isArray(knowledgeMeta?.resistanceClasses) ? knowledgeMeta.resistanceClasses : fallback.resistanceClasses,
      supportEvidence: Array.isArray(knowledgeMeta?.supportEvidence) ? knowledgeMeta.supportEvidence : fallback.supportEvidence,
      keyResistanceGenes: Array.isArray(knowledgeMeta?.keyResistanceGenes) ? knowledgeMeta.keyResistanceGenes : fallback.keyResistanceGenes,
      keyVirulenceGenes: Array.isArray(knowledgeMeta?.keyVirulenceGenes) ? knowledgeMeta.keyVirulenceGenes : fallback.keyVirulenceGenes,
      prioritySerotypeKnowledge,
      mlstKnowledge,
    };
  }
  return fallback;
}

function buildClinicalScene(data) {
  if (isVirusScenarioReport(data)) {
    return buildVirusClinicalScene(data);
  }
  const meta = extractClinicalInterpretation(data);
  const neisseriaAmr = extractNeisseriaAmrInterpretation(data);
  const shouldShowNeisseriaAmr = isNeisseriaMeningitidisSpecies(meta.speciesName);
  const riskClass = meta.riskLevel === "高风险" ? "level-high" : meta.riskLevel === "中风险" ? "level-mid" : "level-low";
  const reportDate = formatDateTime(data?.task?.updated_at || data?.task?.finished_at || data?.task?.created_at || Date.now());
  const conclusion = meta.conclusion || (meta.riskLevel === "高风险"
    ? `本次单菌测序提示 ${meta.speciesName} 为主要检出病原，具有较高临床相关性，可为感染诊疗和感染控制提供辅助依据。`
    : meta.riskLevel === "中风险"
      ? `本次单菌测序提示 ${meta.speciesName} 具有一定临床相关性，建议结合培养、药敏及临床表现综合评估。`
      : `本次单菌测序检出 ${meta.speciesName}，目前分子特征提示整体风险较低，仍需结合临床表现综合判读。`);
  const commonPathogenNote = meta.commonLabel === "属于常见临床致病菌"
    ? "属于常见临床致病菌，具备明确感染相关性。"
    : "并非常见高频临床致病菌，建议结合送检背景和宿主状态审慎判读。";
  const resistanceSummary = meta.resistanceCount > 0
    ? `检出 ${meta.resistanceCount} 条耐药相关记录，重点基因包括 ${escapeHtml(meta.resistanceGenes.join("、") || "未明确")}。`
    : "当前未检出明确关键耐药基因。";
  const virulenceSummary = meta.virulenceCount > 0
    ? `检出 ${meta.virulenceCount} 条毒力相关记录，代表性因子包括 ${escapeHtml(meta.virulenceGenes.join("、") || "未明确")}。`
    : "当前未检出明确毒力因子。";
  const resistanceClassSummary = meta.resistanceClasses.length
    ? meta.resistanceClasses.join("、")
    : "未见明确重点药物类别";
  const neisseriaAmrSummary = shouldShowNeisseriaAmr && neisseriaAmr.available
    ? [neisseriaAmr.headline, ...neisseriaAmr.interpretationItems.slice(0, 3)].filter(Boolean).join("；")
    : "";
  const evidenceMarkup = Array.isArray(meta.evidence) && meta.evidence.length
    ? `
      <div class="clinical-evidence-box">
        <h5>判读依据</h5>
        <ul class="clinical-guidance-list">
          ${meta.evidence.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
        </ul>
      </div>
    `
    : "";
  const recommendationItems = Array.isArray(meta.recommendations) && meta.recommendations.length
    ? meta.recommendations
    : [
        "建议结合培养、药敏试验及感染灶证据，确认检出病原的真实致病性。",
        meta.treatmentHint,
        meta.riskLevel === "高风险" ? "如病情进展快或存在重症基础，建议优先结合感染科意见动态调整诊疗策略。" : "建议结合临床表现与基础疾病情况，综合判断是否需要进一步强化干预。",
      ];
  const mlstHeadline = String(meta?.mlstKnowledge?.headline || "").trim();
  const mlstItems = Array.isArray(meta?.mlstKnowledge?.items) ? meta.mlstKnowledge.items : [];
  const prioritySerotypes = Array.isArray(meta?.prioritySerotypeKnowledge) ? meta.prioritySerotypeKnowledge : [];
  const mlstFragments = [];
  if (mlstHeadline) mlstFragments.push(mlstHeadline);
  mlstItems.forEach((item) => {
    const parts = [];
    if (item?.lineage_text && item.lineage_text !== "-") parts.push(`克隆复合群/Lineage：${item.lineage_text}`);
    if (Array.isArray(item?.regional) && item.regional.length) parts.push(`流行背景：${item.regional.join("；")}`);
    if (item?.interpretation) parts.push(item.interpretation);
    if (parts.length) mlstFragments.push(parts.join("。"));
  });
  const serotypeFragments = prioritySerotypes.map((item) => {
    const parts = [];
    const typedLabel = [item.species, item.serotype].filter(Boolean).join(" ");
    if (typedLabel) parts.push(`关注血清型/血清群：${typedLabel}`);
    if (item.virulence && item.virulence !== "-") parts.push(`毒力相关性：${item.virulence}`);
    if (item.resistance && item.resistance !== "-") parts.push(`耐药相关性：${item.resistance}`);
    if (item.interpretation && item.interpretation !== "-") parts.push(item.interpretation);
    return parts.join("。");
  }).filter(Boolean);
  return `
    <article class="scene-report-card clinical-scene-card clinical-document">
      <header class="clinical-document-header">
        <div class="clinical-document-kicker">Clinical Auxiliary Diagnostic Report</div>
        <div class="clinical-document-headline">
          <div class="clinical-headline-copy">
            <h3 class="clinical-document-title">临床辅助诊断报告</h3>
            <p class="clinical-document-subtitle">基于单菌测序结果形成的辅助判读意见，仅供临床结合培养、药敏及患者表现综合评估。</p>
          </div>
          <div class="clinical-risk-stamp ${riskClass}">
            <span class="clinical-risk-stamp-label">总体风险</span>
            <strong>${escapeHtml(meta.riskLevel)}</strong>
          </div>
        </div>
        <section class="clinical-report-head">
          <div class="clinical-report-head-item">
            <span>报告编号</span>
            <strong>${escapeHtml(meta.taskName)}</strong>
          </div>
          <div class="clinical-report-head-item">
            <span>送检单位</span>
            <strong>${escapeHtml(meta.submittingUnit || "-")}</strong>
          </div>
          <div class="clinical-report-head-item">
            <span>样本类型</span>
            <strong>${escapeHtml(meta.sampleType || "-")}</strong>
          </div>
          <div class="clinical-report-head-item">
            <span>主要检出病原</span>
            <strong>${escapeHtml(meta.speciesName)}</strong>
          </div>
          <div class="clinical-report-head-item">
            <span>报告日期</span>
            <strong>${escapeHtml(reportDate)}</strong>
          </div>
          <div class="clinical-report-head-item clinical-report-head-item-risk ${riskClass}">
            <span>风险等级</span>
            <strong>${escapeHtml(meta.riskLevel)}</strong>
          </div>
        </section>
      </header>

      <section class="clinical-impression-band">
        <div class="clinical-impression-main">
          <span class="clinical-impression-label">临床结论摘要</span>
          <p>${escapeHtml(conclusion)}</p>
        </div>
        <div class="clinical-impression-tags">
          <span>常见致病菌：${escapeHtml(meta.commonLabel.replace("属于", "").replace("并非", "非"))}</span>
          <span>耐药类别：${escapeHtml(resistanceClassSummary)}</span>
          <span>毒力提示：${escapeHtml(meta.virulenceCount > 0 ? "存在相关风险" : "未见明确高风险提示")}</span>
        </div>
      </section>

      <div class="clinical-report-main">
          <section class="clinical-report-section">
            <div class="clinical-report-section-title">
              <span class="num">一</span>
              <div>
                <h4>病原体鉴定结果</h4>
                <p>病原体名称与临床相关性。</p>
              </div>
            </div>
            <div class="clinical-report-section-body">
              <div class="clinical-keyline">
                <strong>${escapeHtml(meta.speciesName)}</strong>
                <span class="${riskClass}">${escapeHtml(meta.riskLevel)}</span>
              </div>
              <p>${escapeHtml(meta.significance)}</p>
              <p>${escapeHtml(commonPathogenNote)}</p>
              <p>可能相关感染类型：${escapeHtml(meta.infectionHint)}。</p>
            </div>
          </section>

          <section class="clinical-report-section">
            <div class="clinical-report-section-title">
              <span class="num">二</span>
              <div>
                <h4>耐药性分析</h4>
                <p>关键耐药基因与可能耐药药物类别。</p>
              </div>
            </div>
            <div class="clinical-report-section-body">
              <p>${resistanceSummary}</p>
              <p>可能涉及的耐药药物类别：${escapeHtml(resistanceClassSummary)}。</p>
              ${neisseriaAmrSummary ? `<p>脑膜炎奈瑟菌耐药突变临床风险提示：${escapeHtml(neisseriaAmrSummary)}。该结果可用于提示潜在耐药风险与经验性用药谨慎方向，但不替代药敏试验和临床用药决策。</p>` : ""}
              <p>${escapeHtml(meta.treatmentHint)}</p>
            </div>
          </section>

          <section class="clinical-report-section">
            <div class="clinical-report-section-title">
              <span class="num">三</span>
              <div>
                <h4>毒力分析</h4>
                <p>毒力因子检出情况与潜在致病风险。</p>
              </div>
            </div>
            <div class="clinical-report-section-body">
              <p>${virulenceSummary}</p>
              <p>${escapeHtml(meta.virulenceHint)}</p>
            </div>
          </section>

          <section class="clinical-report-section">
            <div class="clinical-report-section-title">
              <span class="num">四</span>
              <div>
                <h4>临床意义解读</h4>
                <p>结论优先的临床判读。</p>
              </div>
            </div>
            <div class="clinical-report-section-body">
              <p>${escapeHtml(conclusion)}</p>
              <p>综合判定：当前样本提示 <strong>${escapeHtml(meta.speciesName)}</strong> 具有 <strong class="${riskClass}">${escapeHtml(meta.riskLevel)}</strong> 临床关注度。</p>
              <p>${mlstFragments.length ? `MLST 分型结果显示：${escapeHtml(mlstFragments.join("；"))}` : "当前未形成可用于临床摘要的 MLST 分型背景提示。"} </p>
              <p>${serotypeFragments.length ? `血清型/血清群结果显示：${escapeHtml(serotypeFragments.join("；"))}` : "当前未形成可用于临床摘要的血清型/血清群知识库提示。"} </p>
              <p>建议结合送检部位、基础疾病、炎症指标、培养与药敏结果综合判断。</p>
              ${evidenceMarkup}
            </div>
          </section>

          <section class="clinical-report-section clinical-report-section-emphasis">
            <div class="clinical-report-section-title">
              <span class="num">五</span>
              <div>
                <h4>诊疗建议</h4>
                <p>原则性建议，不替代临床处方。</p>
              </div>
            </div>
            <div class="clinical-report-section-body">
              <ul class="clinical-guidance-list">
                ${recommendationItems.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
              </ul>
            </div>
          </section>
      </div>
    </article>
  `;
}

function buildCdcScene(data) {
  if (isVirusScenarioReport(data)) {
    return buildVirusCdcScene(data);
  }
  const meta = extractCdcInterpretation(data);
  const neisseriaAmr = extractNeisseriaAmrInterpretation(data);
  const shouldShowNeisseriaAmr = isNeisseriaMeningitidisSpecies(meta.speciesName);
  const resistanceMechanism = meta.hasCarbapenemase
    ? "检出碳青霉烯酶相关耐药机制，提示较高公共卫生传播风险。"
    : meta.hasEsbl
      ? "检出 ESBL 相关耐药特征，提示耐药传播风险需重点监测。"
      : meta.resistanceCount > 0
        ? "检出一定数量耐药相关记录，建议结合地区监测策略持续跟踪。"
        : "当前未见明确重点耐药机制提示。";
  const virulenceAssessment = meta.virulenceCount > 0
    ? `检出 ${meta.virulenceCount} 条毒力相关记录，${meta.outbreakPotential}`
    : "当前未检出明显高致病性毒力信号，仍建议结合流行病学背景综合判断。";
  const evidenceMarkup = Array.isArray(meta.supportEvidence) && meta.supportEvidence.length
    ? `
      <details class="cdc-evidence-panel">
        <summary>支持依据</summary>
        <div class="cdc-evidence-list">
          ${meta.supportEvidence.map((item) => `
            <article class="cdc-evidence-item">
              <div class="cdc-evidence-head">
                <strong>${escapeHtml(item?.manual || "--")}</strong>
                <span>${escapeHtml(item?.rule_type || "--")}</span>
              </div>
              <p>${escapeHtml(item?.basis || "--")}</p>
            </article>
          `).join("")}
        </div>
      </details>
    `
    : "";
  const serotypeMarkup = Array.isArray(meta.matchedSerotypeProfiles) && meta.matchedSerotypeProfiles.length
    ? `
      <div class="cdc-marker-block">
        <h5>实际血清型解释</h5>
        <ul>
          ${meta.matchedSerotypeProfiles.map((item) => `<li><strong>${escapeHtml(meta.observedSerotype || item?.pattern || "--")}</strong>：${escapeHtml(item?.meaning || "--")}</li>`).join("")}
        </ul>
      </div>
    `
    : Array.isArray(meta.keySerotypes) && meta.keySerotypes.length
    ? `
      <div class="cdc-marker-block">
        <h5>重点血清型/分型关注</h5>
        <ul>
          ${meta.keySerotypes.map((item) => `<li><strong>${escapeHtml(item?.name || "--")}</strong>：${escapeHtml(item?.meaning || "--")}</li>`).join("")}
        </ul>
      </div>
    `
    : "";
  const resistanceMarkerMarkup = Array.isArray(meta.keyResistanceGenes) && meta.keyResistanceGenes.length
    ? `
      <div class="cdc-marker-block">
        <h5>重点耐药位点</h5>
        <ul>
          ${meta.keyResistanceGenes.map((item) => `
            <li class="cdc-marker-item ${item?.status === "detected" ? "is-detected" : "is-not-detected"}">
              <div class="cdc-marker-line">
                <strong>${escapeHtml(item?.name || item?.label || (Array.isArray(item?.matched_hits) && item.matched_hits[0]) || "--")}</strong>
                <span class="cdc-marker-status">${item?.status === "detected" ? "已检出" : "未检出"}</span>
              </div>
              <span class="cdc-marker-meaning">${escapeHtml(item?.meaning || "--")}</span>
              ${Array.isArray(item?.matched_hits) && item.matched_hits.length ? `<span class="cdc-marker-hit">实际命中：${escapeHtml(item.matched_hits.join("、"))}</span>` : ""}
            </li>
          `).join("")}
        </ul>
      </div>
    `
    : "";
  const virulenceMarkerMarkup = Array.isArray(meta.keyVirulenceGenes) && meta.keyVirulenceGenes.length
    ? `
      <div class="cdc-marker-block">
        <h5>重点毒力位点</h5>
        <ul>
          ${meta.keyVirulenceGenes.map((item) => `
            <li class="cdc-marker-item ${item?.status === "detected" ? "is-detected" : "is-not-detected"}">
              <div class="cdc-marker-line">
                <strong>${escapeHtml(item?.name || item?.label || (Array.isArray(item?.matched_hits) && item.matched_hits[0]) || "--")}</strong>
                <span class="cdc-marker-status">${item?.status === "detected" ? "已检出" : "未检出"}</span>
              </div>
              <span class="cdc-marker-meaning">${escapeHtml(item?.meaning || "--")}</span>
              ${Array.isArray(item?.matched_hits) && item.matched_hits.length ? `<span class="cdc-marker-hit">实际命中：${escapeHtml(item.matched_hits.join("、"))}</span>` : ""}
            </li>
          `).join("")}
        </ul>
      </div>
    `
    : "";
  const neisseriaAmrSummary = shouldShowNeisseriaAmr && neisseriaAmr.available
    ? [neisseriaAmr.headline, ...neisseriaAmr.interpretationItems.slice(0, 4)].filter(Boolean).join("；")
    : "";
  const mlstHeadline = String(meta?.mlstKnowledge?.headline || "").trim();
  const mlstItems = Array.isArray(meta?.mlstKnowledge?.items) ? meta.mlstKnowledge.items : [];
  const prioritySerotypes = Array.isArray(meta?.prioritySerotypeKnowledge) ? meta.prioritySerotypeKnowledge : [];
  const mlstFragments = [];
  if (mlstHeadline) mlstFragments.push(mlstHeadline);
  mlstItems.forEach((item) => {
    const parts = [];
    if (item?.lineage_text && item.lineage_text !== "-") parts.push(`克隆复合群/Lineage：${item.lineage_text}`);
    if (Array.isArray(item?.regional) && item.regional.length) parts.push(`地域分布：${item.regional.join("；")}`);
    if (item?.interpretation) parts.push(`监测解释：${item.interpretation}`);
    if (parts.length) mlstFragments.push(parts.join("。"));
  });
  const serotypeFragments = prioritySerotypes.map((item) => {
    const parts = [];
    const typedLabel = [item.species, item.serotype].filter(Boolean).join(" ");
    if (typedLabel) parts.push(`命中血清型/血清群：${typedLabel}`);
    if (item.panel && item.panel !== "-") parts.push(`监测面板：${item.panel}`);
    if (item.regional && item.regional !== "-") parts.push(`地域分布：${item.regional}`);
    if (item.interpretation && item.interpretation !== "-") parts.push(`监测解释：${item.interpretation}`);
    return parts.join("。");
  }).filter(Boolean);
  return `
    <article class="scene-report-card cdc-scene-card cdc-document">
      <header class="cdc-document-header">
        <div>
          <div class="cdc-document-kicker">Public Health Surveillance Report</div>
          <h3>疾控监测与风险评估报告</h3>
          <p>面向疾控与流行病学专业人员，强调病原体公共卫生意义、传播风险和监测防控建议。</p>
        </div>
        <div class="cdc-risk-badge risk-${escapeHtml(meta.transmissionRisk)}">
          <span>传播风险</span>
          <strong>${escapeHtml(meta.transmissionRisk)}</strong>
        </div>
      </header>
      <section class="cdc-document-grid">
        <section>
          <h4>一、检测结果概述</h4>
          <p>本次单菌测序主要检出病原体为 <strong>${escapeHtml(meta.speciesName)}</strong>，分类学信息定位至 <strong>${escapeHtml(meta.family)}</strong> / <strong>${escapeHtml(meta.genus)}</strong>。当前判定为 <strong>${escapeHtml(meta.reportableLabel)}</strong>。</p>
          <p>归属监测域：<strong>${escapeHtml(meta.surveillanceDomain || "--")}</strong>；国内监测层级：<strong>${escapeHtml(meta.domesticPriorityGroup || "--")}</strong>${meta.monitoringTags ? `；监测体系：<strong>${escapeHtml(meta.monitoringTags)}</strong>` : ""}。</p>
          <p>WHO 2024 优先级分组：<strong>${escapeHtml(meta.whoPriorityGroup || "--")}</strong>${meta.whoPhenotypes ? `；本次支持的重点耐药表型：<strong>${escapeHtml(meta.whoPhenotypes)}</strong>` : ""}。</p>
        </section>
        <section>
          <h4>二、病原体特征分析</h4>
          <p>${escapeHtml(meta.pathogenNote)}</p>
          <p>当前分型/相似性结果：${escapeHtml(meta.typingSummary)}。</p>
          <p>${mlstFragments.length ? `MLST 分型结果提示：${escapeHtml(mlstFragments.join("；"))}` : "当前未形成可用于疾控摘要的 MLST 分型背景提示。"} </p>
          <p>${serotypeFragments.length ? `血清型/血清群监测提示：${escapeHtml(serotypeFragments.join("；"))}` : "当前未形成可用于疾控摘要的血清型/血清群知识库提示。"} </p>
          ${serotypeMarkup}
        </section>
        <section class="cdc-section-wide">
          <h4>三、耐药与毒力评估</h4>
          <div class="cdc-split-panels">
            <article class="cdc-subpanel">
              <h5>耐药评估</h5>
              <p>${escapeHtml(resistanceMechanism)}</p>
              ${neisseriaAmrSummary ? `<p>脑膜炎奈瑟菌耐药突变监测提示：${escapeHtml(neisseriaAmrSummary)}。如与重点谱系、聚集性事件或异常耐药表型同时出现，建议优先纳入同源性比较和耐药传播链监测。</p>` : ""}
              ${meta.domesticResistanceFocus ? `<p>${escapeHtml(meta.domesticResistanceFocus)}</p>` : ""}
              ${resistanceMarkerMarkup}
            </article>
            <article class="cdc-subpanel">
              <h5>毒力评估</h5>
              <p>${escapeHtml(virulenceAssessment)}</p>
              ${virulenceMarkerMarkup}
            </article>
          </div>
        </section>
        <section class="cdc-section-wide">
          <h4>四、流行病学意义</h4>
          <p>${escapeHtml(meta.historyRelation)}</p>
          <p>${escapeHtml(meta.clusterHint)}</p>
        </section>
        <section>
          <h4>五、风险等级判定</h4>
          <p>综合判定本样本传播风险为 <strong>${escapeHtml(meta.transmissionRisk)}</strong>。${escapeHtml(meta.transmissionRisk === "高" ? "提示需重点关注可能的传播链、聚集性事件或重点耐药株扩散。" : meta.transmissionRisk === "中" ? "提示存在一定传播与扩散风险，建议结合地区监测数据动态评估。" : "当前未见明显传播链或聚集性事件信号。")}</p>
        </section>
        <section>
          <h4>六、监测与防控建议</h4>
          <ul class="cdc-guidance-list">
            <li>${escapeHtml(meta.monitoringAdvice)}</li>
            <li>建议结合病例时空分布、样本来源和接触史，评估是否存在局部传播或聚集性事件。</li>
            <li>如纳入法定报告或重点监测范围，建议按规范及时上报并补充分型/同源性证据。</li>
          </ul>
        </section>
      </section>
      ${evidenceMarkup}
    </article>
  `;
}

function extractResearchInterpretation(data) {
  const sections = data?.sections || {};
  const knowledge = sections?.knowledge_interpretation || {};
  const clinical = knowledge?.clinical || {};
  const cdc = knowledge?.cdc || {};
  const dominant = knowledge?.dominant_species || {};
  const pathogen = knowledge?.pathogen_profile || {};
  const task = data?.task || {};
  const events = Array.isArray(knowledge?.matched_event_rules) ? knowledge.matched_event_rules : [];
  const genes = Array.isArray(knowledge?.matched_gene_rules) ? knowledge.matched_gene_rules : [];
  const supportingSpecies = Array.isArray(knowledge?.supporting_species) ? knowledge.supporting_species : [];
  const intraspeciesSignals = Array.isArray(knowledge?.intraspecies_signals) ? knowledge.intraspecies_signals : [];
  const prioritySerotypeKnowledge = extractPrioritySerotypeKnowledge(sections?.priority_serotype || {});
  const mlstKnowledge = sections?.mlst?.knowledge_summary || {};
  const serotypeKnowledge = sections?.serotype?.knowledge_summary || {};
  const speciesName = String(clinical?.speciesName || pathogen?.common_name || pathogen?.species || dominant?.species || "--").trim() || "--";
  const taxid = String(dominant?.taxid || "--").trim() || "--";
  const scientificName = String(dominant?.scientific_name || pathogen?.species || "--").trim() || "--";
  const pathogenType = String(pathogen?.pathogen_type || "").trim();
  const analysisTarget = String(task?.analysis_target || "").trim().toLowerCase();
  const taskSpecies = String(task?.species || "").trim();
  const riskLevel = String(clinical?.riskLevel || "待判读").trim() || "待判读";
  const eventSummaries = events.map((item) => String(item?.summary || "").trim()).filter(Boolean);
  const geneLabels = genes.map((item) => String(item?.report_label || item?.gene_name || "").trim()).filter(Boolean);
  const hitGenes = genes.flatMap((item) => Array.isArray(item?.matched_hits) ? item.matched_hits : []).filter(Boolean);
  const significance = String(clinical?.significance || cdc?.pathogenNote || `当前结果显示 ${speciesName} 为主要研究对象。`).trim();
  const conclusion = String(clinical?.conclusion || `当前结果支持将 ${speciesName} 视为本样本的主要分子流行病学分析对象。`).trim();
  const nextSteps = Array.isArray(clinical?.recommendations) && clinical.recommendations.length
    ? clinical.recommendations
    : [
        "后续可结合分型、覆盖度和样本来源，对当前知识库命中的稳定性进行进一步验证。",
        "若用于论文或项目汇报，可补充关键位点上下游基因组背景与原始注释证据。",
      ];
  return {
    speciesName,
    taxid,
    scientificName,
    riskLevel,
    significance,
    conclusion,
    eventSummaries: dedupeList(eventSummaries, 4),
    geneLabels: dedupeList(geneLabels, 8),
    hitGenes: dedupeList(hitGenes, 10),
    supportingSpecies,
    intraspeciesSignals,
    prioritySerotypeKnowledge,
    mlstKnowledge,
    serotypeKnowledge,
    nextSteps,
    pathogenType,
    analysisTarget,
    taskSpecies,
    isVirusLike: pathogenType === "病毒" || analysisTarget === "virus",
    hasKnowledge: knowledge?.status === "ready",
    confidenceHint: events.length
      ? "当前结论主要由知识库组合事件与重点位点共同支持，更适合解释为分子特征证据，而非单独的因果归因。"
      : "当前结论主要来自物种身份与重点位点命中，仍需结合更多上下文证据提升解释强度。",
  };
}

function extractNeisseriaAmrInterpretation(data) {
  const amrSection = data?.sections?.mlst?.neisseria_amr || {};
  const headline = String(amrSection?.headline || "").trim();
  const interpretationItems = Array.isArray(amrSection?.interpretation_items)
    ? amrSection.interpretation_items.map((item) => String(item || "").trim()).filter(Boolean)
    : [];
  const highlights = Array.isArray(amrSection?.highlights)
    ? amrSection.highlights.map((item) => String(item || "").trim()).filter(Boolean)
    : [];
  const positiveCount = Number.isFinite(Number(amrSection?.positive_count)) ? Number(amrSection.positive_count) : 0;
  const reviewCount = Number.isFinite(Number(amrSection?.review_count)) ? Number(amrSection.review_count) : 0;
  const available = amrSection?.status === "ready" || headline || interpretationItems.length || highlights.length;
  return {
    available: Boolean(available),
    headline,
    interpretationItems,
    highlights,
    positiveCount,
    reviewCount,
  };
}

function extractTbAmrInterpretation(data) {
  const amrSection = data?.sections?.tb_amr || {};
  const headline = String(amrSection?.headline || "").trim();
  const interpretationItems = Array.isArray(amrSection?.interpretation_items)
    ? amrSection.interpretation_items.map((item) => String(item || "").trim()).filter(Boolean)
    : [];
  const focusCalls = Array.isArray(amrSection?.focus_calls) ? amrSection.focus_calls : [];
  const resistanceGrade = amrSection?.resistance_grade && typeof amrSection.resistance_grade === "object"
    ? amrSection.resistance_grade
    : {};
  const available = amrSection?.status === "ready" || headline || interpretationItems.length || focusCalls.length;
  return {
    available: Boolean(available),
    headline,
    interpretationItems,
    focusCalls,
    resistanceGrade,
    positiveCount: Number.isFinite(Number(amrSection?.positive_count)) ? Number(amrSection.positive_count) : 0,
    reviewCount: Number.isFinite(Number(amrSection?.review_count)) ? Number(amrSection.review_count) : 0,
  };
}

function inferResearchAttentionTone(value) {
  const text = String(value || "").trim().toLowerCase();
  if (!text) return "watch";
  if (["高", "high"].some((token) => text.includes(token))) return "high";
  if (["低", "low", "safe", "敏感"].some((token) => text.includes(token))) return "safe";
  return "watch";
}

function renderResearchAttentionStack(cards) {
  const normalized = Array.isArray(cards)
    ? cards
      .filter((item) => item && typeof item === "object")
      .map((item) => ({
        title: String(item.title || "重点提示").trim(),
        body: String(item.body || "").trim(),
        tone: String(item.tone || "watch").trim() || "watch",
      }))
      .filter((item) => item.body)
    : [];
  if (!normalized.length) return "";
  return `
    <div class="research-attention-stack">
      ${normalized.map((card) => `
        <article class="research-attention-card tone-${escapeHtml(card.tone)}">
          <span class="research-attention-kicker">${escapeHtml(card.title)}</span>
          <p>${escapeHtml(card.body)}</p>
        </article>
      `).join("")}
    </div>
  `;
}

function isNeisseriaMeningitidisSpecies(...values) {
  return values.some((value) => {
    const text = String(value || "").trim().toLowerCase();
    if (!text || text === "--") return false;
    return text.includes("neisseria meningitidis") || text.includes("脑膜炎奈瑟");
  });
}

function extractInfluenzaResearchInterpretation(data) {
  const sections = data?.sections || {};
  const serotype = sections?.serotype || {};
  const coverage = sections?.assembly?.coverage || {};
  const summaryCards = Array.isArray(serotype?.summary_cards) ? serotype.summary_cards : [];
  const mutationRows = Array.isArray(serotype?.mutation_table?.rows) ? serotype.mutation_table.rows : [];
  const mutationColumns = Array.isArray(serotype?.mutation_table?.columns) ? serotype.mutation_table.columns : [];
  const segmentRows = Array.isArray(serotype?.segment_manifest?.rows) ? serotype.segment_manifest.rows : [];
  const segmentColumns = Array.isArray(serotype?.segment_manifest?.columns) ? serotype.segment_manifest.columns : [];
  const speciesMetric = getMetricByKey(Array.isArray(data?.overview_metrics) ? data.overview_metrics : [], "virus_taxonomy");
  const influenzaType = String(serotype?.influenza_type || summaryCards.find((item) => item?.label === "流感类型")?.value || "--").trim() || "--";
  const haSubtype = String(serotype?.ha_subtype || summaryCards.find((item) => item?.label === "HA 亚型")?.value || "--").trim() || "--";
  const naSubtype = String(serotype?.na_subtype || summaryCards.find((item) => item?.label === "NA 亚型")?.value || "--").trim() || "--";
  const subtypeCall = String(serotype?.predicted_serotype || summaryCards.find((item) => item?.label === "分型结果")?.value || "--").trim() || "--";
  const scientificName = String(speciesMetric?.items?.[1]?.display || "Influenza virus").trim() || "Influenza virus";
  const familyGenus = String(speciesMetric?.items?.[2]?.display || "Orthomyxoviridae / -").trim();
  const segmentNameIndex = segmentColumns.findIndex((value) => String(value || "").includes("片段") || String(value || "").toLowerCase().includes("segment"));
  const subtypeIndex = segmentColumns.findIndex((value) => String(value || "").toLowerCase().includes("subtype") || String(value || "").includes("亚型"));
  const mutationSegmentIndex = mutationColumns.findIndex((value) => String(value || "").includes("片段"));
  const mutationPositionIndex = mutationColumns.findIndex((value) => String(value || "").includes("位置"));
  const segmentNames = dedupeList(segmentRows.map((row) => String(segmentNameIndex >= 0 ? row?.[segmentNameIndex] : "").trim()).filter(Boolean), 20);
  const subtypeSegments = segmentRows
    .filter((row) => subtypeIndex >= 0 && String(row?.[subtypeIndex] || "").trim() && String(row?.[subtypeIndex] || "").trim() !== "-")
    .map((row) => {
      const name = String(segmentNameIndex >= 0 ? row?.[segmentNameIndex] : "").trim() || "--";
      const subtype = String(row?.[subtypeIndex] || "").trim() || "--";
      return `${name}=${subtype}`;
    });
  const topMutations = mutationRows.slice(0, 8).map((row) => {
    const seg = String(mutationSegmentIndex >= 0 ? row?.[mutationSegmentIndex] : "").trim() || "--";
    const pos = String(mutationPositionIndex >= 0 ? row?.[mutationPositionIndex] : "").trim() || "--";
    const ref = String(row?.[2] || "").trim();
    const alt = String(row?.[3] || "").trim();
    return `${seg}:${pos} ${ref}>${alt}`.replace(/\s+>/, " >");
  });
  const knowledge = extractViralKnowledgeSummary(serotype);
  return {
    speciesName: scientificName,
    influenzaType,
    haSubtype,
    naSubtype,
    subtypeCall,
    familyGenus,
    segmentCount: segmentRows.length,
    segmentNames,
    subtypeSegments,
    mutationCount: mutationRows.length,
    topMutations,
    knowledgeHeadline: knowledge.headline,
    knowledgeFragments: knowledge.fragments,
    knowledgeTemplate: serotype?.report_template && typeof serotype.report_template === "object" ? serotype.report_template : {},
    meanDepth: String(coverage?.mean_depth ?? "--"),
    coverageFraction: formatRate(Number(coverage?.coverage_fraction || 0)),
    coverage10x: formatRate(Number(coverage?.coverage_10x_fraction || 0)),
    coverage100x: formatRate(Number(coverage?.coverage_100x_fraction || 0)),
    segmentCoverageHints: Array.isArray(coverage?.segments)
      ? coverage.segments
          .slice()
          .sort((left, right) => Number(right?.mean_depth || 0) - Number(left?.mean_depth || 0))
          .slice(0, 4)
          .map((item) => `${item?.name || "--"}（深度 ${item?.mean_depth ?? "--"}，覆盖 ${formatRate(Number(item?.coverage_fraction || 0))}）`)
      : [],
  };
}

function extractMonkeypoxResearchInterpretation(data) {
  const sections = data?.sections || {};
  const serotype = sections?.serotype || {};
  const coverage = sections?.assembly?.coverage || {};
  const summaryCards = Array.isArray(serotype?.summary_cards) ? serotype.summary_cards : [];
  const mutationRows = Array.isArray(serotype?.mutation_table?.rows) ? serotype.mutation_table.rows : [];
  const mutationColumns = Array.isArray(serotype?.mutation_table?.columns) ? serotype.mutation_table.columns : [];
  const speciesMetric = getMetricByKey(Array.isArray(data?.overview_metrics) ? data.overview_metrics : [], "virus_taxonomy");
  const clade = String(serotype?.predicted_clade || summaryCards.find((item) => item?.label === "Clade")?.value || "--").trim() || "--";
  const lineage = String(serotype?.predicted_lineage || summaryCards.find((item) => item?.label === "Lineage")?.value || "--").trim() || "--";
  const outbreak = String(serotype?.predicted_outbreak || summaryCards.find((item) => item?.label === "Outbreak")?.value || "--").trim() || "--";
  const speciesName = String(speciesMetric?.items?.[1]?.display || "Monkeypox virus").trim() || "Monkeypox virus";
  const familyGenus = String(speciesMetric?.items?.[2]?.display || "Poxviridae / Orthopoxvirus").trim() || "Poxviridae / Orthopoxvirus";
  const qcMetric = Array.isArray(serotype?.quality_metrics)
    ? serotype.quality_metrics.find((item) => String(item?.label || "").toLowerCase().includes("qc"))
    : null;
  const qualityLabel = String(qcMetric?.value || "--").trim() || "--";
  const geneIndex = mutationColumns.indexOf("基因名");
  const posIndex = mutationColumns.indexOf("位置");
  const hgvsPIndex = mutationColumns.indexOf("HGVS.p");
  const qualityIndex = mutationColumns.indexOf("质量分层");
  const highRows = qualityIndex >= 0
    ? mutationRows.filter((row) => String(row?.[qualityIndex] || "").trim() === "高质量突变")
    : mutationRows;
  const lowRows = qualityIndex >= 0
    ? mutationRows.filter((row) => String(row?.[qualityIndex] || "").trim() === "低质量突变")
    : [];
  const topMutations = highRows.slice(0, 8).map((row) => {
    const gene = String(geneIndex >= 0 ? row?.[geneIndex] : "").trim() || "--";
    const pos = String(posIndex >= 0 ? row?.[posIndex] : "").trim() || "--";
    const hgvsP = String(hgvsPIndex >= 0 ? row?.[hgvsPIndex] : "").trim();
    return hgvsP ? `${gene}:${pos} ${hgvsP}` : `${gene}:${pos}`;
  });
  const knowledge = extractViralKnowledgeSummary(serotype);
  return {
    speciesName,
    familyGenus,
    clade,
    lineage,
    outbreak,
    qualityLabel,
    mutationCount: mutationRows.length,
    highMutationCount: highRows.length,
    lowMutationCount: lowRows.length,
    topMutations,
    knowledgeHeadline: knowledge.headline,
    knowledgeFragments: knowledge.fragments,
    meanDepth: String(coverage?.mean_depth ?? "--"),
    coverageFraction: formatRate(Number(coverage?.coverage_fraction || 0)),
    coverage10x: formatRate(Number(coverage?.coverage_10x_fraction || 0)),
    coverage100x: formatRate(Number(coverage?.coverage_100x_fraction || 0)),
  };
}

function extractViralKnowledgeSummary(section) {
  const knowledge = section?.knowledge_summary && typeof section.knowledge_summary === "object"
    ? section.knowledge_summary
    : {};
  const headline = String(knowledge?.headline || "").trim();
  const items = Array.isArray(knowledge?.items) ? knowledge.items : [];
  const fragments = items.slice(0, 3).map((item) => {
    const serotype = String(item?.serotype || item?.matched_on || "").trim() || "--";
    const interpretation = String(item?.interpretation || "").trim();
    const virulence = Array.isArray(item?.virulence) ? item.virulence.filter(Boolean).slice(0, 2).join(" / ") : "";
    const regional = Array.isArray(item?.regional) ? item.regional.filter(Boolean).slice(0, 2).join(" / ") : "";
    const detail = [interpretation, virulence, regional].filter(Boolean).join("；");
    return detail ? `${serotype}：${detail}` : serotype;
  });
  return { headline, fragments };
}

function normalizeKnowledgeNarrativeText(value) {
  return String(value || "")
    .replace(/[。]+$/g, "")
    .replace(/\s+/g, " ")
    .trim();
}

function isGenericKnowledgeMatchHeadline(value) {
  const text = normalizeKnowledgeNarrativeText(value);
  return /^当前病毒分型结果已命中知识库(中的.*条目|条目)?$/u.test(text);
}

function renderKnowledgeNarrative(meta, fallbackLabel) {
  const headline = isGenericKnowledgeMatchHeadline(meta?.knowledgeHeadline)
    ? ""
    : normalizeKnowledgeNarrativeText(meta?.knowledgeHeadline);
  const rawFragments = Array.isArray(meta?.knowledgeFragments) ? meta.knowledgeFragments.filter(Boolean) : [];
  const fragments = [];
  rawFragments.forEach((fragment) => {
    const cleaned = normalizeKnowledgeNarrativeText(fragment);
    if (!cleaned) return;
    if (headline && (headline.includes(cleaned) || cleaned.includes(headline))) return;
    if (fragments.some((item) => item === cleaned || item.includes(cleaned) || cleaned.includes(item))) return;
    fragments.push(cleaned);
  });
  if (headline && fragments.length) {
    return `结合知识库既有分型资料，本次结果与以下研究背景具有一致性：${escapeHtml(fragments.join("；"))}。${escapeHtml(headline)}`;
  }
  if (headline) {
    return `结合知识库既有分型资料，${escapeHtml(headline)}。`;
  }
  if (fragments.length) {
    return `知识库既往分型资料显示，本次结果更接近以下研究背景：${escapeHtml(fragments.join("；"))}。`;
  }
  return `当前知识库中尚缺少与该${escapeHtml(fallbackLabel)}分型直接对应的背景资料，因此现阶段仍以系统发育位置、历史参考株及原始比对证据作为主要解释依据。`;
}

function extractRespiratoryNextcladeResearchInterpretation(data) {
  const sections = data?.sections || {};
  const serotype = sections?.serotype || {};
  const coverage = sections?.assembly?.coverage || {};
  const summaryCards = Array.isArray(serotype?.summary_cards) ? serotype.summary_cards : [];
  const mutationTable = serotype?.mutation_table && typeof serotype.mutation_table === "object"
    ? serotype.mutation_table
    : { rows: [], columns: [] };
  const mutationRows = Array.isArray(mutationTable?.rows) ? mutationTable.rows : [];
  const mutationColumns = Array.isArray(mutationTable?.columns) ? mutationTable.columns : [];
  const mode = resolveSerotypeMode(serotype);
  const isHadv = mode === "hadv_typing";
  const isHpiv = mode === "hpiv_typing";
  const isHmpv = mode === "hmpv_nextclade";
  const isDenv = mode === "denv_nextclade";
  const isZikav = mode === "zikav_nextclade";
  const isChikv = mode === "chikv_nextclade";
  const isEbola = mode === "ebola_nextclade";
  const isNorovirus = mode === "norovirus_typing";
  const isEnterovirus = mode === "enterovirus_typing";
  const isHiv = mode === "hiv_resistance";
  const isHepatovirus = mode === "hepatovirus_typing";
  const isBandavirus = mode === "bandavirus_typing";
  const isOrthohantavirus = mode === "orthohantavirus_typing";
  const isAstroviridae = mode === "astroviridae_typing";
  const isRhinovirus = mode === "rhinovirus_typing";
  const isSeasonalHcov = mode === "seasonal_hcov_typing";
  const isRotavirus = mode === "rotavirus_typing";
  const hepatovirusBroad = String(serotype?.predicted_group || "").trim().toUpperCase();
  const hepatovirusSpecies = ({ HAV: "Hepatitis A virus", HBV: "Hepatitis B virus", HCV: "Hepatitis C virus", HDV: "Hepatitis D virus", HEV: "Hepatitis E virus" })[hepatovirusBroad] || "Hepatitis virus";
  const hepatovirusLabel = ({ HAV: "甲型肝炎病毒", HBV: "乙型肝炎病毒", HCV: "丙型肝炎病毒", HDV: "丁型肝炎病毒", HEV: "戊型肝炎病毒" })[hepatovirusBroad] || "肝炎病毒";
  const virusShort = isHiv ? "HIV" : (isRotavirus ? "RotaV" : (isNorovirus ? "NoV" : (isEnterovirus ? "EV" : (isHepatovirus ? (hepatovirusBroad || "HepV") : (isBandavirus ? "BandV" : (isOrthohantavirus ? "HTNV" : (isEbola ? "EBOV" : (isAstroviridae ? "AstV" : (isRhinovirus ? "HRV" : (isSeasonalHcov ? "HCoV" : (isChikv ? "CHIKV" : (isZikav ? "ZIKV" : (isDenv ? "DENV" : (isHmpv ? "HMPV" : (isHpiv ? "HPIV" : (isHadv ? "HAdV" : "RSV"))))))))))))))));
  const virusLabel = isHiv ? "HIV" : (isRotavirus ? "轮状病毒" : (isNorovirus ? "诺如病毒" : (isEnterovirus ? "肠道病毒" : (isHepatovirus ? hepatovirusLabel : (isBandavirus ? "班达病毒" : (isOrthohantavirus ? "汉坦病毒" : (isEbola ? "埃博拉病毒" : (isAstroviridae ? "星状病毒" : (isRhinovirus ? "鼻病毒" : (isSeasonalHcov ? "季节性冠状病毒" : (isChikv ? "基孔肯雅病毒" : (isZikav ? "寨卡病毒" : (isDenv ? "登革热病毒" : (isHmpv ? "人偏肺病毒" : (isHpiv ? "人副流感病毒" : (isHadv ? "人腺病毒" : "呼吸道合胞病毒"))))))))))))))));
  const speciesName = isHiv ? "Human immunodeficiency virus" : (isRotavirus ? "Human rotavirus" : (isNorovirus ? "Norovirus" : (isEnterovirus ? "Human enterovirus" : (isHepatovirus ? hepatovirusSpecies : (isBandavirus ? "Bandavirus dabieense" : (isOrthohantavirus ? "Orthohantavirus" : (isEbola ? "Ebola virus" : (isAstroviridae ? "Human astrovirus" : (isRhinovirus ? "Human rhinovirus" : (isSeasonalHcov ? "Seasonal human coronavirus" : (isChikv ? "Chikungunya virus" : (isZikav ? "Zika virus" : (isDenv ? "Dengue virus" : (isHmpv ? "Human metapneumovirus" : (isHpiv ? "Human parainfluenza virus" : (isHadv ? "Human adenovirus" : "Respiratory syncytial virus"))))))))))))))));
  const familyGenus = isHiv ? "Retroviridae / Lentivirus" : (isRotavirus ? "Reoviridae / Rotavirus" : (isNorovirus ? "Caliciviridae / Norovirus" : (isEnterovirus ? "Picornaviridae / Enterovirus" : (isHepatovirus ? "Picornaviridae / Hepatovirus" : (isBandavirus ? "Phenuiviridae / Bandavirus" : (isOrthohantavirus ? "Hantaviridae / Orthohantavirus" : (isEbola ? "Filoviridae / Orthoebolavirus" : (isAstroviridae ? "Astroviridae / Mamastrovirus-Avastrovirus" : (isRhinovirus ? "Picornaviridae / Enterovirus" : (isSeasonalHcov ? "Coronaviridae / Alphacoronavirus-Betacoronavirus" : (isChikv ? "Togaviridae / Alphavirus" : (isZikav ? "Flaviviridae / Flavivirus" : (isDenv ? "Flaviviridae / Flavivirus" : (isHmpv ? "Pneumoviridae / Metapneumovirus" : (isHpiv ? "Paramyxoviridae / Orthorubulavirus-Respirovirus" : (isHadv ? "Adenoviridae / Mastadenovirus" : "Pneumoviridae / Orthopneumovirus"))))))))))))))));
  const clade = String(serotype?.predicted_clade || summaryCards.find((item) => item?.label === "Nextclade Clade" || item?.label === "HPIV 亚型" || item?.label === "HAdV 分型" || item?.label === "双位点分型" || item?.label === "VP1 分型" || item?.label === "子亚型" || item?.label === "HAV 子亚型" || item?.label === "HAV子亚型" || item?.label === "ORF2 分型" || item?.label === "大类分型" || item?.label === "大组分型" || item?.label === "Orthohantavirus 分型")?.value || "--").trim() || "--";
  const group = String(serotype?.predicted_group || summaryCards.find((item) => item?.label === "病毒属" || item?.label === "大亚型" || item?.label === "S片段分型")?.value || "--").trim() || "--";
  const lineage = String(serotype?.predicted_lineage || summaryCards.find((item) => String(item?.label || "").includes("Lineage") || item?.label === "病毒种")?.value || "--").trim() || "--";
  const qcMetric = Array.isArray(serotype?.quality_metrics)
    ? serotype.quality_metrics.find((item) => String(item?.label || "").toLowerCase().includes("qc"))
    : null;
  const qualityLabel = String(qcMetric?.value || "--").trim() || "--";
  const geneIndex = mutationColumns.indexOf("基因名");
  const posIndex = mutationColumns.indexOf("位置");
  const hgvsPIndex = mutationColumns.indexOf("HGVS.p");
  const qualityIndex = mutationColumns.indexOf("质量分层");
  const highRows = qualityIndex >= 0
    ? mutationRows.filter((row) => String(row?.[qualityIndex] || "").trim() === "高质量突变")
    : mutationRows;
  const lowRows = qualityIndex >= 0
    ? mutationRows.filter((row) => String(row?.[qualityIndex] || "").trim() === "低质量突变")
    : [];
  const topMutations = highRows.slice(0, 8).map((row) => {
    const gene = String(geneIndex >= 0 ? row?.[geneIndex] : "").trim() || "--";
    const pos = String(posIndex >= 0 ? row?.[posIndex] : "").trim() || "--";
    const hgvsP = String(hgvsPIndex >= 0 ? row?.[hgvsPIndex] : "").trim();
    return hgvsP ? `${gene}:${pos} ${hgvsP}` : `${gene}:${pos}`;
  });
  const knowledge = extractViralKnowledgeSummary(serotype);
  return {
    virusShort,
    virusLabel,
    speciesName,
    familyGenus,
    isHadv,
    isHiv,
    isNorovirus,
    isEnterovirus,
    isHepatovirus,
    isBandavirus,
    isOrthohantavirus,
    isAstroviridae,
    isRhinovirus,
    isSeasonalHcov,
    isRotavirus,
    clade,
    group,
    lineage,
    isHpiv,
    qualityLabel,
    mutationCount: mutationRows.length,
    highMutationCount: highRows.length,
    lowMutationCount: lowRows.length,
    topMutations,
    knowledgeHeadline: knowledge.headline,
    knowledgeFragments: knowledge.fragments,
    meanDepth: String(coverage?.mean_depth ?? "--"),
    coverageFraction: formatRate(Number(coverage?.coverage_fraction || 0)),
    coverage10x: formatRate(Number(coverage?.coverage_10x_fraction || 0)),
    coverage100x: formatRate(Number(coverage?.coverage_100x_fraction || 0)),
  };
}

function isKnownVirusReport(data) {
  return isSarsCov2NextcladeReport(data)
    || isHmpvNextcladeReport(data)
    || isDenvNextcladeReport(data)
    || isZikavNextcladeReport(data)
    || isChikvNextcladeReport(data)
    || isEbolaNextcladeReport(data)
    || isHpivTypingReport(data)
    || isHadvTypingReport(data)
    || isNorovirusTypingReport(data)
    || isEnterovirusTypingReport(data)
    || isHepatovirusTypingReport(data)
    || isBandavirusTypingReport(data)
    || isOrthohantavirusTypingReport(data)
    || isOrthoebolavirusTypingReport(data)
    || isAstroviridaeTypingReport(data)
    || isRhinovirusTypingReport(data)
    || isSeasonalHcovTypingReport(data)
    || isRotavirusTypingReport(data)
    || isRsvNextcladeReport(data)
    || isMonkeypoxNextcladeReport(data)
    || isInfluenzaTypingReport(data)
    || isHivTypingReport(data);
}

function isSingleBacteriaReport(data) {
  return String(data?.task?.method || "").trim() !== "meta"
    && String(data?.task?.analysis_target || "bacteria").trim() === "bacteria";
}

function isVirusScenarioReport(data) {
  return String(data?.task?.method || "").trim() !== "meta"
    && (String(data?.task?.analysis_target || "").trim() === "virus" || isKnownVirusReport(data));
}

function supportsReportScenario(data) {
  return isSingleBacteriaReport(data) || isVirusScenarioReport(data);
}

function formatVirusTypingResult(meta) {
  if (!meta) return "--";
  if (meta.kind === "influenza") {
    return [meta.influenzaType, meta.haSubtype, meta.naSubtype].filter((value) => value && value !== "--").join(" / ") || meta.subtypeCall || "--";
  }
  if (meta.kind === "monkeypox") {
    return [meta.clade, meta.lineage, meta.outbreak].filter((value) => value && value !== "--").join(" / ") || "--";
  }
  if (meta.isHepatovirus) {
    return [meta.group, meta.clade].filter((value) => value && value !== "--").join(" / ") || "--";
  }
  if (meta.isRotavirus || meta.isAstroviridae || meta.isBandavirus || meta.isOrthohantavirus) {
    return [meta.clade, meta.group, meta.lineage].filter((value) => value && value !== "--").join(" / ") || "--";
  }
  return [meta.clade, meta.lineage].filter((value) => value && value !== "--").join(" / ") || meta.clade || "--";
}

function applyKnowledgeVirusReportTemplate(fallback, knowledgeTemplate) {
  if (!knowledgeTemplate || typeof knowledgeTemplate !== "object") return fallback;
  const status = String(knowledgeTemplate.status || "").trim();
  if (status && status !== "ready") return fallback;
  const readText = (camelKey, snakeKey) => String(knowledgeTemplate?.[camelKey] || knowledgeTemplate?.[snakeKey] || "").trim();
  const readList = (camelKey, snakeKey) => {
    const value = Array.isArray(knowledgeTemplate?.[camelKey])
      ? knowledgeTemplate[camelKey]
      : (Array.isArray(knowledgeTemplate?.[snakeKey]) ? knowledgeTemplate[snakeKey] : []);
    return value.map((item) => String(item || "").trim()).filter(Boolean);
  };
  return {
    ...fallback,
    templateSource: readText("source", "source") || fallback.templateSource || "",
    templateId: readText("id", "id") || fallback.templateId || "",
    category: readText("category", "category") || fallback.category,
    clinicalRisk: readText("clinicalRisk", "clinical_risk") || fallback.clinicalRisk,
    cdcRisk: readText("cdcRisk", "cdc_risk") || fallback.cdcRisk,
    evidenceBasis: readText("evidenceBasis", "evidence_basis") || fallback.evidenceBasis,
    clinicalMeaning: readText("clinicalMeaning", "clinical_meaning") || fallback.clinicalMeaning,
    cdcMeaning: readText("cdcMeaning", "cdc_meaning") || fallback.cdcMeaning,
    clinicalRecommendations: readList("clinicalRecommendations", "clinical_recommendations").length
      ? readList("clinicalRecommendations", "clinical_recommendations")
      : fallback.clinicalRecommendations,
    cdcRecommendations: readList("cdcRecommendations", "cdc_recommendations").length
      ? readList("cdcRecommendations", "cdc_recommendations")
      : fallback.cdcRecommendations,
  };
}

function resolveVirusReportTemplate(meta, knowledgeText = "") {
  const isVectorBorne = ["DENV", "ZIKV", "CHIKV"].includes(String(meta?.virusShort || "").trim());
  const isRespiratory = meta?.kind === "influenza"
    || meta?.isHpiv
    || meta?.isHadv
    || meta?.isRhinovirus
    || meta?.isSeasonalHcov
    || ["RSV", "HMPV"].includes(String(meta?.virusShort || "").trim());
  const isEnteric = meta?.isNorovirus || meta?.isEnterovirus || meta?.isRotavirus || meta?.isAstroviridae;
  const isHepatitis = meta?.isHepatovirus;
  const isNaturalFocus = meta?.isBandavirus || meta?.isOrthohantavirus;
  const isBloodborne = meta?.isHiv || isHepatitis;
  const base = {
    category: "general",
    clinicalRisk: "低-中风险",
    cdcRisk: "低",
    evidenceBasis: "参考分型、覆盖度、质量指标、突变位点与 IGV 比对证据",
    clinicalMeaning: knowledgeText || `当前结果提示样本中存在 ${meta?.virusLabel || "病毒"} 相关分子证据；是否构成活动性感染需结合样本类型、采样时间、临床表现和必要的复核检测。`,
    cdcMeaning: "该病毒结果可作为专项监测线索，建议结合时空聚集、样本来源和历史同型别记录判断公共卫生处置优先级。",
    clinicalRecommendations: [
      "建议结合样本类型、采样时间、症状体征、影像学/实验室指标及必要的复核检测综合判断临床相关性。",
      "若用于临床处置，本报告提供的是病原分型和分子证据，不替代法定诊断标准、病原学确认或临床医嘱。",
      "如结果与临床表现不一致，建议复核原始 reads、覆盖度和关键位点，并结合其他检测方法确认。",
    ],
    cdcRecommendations: [
      "建议结合病例发病时间、采样地点、旅行/暴露史和密切接触信息，判断输入性、聚集性或本地传播可能。",
      "建议将分型结果与本地历史序列或上级平台数据进行比对，必要时补充系统发育或同源性分析。",
      "若涉及法定报告、重点传染病或聚集性事件，应按现行规范及时上报并补充流调证据。",
    ],
  };
  if (meta?.kind === "influenza") {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "respiratory",
      clinicalRisk: meta?.subtypeCall && meta.subtypeCall !== "--" ? "中风险" : "待复核",
      cdcRisk: meta?.subtypeCall && meta.subtypeCall !== "--" ? "中" : "待评估",
      evidenceBasis: "流感类型、HA/NA 分型、segment 组成、覆盖度与突变注释",
      clinicalMeaning: "流感病毒检出具有明确呼吸道感染相关性，但病情轻重、传染期和治疗决策仍需结合症状、采样时间、抗原/核酸复核及基础疾病综合判断。",
      cdcMeaning: "流感结果应纳入呼吸道传染病季节性监测，重点关注亚型变化、聚集性病例和重症病例比例。",
      clinicalRecommendations: [
        "建议结合发病时间窗、重症风险因素和当地诊疗规范评估抗病毒治疗时机。",
        "建议结合呼吸道症状、抗原/核酸复核结果和采样质量判断检出结果的临床相关性。",
        "如为重症、聚集性或特殊人群样本，建议优先复核分型、segment 覆盖度和关键突变位点。",
      ],
      cdcRecommendations: [
        "建议纳入流感季节性监测，关注亚型变化、重症比例和聚集性病例。",
        "建议与本地哨点监测、疫苗株背景和历史序列进行比对，判断是否存在异常谱系变化。",
        "如出现学校、养老机构或医院聚集性病例，建议补充分型/同源性证据并按规范处置。",
      ],
    }, meta?.knowledgeTemplate);
  }
  if (meta?.kind === "monkeypox") {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "contact-transmitted",
      clinicalRisk: meta?.clade && meta.clade !== "--" ? "中风险" : "待复核",
      cdcRisk: meta?.clade && meta.clade !== "--" ? "中" : "待评估",
      evidenceBasis: "Nextclade clade/lineage、覆盖度、QC 与突变位点",
      clinicalMeaning: "猴痘病毒检出需结合皮疹、发热、暴露史与采样部位判断临床相关性；分型结果可辅助追踪传播背景，但不替代临床诊断流程。",
      cdcMeaning: "猴痘结果具有公共卫生追踪价值，建议结合个案调查、接触者管理和本地历史序列开展传播链评估。",
      clinicalRecommendations: [
        "建议结合皮疹部位、病程阶段、暴露史和采样类型评估检出结果的临床意义。",
        "建议必要时复核关键位点和覆盖度，避免因低覆盖或混样造成谱系解释偏差。",
        "临床处置仍应依据现行诊疗规范和感染防控要求执行。",
      ],
      cdcRecommendations: [
        "建议结合个案调查、接触者管理和活动轨迹评估传播链。",
        "建议将 clade/lineage 与本地及上级平台历史序列比对，判断是否为既有传播链延续。",
        "如存在聚集性或跨区域关联，应补充同源性分析和暴露网络信息。",
      ],
    }, meta?.knowledgeTemplate);
  }
  if (isHepatitis) {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "hepatitis",
      clinicalRisk: "中风险",
      cdcRisk: "中",
      evidenceBasis: "肝炎病毒 broad 大亚型、子亚型/基因型参考竞争、覆盖度与突变注释",
      clinicalMeaning: knowledgeText || "肝炎病毒分型结果可辅助判断病毒类别与分子流行病学背景；临床诊断仍需结合肝功能、血清学标志物、病毒载量、病程阶段和既往感染/免疫史综合判断。",
      cdcMeaning: "肝炎病毒结果应重点关注传播途径、感染来源和同型别聚集情况；不同大亚型对应的流调问题不同，不应只按通用病毒阳性处理。",
      clinicalRecommendations: [
        "建议结合肝功能、血清学标志物、病毒载量和临床病程判断活动性感染与疾病阶段。",
        "建议核对 broad 大亚型与子亚型/基因型是否一致，必要时复核覆盖度和关键参考株选择。",
        "如涉及治疗或随访，应由临床结合指南、既往感染史和免疫状态综合决策。",
      ],
      cdcRecommendations: [
        "建议按 HAV/HBV/HCV/HDV/HEV 不同传播特点分别补充流调信息，不要混用同一处置路径。",
        "建议关注同一单位、家庭或共同暴露场景中的同型别聚集信号。",
        "如用于暴发研判，应补充采样时间、暴露史和更高分辨率的序列比较证据。",
      ],
    }, meta?.knowledgeTemplate);
  }
  if (isNaturalFocus) {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "natural-focus",
      clinicalRisk: "中风险",
      cdcRisk: "中",
      evidenceBasis: meta?.isBandavirus
        ? "Bandavirus 大亚型、A_F/CJ 三片段分型、重配提示与覆盖度"
        : "Orthohantavirus broad 分型、S 片段定型、L/M/S 片段参考选择与覆盖度",
      clinicalMeaning: knowledgeText || "该类病毒结果应结合发热、出血倾向、肾损伤或血小板变化等临床表现判断；分段分型结果可辅助判断自然疫源背景和潜在暴露来源。",
      cdcMeaning: "该类病毒具备自然疫源性监测意义，建议结合病例暴露史、地域来源、媒介/宿主线索和本地监测资料综合研判。",
      clinicalRecommendations: [
        "建议结合发热、血小板、肾功能、出血倾向及流行病学暴露史综合判断临床相关性。",
        "建议复核分段分型是否一致，若 L/M/S 或 A_F/CJ 证据不一致，应谨慎解释潜在重配或混合信号。",
        "如病情进展或暴露史明确，建议结合规范检测和临床专科意见动态评估。",
      ],
      cdcRecommendations: [
        "建议补充病例居住地、活动地、野外/农田暴露、动物或媒介接触信息。",
        "建议结合宿主和媒介监测资料判断是否存在自然疫源地活跃信号。",
        "如同一区域出现多例同型别结果，建议开展时空聚集和传播风险复核。",
      ],
    }, meta?.knowledgeTemplate);
  }
  if (isVectorBorne) {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "vector-borne",
      clinicalRisk: "低-中风险",
      cdcRisk: "中",
      clinicalMeaning: knowledgeText || "虫媒病毒检出需结合发热、皮疹、关节痛、出血表现、旅行史和采样时间窗判断临床相关性；分型结果更适合用于输入来源和传播背景分析。",
      cdcMeaning: "该结果具备虫媒病毒监测意义，建议结合旅行史、媒介密度、病例时空分布和本地输入/本地传播背景研判。",
      clinicalRecommendations: [
        "建议结合旅行史、蚊媒暴露史、发病时间窗和血清学/核酸复核结果判断临床意义。",
        "建议关注采样时间对核酸检出率的影响，必要时补充血清学或复采证据。",
        "如出现重症表现，应结合当地诊疗规范和实验室确认结果及时评估。",
      ],
      cdcRecommendations: [
        "建议核查旅行史、活动轨迹和发病地，区分输入病例与本地传播风险。",
        "建议结合媒介密度、季节和周边病例分布判断是否需要强化媒介控制。",
        "如出现同区域同时间多例，应补充分型/系统发育比较和现场流调证据。",
      ],
    }, meta?.knowledgeTemplate);
  }
  if (isBloodborne) {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "bloodborne",
      clinicalRisk: "中风险",
      cdcRisk: "中",
      clinicalMeaning: knowledgeText || "血源性或慢性感染相关病毒结果应结合确认试验、病毒载量、免疫状态和既往诊疗史综合解释；测序分型可辅助耐药、传播背景和随访管理。",
      cdcMeaning: "该结果可用于传播网络和重点人群监测线索，但不应替代确认试验、个案管理和规范报告流程。",
      clinicalRecommendations: [
        "建议结合确认试验、病毒载量、免疫状态和既往治疗史判断临床意义。",
        "如涉及耐药或亚型解释，应复核覆盖度、关键位点和参考选择结果。",
        "诊疗决策应依据现行临床指南和专科评估，不以单一测序报告替代。",
      ],
      cdcRecommendations: [
        "建议结合个案管理、传播风险评估和重点人群监测资料进行解释。",
        "如用于传播网络分析，应补充匿名化流调信息和更高分辨率序列比较。",
        "涉及报告管理时，应按现行规范和确认试验结果执行。",
      ],
    }, meta?.knowledgeTemplate);
  }
  if (isRespiratory) {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "respiratory",
      clinicalMeaning: knowledgeText || "呼吸道病毒检出需结合症状、采样部位、病程阶段和共感染背景判断临床意义；分型结果可辅助判断流行株背景和院内/社区传播线索。",
      cdcMeaning: "呼吸道病毒结果适合纳入季节性和聚集性监测，重点关注同型别病例聚集、特殊机构暴发和重症病例。",
      clinicalRecommendations: [
        "建议结合呼吸道症状、病程阶段、采样质量和共感染证据判断临床相关性。",
        "如为重症或免疫低下患者，建议复核覆盖度和关键变异位点。",
        "抗病毒或感染控制决策应结合当地诊疗规范、病原确认和患者风险因素。",
      ],
      cdcRecommendations: [
        "建议关注同型别病例在学校、养老机构、医院等场景中的聚集信号。",
        "建议与同期本地呼吸道病原监测数据对照，判断是否存在流行株变化。",
        "如出现异常重症或聚集性事件，应补充分型和同源性证据。",
      ],
    }, meta?.knowledgeTemplate);
  }
  if (isEnteric) {
    return applyKnowledgeVirusReportTemplate({
      ...base,
      category: "enteric",
      clinicalMeaning: knowledgeText || "肠道病毒或胃肠炎相关病毒检出需结合腹泻、呕吐、发热、采样时间和暴露史判断临床相关性；分型结果可辅助食品、水源或机构聚集事件研判。",
      cdcMeaning: "该结果适合纳入肠道传染病或胃肠炎聚集性监测，重点关注共同暴露、机构传播和同型别聚集。",
      clinicalRecommendations: [
        "建议结合胃肠道症状、采样时间、脱水程度和共感染结果判断临床意义。",
        "如结果用于个案诊疗，需结合病程和其他病原检测排除偶然携带或残留核酸。",
        "建议复核分型和覆盖度，尤其是用于聚集性事件解释时。",
      ],
      cdcRecommendations: [
        "建议补充共同就餐、水源、托幼/学校/养老机构暴露史。",
        "建议关注同型别病例聚集，并结合环境或食品样本结果进行综合判断。",
        "如涉及暴发调查，应补充采样时间轴和序列同源性比较。",
      ],
    }, meta?.knowledgeTemplate);
  }
  return applyKnowledgeVirusReportTemplate(base, meta?.knowledgeTemplate);
}

function extractVirusScenarioInterpretation(data) {
  const task = data?.task || {};
  const serotype = data?.sections?.serotype || {};
  const coverage = data?.sections?.assembly?.coverage || {};
  const knowledgeTemplate = serotype?.report_template && typeof serotype.report_template === "object" ? serotype.report_template : {};
  if (isInfluenzaTypingReport(data)) {
    const meta = extractInfluenzaResearchInterpretation(data);
    const typingResult = formatVirusTypingResult({ ...meta, kind: "influenza" });
    const template = resolveVirusReportTemplate({ ...meta, kind: "influenza", typingResult, knowledgeTemplate });
    return {
      ...meta,
      ...template,
      kind: "influenza",
      virusLabel: "流感病毒",
      typingResult,
      taskName: String(task?.name || task?.id || "-"),
      sampleType: String(task?.sample_type || task?.inputtype || "-"),
      submittingUnit: String(task?.owner_group || task?.owner || "-"),
    };
  }
  if (isMonkeypoxNextcladeReport(data)) {
    const meta = extractMonkeypoxResearchInterpretation(data);
    const typingResult = formatVirusTypingResult({ ...meta, kind: "monkeypox" });
    const template = resolveVirusReportTemplate({ ...meta, kind: "monkeypox", typingResult, knowledgeTemplate });
    return {
      ...meta,
      ...template,
      kind: "monkeypox",
      virusLabel: "猴痘病毒",
      typingResult,
      taskName: String(task?.name || task?.id || "-"),
      sampleType: String(task?.sample_type || task?.inputtype || "-"),
      submittingUnit: String(task?.owner_group || task?.owner || "-"),
    };
  }
  const meta = extractRespiratoryNextcladeResearchInterpretation(data);
  const typingResult = formatVirusTypingResult(meta);
  const knowledgeText = normalizeKnowledgeNarrativeText(meta.knowledgeHeadline || meta.knowledgeFragments?.[0] || "");
  const template = resolveVirusReportTemplate({ ...meta, kind: "generic-virus", typingResult, knowledgeTemplate: meta.knowledgeTemplate || knowledgeTemplate }, knowledgeText);
  return {
    ...meta,
    ...template,
    kind: "generic-virus",
    typingResult,
    taskName: String(task?.name || task?.id || "-"),
    sampleType: String(task?.sample_type || task?.inputtype || "-"),
    submittingUnit: String(task?.owner_group || task?.owner || "-"),
    meanDepth: String(coverage?.mean_depth ?? meta.meanDepth ?? "--"),
  };
}

function buildVirusClinicalScene(data) {
  const meta = extractVirusScenarioInterpretation(data);
  const reportDate = formatDateTime(data?.task?.updated_at || data?.task?.finished_at || data?.task?.created_at || Date.now());
  const riskClass = meta.clinicalRisk === "中风险" ? "level-mid" : meta.clinicalRisk === "低-中风险" ? "level-low" : meta.clinicalRisk === "低风险" ? "level-low" : "level-mid";
  const mutationSummary = Number(meta.mutationCount || 0) > 0
    ? `当前整理 ${meta.mutationCount} 个突变位点，高质量 ${meta.highMutationCount ?? "--"} 个、低质量 ${meta.lowMutationCount ?? "--"} 个${Array.isArray(meta.topMutations) && meta.topMutations.length ? `；代表性位点包括 ${meta.topMutations.join("；")}` : ""}。`
    : "当前未整理出可用于临床摘要的重点突变位点。";
  const recommendations = Array.isArray(meta.clinicalRecommendations) && meta.clinicalRecommendations.length
    ? meta.clinicalRecommendations
    : [
        "建议结合样本类型、采样时间、症状体征、影像学/实验室指标及必要的复核检测综合判断临床相关性。",
        "若用于临床处置，本报告提供的是病原分型和分子证据，不替代法定诊断标准、病原学确认或临床医嘱。",
        "如结果与临床表现不一致，建议复核原始 reads、覆盖度和关键位点，并结合其他检测方法确认。",
      ];
  return `
    <article class="scene-report-card clinical-scene-card clinical-document">
      <header class="clinical-document-header">
        <div class="clinical-document-kicker">Viral Clinical Report</div>
        <div class="clinical-document-headline">
          <div class="clinical-headline-copy">
            <h3 class="clinical-document-title">病毒检测医院报告</h3>
            <p class="clinical-document-subtitle">面向临床阅读，突出病毒检出结论、分型证据、临床意义和复核建议。</p>
          </div>
          <div class="clinical-risk-stamp ${riskClass}">
            <span class="clinical-risk-stamp-label">临床关注</span>
            <strong>${escapeHtml(meta.clinicalRisk)}</strong>
          </div>
        </div>
        <section class="clinical-report-head">
          <div class="clinical-report-head-item"><span>报告编号</span><strong>${escapeHtml(meta.taskName || data?.task?.name || data?.task?.id || "-")}</strong></div>
          <div class="clinical-report-head-item"><span>送检单位</span><strong>${escapeHtml(meta.submittingUnit || data?.task?.owner_group || data?.task?.owner || "-")}</strong></div>
          <div class="clinical-report-head-item"><span>样本类型</span><strong>${escapeHtml(meta.sampleType || data?.task?.sample_type || "-")}</strong></div>
          <div class="clinical-report-head-item"><span>主要检出病毒</span><strong>${escapeHtml(meta.speciesName || meta.virusLabel)}</strong></div>
          <div class="clinical-report-head-item"><span>分型结果</span><strong>${escapeHtml(meta.typingResult)}</strong></div>
          <div class="clinical-report-head-item"><span>报告日期</span><strong>${escapeHtml(reportDate)}</strong></div>
        </section>
      </header>
      <section class="clinical-impression-band">
        <div class="clinical-impression-main">
          <span class="clinical-impression-label">临床结论摘要</span>
          <p>本次病毒测序结果提示 <strong>${escapeHtml(meta.speciesName || meta.virusLabel)}</strong> 相关分子证据，核心分型为 <strong>${escapeHtml(meta.typingResult)}</strong>。</p>
        </div>
        <div class="clinical-impression-tags">
          <span>检测对象：${escapeHtml(meta.virusLabel || "病毒")}</span>
          <span>覆盖深度：${escapeHtml(meta.meanDepth || "--")}</span>
          <span>证据属性：分子分型</span>
        </div>
      </section>
      <div class="clinical-report-main">
        <section class="clinical-report-section">
          <div class="clinical-report-section-title"><span class="num">一</span><div><h4>病毒检出结果</h4><p>主要病毒及分型结论。</p></div></div>
          <div class="clinical-report-section-body">
            <div class="clinical-keyline"><strong>${escapeHtml(meta.speciesName || meta.virusLabel)}</strong><span class="${riskClass}">${escapeHtml(meta.typingResult)}</span></div>
            <p>分类学背景：${escapeHtml(meta.familyGenus || "--")}。</p>
            <p>主要证据：${escapeHtml(meta.evidenceBasis)}。</p>
          </div>
        </section>
        <section class="clinical-report-section">
          <div class="clinical-report-section-title"><span class="num">二</span><div><h4>临床意义</h4><p>结果边界和解释口径。</p></div></div>
          <div class="clinical-report-section-body">
            <p>${escapeHtml(meta.clinicalMeaning)}</p>
            <p>${escapeHtml(mutationSummary)}</p>
          </div>
        </section>
        <section class="clinical-report-section clinical-report-section-emphasis">
          <div class="clinical-report-section-title"><span class="num">三</span><div><h4>复核与处置建议</h4><p>原则性建议，不替代临床诊疗。</p></div></div>
          <div class="clinical-report-section-body">
            <ul class="clinical-guidance-list">
              ${recommendations.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
            </ul>
          </div>
        </section>
      </div>
    </article>
  `;
}

function buildVirusCdcScene(data) {
  const meta = extractVirusScenarioInterpretation(data);
  const riskLabel = meta.cdcRisk || "待评估";
  const riskCopy = riskLabel === "中"
    ? "提示存在需要持续监测的传播或输入风险，建议结合流调资料动态评估。"
    : riskLabel === "低"
      ? "当前分子结果更适合作为常规监测线索，暂未单独提示高优先级公共卫生处置。"
      : "当前分型或覆盖度证据仍需复核，暂不宜单独据此作出公共卫生风险定级。";
  const highQualityMutationText = Number(meta.mutationCount || 0) > 0
    ? `突变位点共 ${meta.mutationCount} 个，高质量 ${meta.highMutationCount ?? "--"} 个；${Array.isArray(meta.topMutations) && meta.topMutations.length ? `代表性位点：${meta.topMutations.join("；")}。` : ""}`
    : "当前未形成可用于监测摘要的突变位点列表。";
  const recommendations = Array.isArray(meta.cdcRecommendations) && meta.cdcRecommendations.length
    ? meta.cdcRecommendations
    : [
        "建议结合病例发病时间、采样地点、旅行/暴露史和密切接触信息，判断输入性、聚集性或本地传播可能。",
        "建议将分型结果与本地历史序列或上级平台数据进行比对，必要时补充系统发育或同源性分析。",
        "若涉及法定报告、重点传染病或聚集性事件，应按现行规范及时上报并补充流调证据。",
      ];
  return `
    <article class="scene-report-card cdc-scene-card cdc-document">
      <header class="cdc-document-header">
        <div>
          <div class="cdc-document-kicker">Viral Public Health Report</div>
          <h3>病毒疾控监测与风险评估报告</h3>
          <p>面向疾控与流行病学场景，突出病毒分型、传播监测价值、证据边界和后续处置建议。</p>
        </div>
        <div class="cdc-risk-badge risk-${escapeHtml(riskLabel)}">
          <span>监测风险</span>
          <strong>${escapeHtml(riskLabel)}</strong>
        </div>
      </header>
      <section class="cdc-document-grid">
        <section>
          <h4>一、检测结果概述</h4>
          <p>本次样本主要检出 <strong>${escapeHtml(meta.speciesName || meta.virusLabel)}</strong>，分型/谱系结果为 <strong>${escapeHtml(meta.typingResult)}</strong>。</p>
          <p>分类学背景：<strong>${escapeHtml(meta.familyGenus || "--")}</strong>；证据基础：<strong>${escapeHtml(meta.evidenceBasis)}</strong>。</p>
        </section>
        <section>
          <h4>二、监测意义</h4>
          <p>${escapeHtml(meta.cdcMeaning)}</p>
          <p>${escapeHtml(highQualityMutationText)}</p>
        </section>
        <section class="cdc-section-wide">
          <h4>三、风险等级判定</h4>
          <p>综合判定当前病毒监测风险为 <strong>${escapeHtml(riskLabel)}</strong>。${escapeHtml(riskCopy)}</p>
        </section>
        <section class="cdc-section-wide">
          <h4>四、流行病学处置建议</h4>
          <ul class="cdc-guidance-list">
            ${recommendations.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
          </ul>
        </section>
      </section>
    </article>
  `;
}

function buildResearchScene(data) {
  if (isInfluenzaTypingReport(data)) {
    const meta = extractInfluenzaResearchInterpretation(data);
    const attentionMarkup = renderResearchAttentionStack([
      {
        title: "重点提示",
        tone: inferResearchAttentionTone(meta?.qualityLabel),
        body: `当前样本核心分型结论为 ${meta.subtypeCall}，HA 亚型为 ${meta.haSubtype}，NA 亚型为 ${meta.naSubtype}。`,
      },
    ]);
    return `
      <article class="scene-report-card research-scene-card">
        <header class="scene-report-head">
          <span class="scene-report-kicker">Research Interpretation</span>
          <h3>科研判读摘要</h3>
        </header>
        <div class="scene-report-sections">
          <section>
            <h4>主导研究对象</h4>
            <p>当前样本主导病毒为 <strong>${escapeHtml(meta.speciesName)}</strong>，分型结果为 <strong>${escapeHtml(meta.subtypeCall)}</strong>，属于 <strong>${escapeHtml(meta.influenzaType)}</strong>。</p>
            <p>分类学背景归属于 <strong>${escapeHtml(meta.familyGenus)}</strong>，当前结果主要反映流感分型、segment 组成及变异谱特征，不适宜套用细菌 MLST 或血清群框架进行解释。</p>
          </section>
          <section>
            <h4>分子特征提示</h4>
            ${attentionMarkup}
            <p>本次结果支持 <strong>${escapeHtml(meta.subtypeCall)}</strong> 作为当前样本的核心分型结论，其中 HA 亚型为 <strong>${escapeHtml(meta.haSubtype)}</strong>，NA 亚型为 <strong>${escapeHtml(meta.naSubtype)}</strong>。</p>
            <p>最终参考集合由 ${escapeHtml(String(meta.segmentCount))} 个 segment 组成${meta.subtypeSegments.length ? `，其中重点 subtype 片段包括 ${escapeHtml(meta.subtypeSegments.join("；"))}` : ""}。</p>
            <p>覆盖度概览显示整体平均深度约为 <strong>${escapeHtml(meta.meanDepth)}</strong>，原始覆盖度约 <strong>${escapeHtml(meta.coverageFraction)}</strong>，10x 覆盖约 <strong>${escapeHtml(meta.coverage10x)}</strong>，100x 覆盖约 <strong>${escapeHtml(meta.coverage100x)}</strong>。</p>
            <p>${meta.mutationCount ? `当前共识别到 ${escapeHtml(String(meta.mutationCount))} 个突变位点，代表性位点包括 ${escapeHtml(meta.topMutations.join("；"))}。` : "当前未整理出可用于科研摘要的代表性突变位点。"} </p>
            <p>${renderKnowledgeNarrative(meta, "流感")}</p>
          </section>
          <section>
            <h4>证据边界</h4>
            <p>上述结论主要建立在 IRMA 风格参考集合初筛、HA/NA 最优亚型选择、segment 覆盖情况和突变位点表之上，因此其证据属性更偏向流感分型及分子流行病学归属。</p>
            <p>对于毒力、宿主适应性或传播背景的进一步解释，仍需结合系统发育位置、关键氨基酸位点注释及历史株背景综合判断。</p>
          </section>
          <section>
            <h4>科研提示</h4>
            <p>${meta.segmentCoverageHints.length ? `深度较高的 segment 主要包括：${escapeHtml(meta.segmentCoverageHints.join("；"))}。` : "各 segment 的覆盖均一性可作为后续比较分析的重要参照。"} 在结果展示层面，分型结论、8 segment 组成、覆盖度分段特征及代表性突变位点构成核心信息。</p>
          </section>
        </div>
      </article>
    `;
  }
  if (isMonkeypoxNextcladeReport(data)) {
    const meta = extractMonkeypoxResearchInterpretation(data);
    const attentionMarkup = renderResearchAttentionStack([
      {
        title: "重点提示",
        tone: inferResearchAttentionTone(meta?.qualityLabel),
        body: `当前样本核心谱系结论为 ${meta.clade}${meta.lineage !== "--" ? ` / ${meta.lineage}` : ""}${meta.outbreak !== "--" ? `，outbreak 标记为 ${meta.outbreak}` : ""}。`,
      },
    ]);
    return `
      <article class="scene-report-card research-scene-card">
        <header class="scene-report-head">
          <span class="scene-report-kicker">Research Interpretation</span>
          <h3>科研判读摘要</h3>
        </header>
        <div class="scene-report-sections">
          <section>
            <h4>主导研究对象</h4>
            <p>当前样本主导病毒为 <strong>${escapeHtml(meta.speciesName)}</strong>，Nextclade 分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""}。</p>
            <p>分类学背景归属于 <strong>${escapeHtml(meta.familyGenus)}</strong>，当前结果更适合从 clade、lineage、覆盖度及变异谱角度进行解释，不适宜沿用细菌物种鉴定、MLST 或耐药基因框架。</p>
          </section>
          <section>
            <h4>分子特征提示</h4>
            ${attentionMarkup}
            <p>本次结果支持 <strong>${escapeHtml(meta.clade)}</strong>${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""} 作为当前样本的核心谱系结论${meta.outbreak !== "--" ? `，对应 outbreak 标记为 <strong>${escapeHtml(meta.outbreak)}</strong>` : ""}。</p>
            <p>覆盖度概览显示整体平均深度约为 <strong>${escapeHtml(meta.meanDepth)}</strong>，1x 覆盖约 <strong>${escapeHtml(meta.coverageFraction)}</strong>，10x 覆盖约 <strong>${escapeHtml(meta.coverage10x)}</strong>，100x 覆盖约 <strong>${escapeHtml(meta.coverage100x)}</strong>。</p>
            <p>${meta.mutationCount ? `当前共整理出 ${escapeHtml(String(meta.mutationCount))} 个突变位点，其中高质量 ${escapeHtml(String(meta.highMutationCount))} 个、低质量 ${escapeHtml(String(meta.lowMutationCount))} 个${meta.topMutations.length ? `；代表性高质量位点包括 ${escapeHtml(meta.topMutations.join("；"))}` : ""}。` : "当前未整理出可用于科研摘要的猴痘突变位点。"} </p>
            <p>${renderKnowledgeNarrative(meta, "猴痘病毒")}</p>
          </section>
          <section>
            <h4>证据边界</h4>
            <p>上述结论主要建立在 hMPXV Nextclade 数据集的 clade/lineage 判定、比对覆盖度及 SnpEff 位点注释之上，因此其证据属性更偏向猴痘分子流行病学及谱系归属。</p>
            <p>${meta.qualityLabel !== "--" ? `本次 Nextclade 质控状态为 <strong>${escapeHtml(meta.qualityLabel)}</strong>，` : ""}对于传播背景、跨地区输入关系或 outbreak 聚类的进一步解释，仍需结合系统发育树、时间地点信息及关键位点原始比对证据综合判断。</p>
          </section>
          <section>
            <h4>科研提示</h4>
            <p>在结果展示层面，clade/lineage 判定、全基因组覆盖度、突变位点表及 IGV 关键位点核查结果构成主要证据；若涉及 Nanopore 等测序技术背景，homopolymer 区域 indel 的低质量处理方式亦应同步说明。</p>
          </section>
        </div>
      </article>
    `;
  }
  if (isChikvNextcladeReport(data)) {
    const meta = extractRespiratoryNextcladeResearchInterpretation(data);
    const attentionMarkup = renderResearchAttentionStack([
      {
        title: "重点提示",
        tone: inferResearchAttentionTone(meta?.qualityLabel),
        body: `当前样本核心谱系结论为 ${meta.clade}${meta.lineage !== "--" ? ` / ${meta.lineage}` : ""}。`,
      },
    ]);
    return `
      <article class="scene-report-card research-scene-card">
        <header class="scene-report-head">
          <span class="scene-report-kicker">Research Interpretation</span>
          <h3>科研判读摘要</h3>
        </header>
        <div class="scene-report-sections">
          <section>
            <h4>主导研究对象</h4>
            <p>当前样本主导病毒为 <strong>${escapeHtml(meta.speciesName)}</strong>，Nextclade 分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""}。</p>
            <p>分类学背景归属于 <strong>${escapeHtml(meta.familyGenus)}</strong>，当前结果主要围绕 CHIKV 的 clade、lineage、覆盖度、突变谱与参考株比对证据组织，不适宜套用细菌 MLST、血清群或肝炎病毒 broad/subtype 双层模板进行解释。</p>
          </section>
          <section>
            <h4>分子特征提示</h4>
            ${attentionMarkup}
            <p>本次结果支持 <strong>${escapeHtml(meta.clade)}</strong>${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""} 作为当前样本的核心谱系结论，可作为后续流行病学背景和参考株亲缘位置解释的主要起点。</p>
            <p>覆盖度概览显示整体平均深度约为 <strong>${escapeHtml(meta.meanDepth)}</strong>，1x 覆盖约 <strong>${escapeHtml(meta.coverageFraction)}</strong>，10x 覆盖约 <strong>${escapeHtml(meta.coverage10x)}</strong>，100x 覆盖约 <strong>${escapeHtml(meta.coverage100x)}</strong>。</p>
            <p>${meta.mutationCount ? `当前共整理出 ${escapeHtml(String(meta.mutationCount))} 个突变位点，其中高质量 ${escapeHtml(String(meta.highMutationCount))} 个、低质量 ${escapeHtml(String(meta.lowMutationCount))} 个${meta.topMutations.length ? `；代表性高质量位点包括 ${escapeHtml(meta.topMutations.join("；"))}` : ""}。` : "当前未整理出可用于科研摘要的 CHIKV 突变位点。"} </p>
            <p>${renderKnowledgeNarrative(meta, "基孔肯雅病毒")}</p>
          </section>
          <section>
            <h4>证据边界</h4>
            <p>上述结论主要建立在 CHIKV Nextclade 数据集分型、固定参考比对、覆盖度统计、突变位点表与系统发育树结果之上，因此其证据属性更偏向 CHIKV 的谱系归属与分子流行病学解释。</p>
            <p>${meta.qualityLabel !== "--" ? `当前 Nextclade 质控状态为 <strong>${escapeHtml(meta.qualityLabel)}</strong>。` : ""}对于输入来源、跨地区传播链、暴发关联或功能影响的进一步解释，仍需结合时间地点信息、系统树位置与关键位点原始比对证据综合判断。</p>
          </section>
          <section>
            <h4>科研提示</h4>
            <p>在结果展示层面，CHIKV 分型总表、质量指标、突变位点表、IGV 关键位点核查结果及系统发育树共同构成主要证据；其中系统树和 IGV 更适合作为疾控场景下复核谱系归属与关键位点可靠性的第二视角。</p>
          </section>
        </div>
      </article>
    `;
  }
  if (isAstroviridaeTypingReport(data)) {
    const meta = extractRespiratoryNextcladeResearchInterpretation(data);
    const attentionMarkup = renderResearchAttentionStack([
      {
        title: "重点提示",
        tone: inferResearchAttentionTone(meta?.qualityLabel),
        body: `当前样本核心分型结论为 ${meta.clade}${meta.group !== "--" ? ` / ${meta.group}` : ""}${meta.lineage !== "--" ? ` / ${meta.lineage}` : ""}。`,
      },
    ]);
    return `
      <article class="scene-report-card research-scene-card">
        <header class="scene-report-head">
          <span class="scene-report-kicker">Research Interpretation</span>
          <h3>科研判读摘要</h3>
        </header>
        <div class="scene-report-sections">
          <section>
            <h4>主导研究对象</h4>
            <p>当前样本主导病毒为 <strong>${escapeHtml(meta.speciesName)}</strong>，ORF2 分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>${meta.group !== "--" ? `，病毒属归为 <strong>${escapeHtml(meta.group)}</strong>` : ""}${meta.lineage !== "--" ? `，病毒种为 <strong>${escapeHtml(meta.lineage)}</strong>` : ""}。</p>
            <p>分类学背景归属于 <strong>${escapeHtml(meta.familyGenus)}</strong>，当前结果主要反映 ORF2 分型、属种归属、最优参考选择以及 ORF2 系统树支持情况，不适宜沿用班达病毒的大亚型或分段重配框架进行解释。</p>
          </section>
          <section>
            <h4>分子特征提示</h4>
            ${attentionMarkup}
            <p>本次结果支持 <strong>${escapeHtml(meta.clade)}</strong>${meta.group !== "--" ? ` / <strong>${escapeHtml(meta.group)}</strong>` : ""}${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""} 作为当前样本的核心星状病毒分型与分类归属结论。</p>
            <p>覆盖度概览显示整体平均深度约为 <strong>${escapeHtml(meta.meanDepth)}</strong>，1x 覆盖约 <strong>${escapeHtml(meta.coverageFraction)}</strong>，10x 覆盖约 <strong>${escapeHtml(meta.coverage10x)}</strong>，100x 覆盖约 <strong>${escapeHtml(meta.coverage100x)}</strong>。</p>
            <p>结合最优全基因组参考与 VADR 注释，当前结果可对 ORF2 所在结构区段进行交叉核对，并通过 ORF2 系统发育树对初筛分型结论作进一步验证。</p>
            <p>${meta.mutationCount ? `当前共整理出 ${escapeHtml(String(meta.mutationCount))} 个突变位点，其中高质量 ${escapeHtml(String(meta.highMutationCount))} 个、低质量 ${escapeHtml(String(meta.lowMutationCount))} 个${meta.topMutations.length ? `；代表性高质量位点包括 ${escapeHtml(meta.topMutations.join("；"))}` : ""}。` : "当前未整理出可用于科研摘要的星状病毒突变位点。"} </p>
            <p>${renderKnowledgeNarrative(meta, meta.virusLabel)}</p>
          </section>
          <section>
            <h4>证据边界</h4>
            <p>上述结论主要建立在 ORF2 参考库分型、候选参考竞争、VADR 注释、ORF2 系统树及覆盖度统计之上，因此其证据属性更偏向 ORF2 分类归属及分子流行病学解释。</p>
            <p>对于宿主适应性、传播链或功能影响的进一步解释，仍需结合完整基因组背景、关键氨基酸位点及历史参考株资料综合判断。</p>
          </section>
          <section>
            <h4>科研提示</h4>
            <p>在结果展示层面，ORF2 分型、病毒属/种判定、最优参考株、覆盖度概览、ORF2 系统发育树及 IGV 关键位点证据构成主要信息。</p>
          </section>
        </div>
      </article>
    `;
  }
  if (isHivTypingReport(data)) {
    const serotype = data?.sections?.serotype || {};
    const hivKnowledge = extractViralKnowledgeSummary(serotype);
    const summaryCards = Array.isArray(serotype?.summary_cards) ? serotype.summary_cards : [];
    const mutationPanels = Array.isArray(serotype?.mutation_panels) ? serotype.mutation_panels : [];
    const resistanceTable = serotype?.resistance_table && typeof serotype.resistance_table === "object"
      ? serotype.resistance_table
      : { rows: [], columns: [] };
    const resistanceRows = Array.isArray(resistanceTable?.rows) ? resistanceTable.rows : [];
    const subtype = String(serotype?.predicted_clade || "--").trim() || "--";
    const broadType = String(serotype?.predicted_group || "--").trim() || "--";
    const recombination = String(summaryCards.find((item) => item?.label === "重组判定")?.value || "--").trim() || "--";
    const representativeReference = String(summaryCards.find((item) => item?.label === "代表株参考")?.value || "--").trim() || "--";
    const candidateParents = String(mutationPanels.find((item) => item?.label === "候选父本")?.value || "--").trim() || "--";
    const resistanceClasses = ["NRTI 最高等级", "NNRTI 最高等级", "PI 最高等级", "INSTI 最高等级"]
      .map((label) => {
        const value = String(summaryCards.find((item) => item?.label === label)?.value || "--").trim() || "--";
        return `${label.replace(" 最高等级", "")}: ${value}`;
      });
    const attentionMarkup = renderResearchAttentionStack([
      {
        title: "重点提示",
        tone: inferResearchAttentionTone(resistanceClasses.join(" ")),
        body: `当前样本核心分型结论为 ${subtype !== "--" ? subtype : broadType}，重组判定为 ${recombination}。`,
      },
    ]);
    return `
      <article class="scene-report-card research-scene-card">
        <header class="scene-report-head">
          <span class="scene-report-kicker">Research Interpretation</span>
          <h3>科研判读摘要</h3>
        </header>
        <div class="scene-report-sections">
          <section>
            <h4>主导研究对象</h4>
            <p>当前样本主导病毒为 <strong>Human immunodeficiency virus</strong>，broad 分型结果指向 <strong>${escapeHtml(broadType)}</strong>${subtype !== "--" && subtype !== broadType ? `，进一步收敛到 <strong>${escapeHtml(subtype)}</strong> 子亚型背景` : ""}。</p>
            <p>分类学背景归属于 <strong>Retroviridae / Lentivirus</strong>，当前页面主要整合 HIV-1/HIV-2 broad typing、代表株覆盖度筛选、REGA-like 子亚型/重组证据和 HIVDB 耐药解释。</p>
            <p>${renderKnowledgeNarrative({
              knowledgeHeadline: hivKnowledge.headline,
              knowledgeFragments: hivKnowledge.fragments,
            }, "HIV 分型")}</p>
          </section>
          <section>
            <h4>分子特征提示</h4>
            ${attentionMarkup}
            <p>本次结果支持 <strong>${escapeHtml(subtype !== "--" ? subtype : broadType)}</strong> 作为当前样本的核心分型结论，重组判定为 <strong>${escapeHtml(recombination)}</strong>${representativeReference !== "--" ? `，一致性生成阶段采用的代表株参考为 <strong>${escapeHtml(representativeReference)}</strong>` : ""}。</p>
            <p>${candidateParents !== "--" ? `候选父本组成提示为 ${escapeHtml(candidateParents)}。` : "当前未整理出明确的候选父本摘要。"} 这一层证据更偏向 HIV 分子流行病学与重组背景判读。</p>
            <p>耐药解释层面，当前各药物类别的最高等级分别为 ${escapeHtml(resistanceClasses.join("；"))}。</p>
            <p>${resistanceRows.length ? `当前共输出 ${escapeHtml(String(resistanceRows.length))} 条 HIVDB 药物解释记录；建议优先结合 PI / NNRTI / INSTI 的高等级条目和对应触发规则进行判读。` : "当前未输出 HIVDB 药物解释记录。"} </p>
          </section>
          <section>
            <h4>证据边界</h4>
            <p>上述结论建立在 broad 参考覆盖度筛选、子亚型代表株竞争、bootscan 窗口支持和 HIVDB XML 打分规则之上，因此其证据属性主要是 HIV 型别归属、重组提示与药物耐药风险。</p>
            <p>对于传播链、感染来源或复杂重组事件的进一步解释，仍需结合更完整的参考集合、系统树位置以及关键断点附近的原始比对证据综合判断。</p>
          </section>
          <section>
            <h4>科研提示</h4>
            <p>在结果展示层面，HIV 大亚型、子亚型、重组判定、代表株筛选、bootscan 曲线、PR/RT/IN 突变与 HIVDB 药物解释共同构成主要证据链。</p>
          </section>
        </div>
      </article>
    `;
  }
  if (isRsvNextcladeReport(data) || isHmpvNextcladeReport(data) || isDenvNextcladeReport(data) || isZikavNextcladeReport(data) || isChikvNextcladeReport(data) || isEbolaNextcladeReport(data) || isHpivTypingReport(data) || isHadvTypingReport(data) || isNorovirusTypingReport(data) || isEnterovirusTypingReport(data) || isHepatovirusTypingReport(data) || isBandavirusTypingReport(data) || isOrthohantavirusTypingReport(data) || isOrthoebolavirusTypingReport(data) || isAstroviridaeTypingReport(data) || isRhinovirusTypingReport(data) || isSeasonalHcovTypingReport(data) || isRotavirusTypingReport(data) || isHivTypingReport(data)) {
    const meta = extractRespiratoryNextcladeResearchInterpretation(data);
    const attentionMarkup = renderResearchAttentionStack([
      {
        title: "重点提示",
        tone: inferResearchAttentionTone(meta?.qualityLabel),
        body: `${meta.isHepatovirus
          ? `当前样本核心分型结论为 ${meta.group !== "--" ? meta.group : "肝炎病毒"}${meta.clade !== "--" && meta.clade !== meta.group ? ` / ${meta.clade}` : ""}。`
          : `当前样本核心分型结论为 ${meta.clade}${meta.isRotavirus ? (meta.lineage !== "--" ? ` / ${meta.lineage}` : "") : (meta.isAstroviridae ? `${meta.group !== "--" ? ` / ${meta.group}` : ""}${meta.lineage !== "--" ? ` / ${meta.lineage}` : ""}` : (meta.isBandavirus ? `${meta.group !== "--" ? ` / ${meta.group}` : ""}${meta.lineage !== "--" ? ` / ${meta.lineage}` : ""}` : (meta.isOrthohantavirus ? `${meta.group !== "--" ? ` / ${meta.group}` : ""}` : (!meta.isHpiv && meta.lineage !== "--" ? ` / ${meta.lineage}` : ""))))}。`}`,
      },
    ]);
    return `
      <article class="scene-report-card research-scene-card">
        <header class="scene-report-head">
          <span class="scene-report-kicker">Research Interpretation</span>
          <h3>科研判读摘要</h3>
        </header>
        <div class="scene-report-sections">
          <section>
            <h4>主导研究对象</h4>
            <p>当前样本主导病毒为 <strong>${escapeHtml(meta.speciesName)}</strong>，${meta.isRotavirus ? `大组分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>${meta.lineage !== "--" ? `，组合分型为 <strong>${escapeHtml(meta.lineage)}</strong>` : ""}。` : (meta.isHpiv ? `最优参考分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>。` : (meta.isHepatovirus ? `肝炎病毒大亚型结果指向 <strong>${escapeHtml(meta.group)}</strong>${meta.clade !== "--" && meta.clade !== meta.group ? `，子亚型/基因型为 <strong>${escapeHtml(meta.clade)}</strong>` : ""}。` : (meta.isAstroviridae ? `ORF2 分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>${meta.group !== "--" ? `，病毒属归为 <strong>${escapeHtml(meta.group)}</strong>` : ""}${meta.lineage !== "--" ? `，病毒种为 <strong>${escapeHtml(meta.lineage)}</strong>` : ""}。` : (meta.isBandavirus ? `Bandavirus 大亚型结果指向 <strong>${escapeHtml(meta.group)}</strong>${meta.clade !== "--" ? `，A_F(LMS) 为 <strong>${escapeHtml(meta.clade)}</strong>` : ""}${meta.lineage !== "--" ? `，CJ(LMS) 为 <strong>${escapeHtml(meta.lineage)}</strong>` : ""}。` : (meta.isOrthohantavirus ? `Orthohantavirus 分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>${meta.group !== "--" ? `，S 片段分型为 <strong>${escapeHtml(meta.group)}</strong>` : ""}。` : `Nextclade 分型结果指向 <strong>${escapeHtml(meta.clade)}</strong>${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""}。`)))))}</p>
            <p>分类学背景归属于 <strong>${escapeHtml(meta.familyGenus)}</strong>，当前结果主要反映 ${escapeHtml(meta.virusShort)} 的${meta.isRotavirus ? "大组判定、G/P 组合分型、覆盖度及参考选择特征" : (meta.isHpiv ? "参考亚型、覆盖度及变异谱特征" : (meta.isHepatovirus ? "大亚型判定、子亚型/基因型筛选、覆盖度及参考选择特征" : (meta.isAstroviridae ? "ORF2 分型、属种归属、最优参考及 ORF2 系统树特征" : (meta.isBandavirus ? "大亚型、A_F/CJ 三片段分型、重配提示及参考选择特征" : (meta.isOrthohantavirus ? "broad 分型、S 片段定型、L/M/S 三片段参考选择及自然疫源病解释背景" : "clade、lineage、覆盖度及变异谱特征")))))}。</p>
          </section>
          <section>
            <h4>分子特征提示</h4>
            ${attentionMarkup}
            <p>${meta.isHepatovirus
              ? `本次结果支持 <strong>${escapeHtml(meta.group !== "--" ? meta.group : "肝炎病毒")}</strong>${meta.clade !== "--" && meta.clade !== meta.group ? ` 大亚型背景下的 <strong>${escapeHtml(meta.clade)}</strong> 子亚型/基因型归属` : " 大亚型归属"}，当前样本更接近对应大亚型完整基因组参考集合。`
              : `本次结果支持 <strong>${escapeHtml(meta.clade)}</strong>${meta.isRotavirus ? (meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : "") : (meta.isAstroviridae ? `${meta.group !== "--" ? ` / <strong>${escapeHtml(meta.group)}</strong>` : ""}${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""}` : (meta.isBandavirus ? `${meta.group !== "--" ? ` / <strong>${escapeHtml(meta.group)}</strong>` : ""}${meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""}` : (meta.isOrthohantavirus ? `${meta.group !== "--" ? ` / <strong>${escapeHtml(meta.group)}</strong>` : ""}` : (!meta.isHpiv && meta.lineage !== "--" ? ` / <strong>${escapeHtml(meta.lineage)}</strong>` : ""))))} 作为当前样本的核心${meta.isRotavirus ? "轮状病毒分型" : (meta.isHpiv ? "参考分型" : (meta.isAstroviridae ? "星状病毒分型与分类归属" : (meta.isBandavirus ? "班达病毒分型与分段归属" : (meta.isOrthohantavirus ? "汉坦病毒分型与片段归属" : "谱系"))))}结论。`}</p>
            <p>覆盖度概览显示整体平均深度约为 <strong>${escapeHtml(meta.meanDepth)}</strong>，1x 覆盖约 <strong>${escapeHtml(meta.coverageFraction)}</strong>，10x 覆盖约 <strong>${escapeHtml(meta.coverage10x)}</strong>，100x 覆盖约 <strong>${escapeHtml(meta.coverage100x)}</strong>。</p>
            ${meta.isHepatovirus ? `<p>肝炎病毒结果采用“先 broad 大亚型、后对应大亚型子亚型/基因型”的两级筛选策略；当前结论已经过 broad typing 与对应参考库竞争两轮约束，更适合从分子流行病学和参考株亲缘背景角度解释，而不是沿用细菌 MLST 或血清群模板。</p>` : ""}
            ${meta.isAstroviridae ? `<p>结合最优全基因组参考与 VADR 注释，当前结果可对 ORF2 所在结构区段进行交叉核对，并通过 ORF2 系统发育树对初筛分型结论作进一步验证。</p>` : ""}
            ${meta.isBandavirus ? `<p>Bandavirus 结果可结合参考筛选结果、A_F 三片段最优参考与 CJ 三片段 genotype 共同评估分段一致性；若 L/M/S 判定不一致，则需进一步考虑潜在重组或重配信号。</p>` : ""}
            ${meta.isOrthohantavirus ? `<p>Orthohantavirus 结果可结合 broad 分型支持摘要、S 片段分型以及 L/M/S 三片段最优参考共同评估分型稳定性；若型别落在 HTNV、SEOV、DOBV、PUUV 或 AMRV 等条目时，应优先从 HFRS 背景解释其自然疫源病意义。</p>` : ""}
            <p>${meta.mutationCount ? `当前共整理出 ${escapeHtml(String(meta.mutationCount))} 个突变位点，其中高质量 ${escapeHtml(String(meta.highMutationCount))} 个、低质量 ${escapeHtml(String(meta.lowMutationCount))} 个${meta.topMutations.length ? `；代表性高质量位点包括 ${escapeHtml(meta.topMutations.join("；"))}` : ""}。` : `当前未整理出可用于科研摘要的 ${escapeHtml(meta.virusShort)} 突变位点。`}</p>
            <p>${renderKnowledgeNarrative(meta, meta.virusLabel)}</p>
          </section>
          <section>
            <h4>证据边界</h4>
            <p>上述结论主要建立在参考株比对、${meta.isRotavirus ? "大组覆盖度比较与 VP4/VP7 组合分型" : (meta.isHpiv ? "最优参考选择" : (meta.isHepatovirus ? "肝炎病毒 broad 大亚型筛选、子亚型/基因型参考竞争与 GFF 注释一致性验证" : (meta.isAstroviridae ? "ORF2 参考库分型、候选参考竞争、VADR 注释与 ORF2 系统树" : (meta.isBandavirus ? "Bandavirus 大亚型筛选、SFTSV 的 A_F/CJ 三片段分型及分段参考竞争" : (meta.isOrthohantavirus ? "Orthohantavirus broad 分型、S 片段定型、L/M/S 三片段参考竞争与知识库型别关联" : "Nextclade 数据集分型")))))}、覆盖度统计及突变位点表之上，因此其证据属性更偏向 ${escapeHtml(meta.virusShort)} 的${meta.isRotavirus ? "分型归属及分子流行病学" : (meta.isHpiv ? "参考分型及分子流行病学" : (meta.isHepatovirus ? "肝炎病毒型别归属及分子流行病学" : (meta.isAstroviridae ? "ORF2 分类归属及分子流行病学" : (meta.isBandavirus ? "片段分型归属及分子流行病学" : (meta.isOrthohantavirus ? "型别归属、分段证据整合及自然疫源病背景解释" : "谱系归属及分子流行病学")))))}。</p>
            <p>对于传播链、宿主适应或功能影响的进一步解释，仍需结合系统发育树、关键氨基酸位点及历史参考株背景综合判断。</p>
          </section>
          <section>
            <h4>科研提示</h4>
            <p>${!meta.isHpiv && !meta.isHepatovirus && !meta.isAstroviridae && !meta.isBandavirus && !meta.isOrthohantavirus && meta.qualityLabel !== "--" ? `当前 Nextclade 质控状态为 ${escapeHtml(meta.qualityLabel)}。` : ""}${meta.isHepatovirus ? "在结果展示层面，肝炎病毒大亚型、子亚型/基因型、最优参考株、覆盖度概览、突变表与 IGV 关键位点证据构成主要信息。" : (meta.isAstroviridae ? "在结果展示层面，ORF2 分型、病毒属/种判定、最优参考株、覆盖度概览、ORF2 系统发育树及 IGV 关键位点证据构成主要信息。" : (meta.isBandavirus ? "在结果展示层面，Bandavirus 大亚型、A_F(LMS)、CJ(LMS)、重配提示、三片段最优参考及覆盖度概览构成主要信息。" : (meta.isOrthohantavirus ? "在结果展示层面，Orthohantavirus 分型、S 片段分型、broad 支持摘要、L/M/S 三片段最优参考、知识库关联结论及 IGV 关键位点证据构成主要信息。" : "在结果展示层面，分型总表、覆盖度概览、高低质量突变表及 IGV 关键位点证据构成主要信息。")))}</p>
          </section>
        </div>
      </article>
    `;
  }
  const meta = extractResearchInterpretation(data);
  if (meta.isVirusLike) {
    const serotype = data?.sections?.serotype || {};
    const mode = resolveSerotypeMode(serotype);
    const summaryCards = Array.isArray(serotype?.summary_cards) ? serotype.summary_cards : [];
    const typingLabels = [];
    const predictedGroup = String(serotype?.predicted_group || "").trim();
    const predictedClade = String(serotype?.predicted_clade || "").trim();
    const predictedLineage = String(serotype?.predicted_lineage || "").trim();
    const referenceName = String(serotype?.reference_name || serotype?.selected_reference || "").trim();
    if (predictedGroup && predictedGroup !== "--") typingLabels.push(predictedGroup);
    if (predictedClade && predictedClade !== "--" && !typingLabels.includes(predictedClade)) typingLabels.push(predictedClade);
    if (predictedLineage && predictedLineage !== "--" && !typingLabels.includes(predictedLineage)) typingLabels.push(predictedLineage);
    if (!typingLabels.length) {
      summaryCards.forEach((item) => {
        const label = String(item?.label || "").trim();
        const value = String(item?.value || "").trim();
        if (!label || !value || value === "--") return;
        if (["大亚型", "子亚型", "HAV 子亚型", "HAV子亚型", "Lineage", "双位点分型", "VP1 分型", "HAdV 分型", "Orthohantavirus 分型", "大组分型", "组合分型", "ORF2 分型"].includes(label)) {
          typingLabels.push(value);
        }
      });
    }
    const supportSpeciesMarkup = meta.supportingSpecies.length
      ? `<p>背景共存物种：${escapeHtml(meta.supportingSpecies.map((item) => `${item.species}（${Number(item.ratio || 0).toFixed(2)}%）`).join("；"))}</p>`
      : `<p>当前未见占比达到重点阈值的次高异种背景，后续分析可优先聚焦主导病毒。</p>`;
    const intraspeciesMarkup = meta.intraspeciesSignals.length
      ? `<p>同种内分型/亚群信号：${escapeHtml(meta.intraspeciesSignals.map((item) => `${item.label}（${Number(item.ratio || 0).toFixed(2)}%）`).join("；"))}</p>`
      : "";
    const attentionMarkup = renderResearchAttentionStack([
      {
        title: "重点提示",
        tone: inferResearchAttentionTone(meta?.riskLevel),
        body: meta.conclusion,
      },
    ]);
    return `
      <article class="scene-report-card research-scene-card">
        <header class="scene-report-head">
          <span class="scene-report-kicker">Research Interpretation</span>
          <h3>科研判读摘要</h3>
        </header>
        <div class="scene-report-sections">
          <section>
            <h4>主导研究对象</h4>
            <p>当前样本主导病毒为 <strong>${escapeHtml(meta.speciesName)}</strong>（NCBI TaxID: ${escapeHtml(meta.taxid)}；学名：${escapeHtml(meta.scientificName)}）。知识库风险层级评估为 <strong>${escapeHtml(meta.riskLevel)}</strong>。</p>
            <p>${escapeHtml(meta.significance)}</p>
          </section>
          <section>
            <h4>分子特征提示</h4>
            ${attentionMarkup}
            <p>${escapeHtml(meta.conclusion)}</p>
            <p>${typingLabels.length ? `当前病毒分型结果支持：${escapeHtml(typingLabels.join(" / "))}。` : "当前未形成稳定的病毒分型标签，建议优先结合覆盖度、最优参考和关键变异位点综合判断。"}</p>
            <p>${referenceName && referenceName !== "--" ? `当前最优参考序列为 <strong>${escapeHtml(referenceName)}</strong>。` : "当前未提供稳定的最优参考序列名称。"}</p>
            <p>${renderKnowledgeNarrative({ knowledgeHeadline: serotype?.knowledge_summary?.headline, knowledgeFragments: Array.isArray(serotype?.knowledge_summary?.items) ? serotype.knowledge_summary.items.slice(0, 3).map((item) => {
              const sero = String(item?.serotype || item?.matched_on || "").trim() || "--";
              const interpretation = String(item?.interpretation || "").trim();
              return interpretation ? `${sero}：${interpretation}` : sero;
            }) : [] }, meta.speciesName)}</p>
            <p>${meta.eventSummaries.length ? `组合事件结果显示：${escapeHtml(meta.eventSummaries.join("；"))}。` : "当前未命中明确的病毒组合事件规则，结果更适合解释为分型、覆盖度与关键位点层面的分子特征。"} </p>
          </section>
          <section>
            <h4>证据边界</h4>
            <p>${escapeHtml(meta.confidenceHint)}</p>
            ${intraspeciesMarkup}
            ${supportSpeciesMarkup}
          </section>
          <section>
            <h4>科研提示</h4>
            <p>${escapeHtml(meta.nextSteps.join(" "))}</p>
          </section>
        </div>
      </article>
    `;
  }
  const neisseriaAmr = extractNeisseriaAmrInterpretation(data);
  const tbAmr = extractTbAmrInterpretation(data);
  const shouldShowNeisseriaAmr = isNeisseriaMeningitidisSpecies(meta.speciesName, meta.scientificName);
  const shouldShowTbAmr = String(meta.speciesName || "").includes("结核分枝杆菌")
    || String(meta.scientificName || "").toLowerCase().includes("mycobacterium tuberculosis");
  const dominantScientific = String(meta.scientificName || "").trim().toLowerCase();
  const dominantDisplay = String(meta.speciesName || "").trim().toLowerCase();
  const extraIntraspeciesSignals = [];
  const filteredSupportingSpecies = [];
  (Array.isArray(meta.supportingSpecies) ? meta.supportingSpecies : []).forEach((item) => {
    const label = String(item?.species || "").trim();
    const lowerLabel = label.toLowerCase();
    if (
      (dominantScientific && lowerLabel.includes(dominantScientific)) ||
      (dominantDisplay && dominantDisplay !== "--" && lowerLabel.includes(dominantDisplay))
    ) {
      extraIntraspeciesSignals.push({
        label,
        ratio: Number(item?.ratio || 0),
      });
      return;
    }
    filteredSupportingSpecies.push(item);
  });
  const mergedIntraspeciesSignals = [
    ...(Array.isArray(meta.intraspeciesSignals) ? meta.intraspeciesSignals : []),
    ...extraIntraspeciesSignals,
  ];
  const intraspeciesMarkup = mergedIntraspeciesSignals.length
    ? `<p>同种内分型/亚群信号：${escapeHtml(mergedIntraspeciesSignals.map((item) => `${item.label}（${Number(item.ratio || 0).toFixed(2)}%）`).join("；"))}</p>`
    : "";
  const supportSpeciesMarkup = filteredSupportingSpecies.length
    ? `<p>背景共存物种：${escapeHtml(filteredSupportingSpecies.map((item) => `${item.species}（${Number(item.ratio || 0).toFixed(2)}%）`).join("；"))}</p>`
    : `<p>当前未见占比达到重点阈值的次高异种背景，后续分析可优先聚焦主导物种。</p>`;
  const mlstHeadline = String(meta?.mlstKnowledge?.headline || "").trim();
  const mlstItems = Array.isArray(meta?.mlstKnowledge?.items) ? meta.mlstKnowledge.items : [];
  const serotypeKnowledgeHeadline = String(meta?.serotypeKnowledge?.headline || "").trim();
  const serotypeKnowledgeItems = Array.isArray(meta?.serotypeKnowledge?.items) ? meta.serotypeKnowledge.items : [];
  const prioritySerotypes = Array.isArray(meta?.prioritySerotypeKnowledge) ? meta.prioritySerotypeKnowledge : [];
  const serotypeFragments = prioritySerotypes.map((item) => {
    const fragments = [];
    const typedLabel = [item.species, item.serotype].filter(Boolean).join(" ");
    if (typedLabel) fragments.push(`关注血清型命中：${typedLabel}`);
    if (item.panel && item.panel !== "-") fragments.push(`监测面板：${item.panel}`);
    if (item.virulence && item.virulence !== "-") fragments.push(`毒力关联：${item.virulence}`);
    if (item.resistance && item.resistance !== "-") fragments.push(`耐药关联：${item.resistance}`);
    if (item.regional && item.regional !== "-") fragments.push(`地域分布：${item.regional}`);
    if (item.interpretation && item.interpretation !== "-") fragments.push(`解释提示：${item.interpretation}`);
    return fragments.join("。");
  }).filter(Boolean);
  const mlstFragments = [];
  if (mlstHeadline) mlstFragments.push(mlstHeadline);
  mlstItems.forEach((item) => {
    const fragments = [];
    if (item?.lineage_text && item.lineage_text !== "-") fragments.push(`克隆复合群/Lineage：${item.lineage_text}`);
    if (Array.isArray(item?.virulence) && item.virulence.length) fragments.push(`毒力关联：${item.virulence.join("；")}`);
    if (Array.isArray(item?.resistance) && item.resistance.length) fragments.push(`耐药关联：${item.resistance.join("；")}`);
    if (Array.isArray(item?.regional) && item.regional.length) fragments.push(`地域分布：${item.regional.join("；")}`);
    if (item?.interpretation) fragments.push(`解释提示：${item.interpretation}`);
    if (fragments.length) mlstFragments.push(fragments.join("。"));
  });
  const serotypeKnowledgeFragments = [];
  if (serotypeKnowledgeHeadline) serotypeKnowledgeFragments.push(serotypeKnowledgeHeadline);
  serotypeKnowledgeItems.forEach((item) => {
    const fragments = [];
    const label = String(item?.serotype || item?.matched_on || "").trim();
    if (label) fragments.push(label);
    if (Array.isArray(item?.regional) && item.regional.length) fragments.push(`地域分布：${item.regional.join("；")}`);
    if (Array.isArray(item?.key_markers) && item.key_markers.length) fragments.push(`关键标记：${item.key_markers.join("；")}`);
    if (item?.interpretation) fragments.push(`解释提示：${item.interpretation}`);
    if (fragments.length) serotypeKnowledgeFragments.push(fragments.join("。"));
  });
  const neisseriaAmrFragments = [];
  if (shouldShowNeisseriaAmr && neisseriaAmr.available) {
    if (neisseriaAmr.headline) neisseriaAmrFragments.push(neisseriaAmr.headline);
    if (neisseriaAmr.interpretationItems.length) {
      neisseriaAmrFragments.push(...neisseriaAmr.interpretationItems.slice(0, 4));
    } else if (neisseriaAmr.highlights.length) {
      neisseriaAmrFragments.push(`重点命中位点：${neisseriaAmr.highlights.slice(0, 3).join("；")}`);
    }
  }
  const tbAmrFragments = [];
  const molecularAttentionCards = [
    {
      tone: inferResearchAttentionTone(meta?.riskLevel),
      title: "重点提示",
      body: meta.conclusion,
    },
  ];
  if (shouldShowTbAmr && tbAmr.available) {
    const gradeLabel = String(tbAmr?.resistanceGrade?.label || "").trim();
    const gradeReason = String(tbAmr?.resistanceGrade?.reason || "").trim();
    if (gradeLabel) {
      tbAmrFragments.push(`结核耐药分级：${gradeLabel}${gradeReason ? `（${gradeReason}）` : ""}`);
      if (gradeLabel !== "敏感") {
        molecularAttentionCards.push({
          tone: String(tbAmr?.resistanceGrade?.tone || "").trim() || "watch",
          title: "需要高度关注",
          body: `结核耐药分级为 ${gradeLabel}${gradeReason ? `，${gradeReason}` : ""}`,
        });
      }
    }
    if (tbAmr.headline && !gradeLabel) tbAmrFragments.push(tbAmr.headline);
    if (tbAmr.focusCalls.length) {
      tbAmrFragments.push(`重点耐药药物：${tbAmr.focusCalls.map((item) => {
        const drug = String(item?.drug || "").trim();
        const verdict = String(item?.verdict || "").trim();
        const mutations = Array.isArray(item?.mutations) ? item.mutations.filter(Boolean).slice(0, 3).join("、") : "";
        return [drug, verdict ? `：${verdict}` : "", mutations ? `（${mutations}）` : ""].join("");
      }).filter(Boolean).join("；")}`);
    } else if (tbAmr.interpretationItems.length) {
      tbAmrFragments.push(...tbAmr.interpretationItems.slice(0, 2));
    }
  }
  const molecularFeatureParagraphs = [`<p>${escapeHtml(meta.conclusion)}</p>`];
  if (meta.eventSummaries.length) {
    molecularFeatureParagraphs.push(`<p>组合事件结果显示：${escapeHtml(meta.eventSummaries.join("；"))}。</p>`);
  }
  if (meta.geneLabels.length || meta.hitGenes.length) {
    molecularFeatureParagraphs.push(`<p>${meta.geneLabels.length ? `重点规则命中包括：${escapeHtml(meta.geneLabels.join("、"))}。` : ""}${meta.hitGenes.length ? `${meta.geneLabels.length ? " " : ""}实际命中基因包括：${escapeHtml(meta.hitGenes.join("、"))}。` : ""}</p>`);
  }
  if (mlstFragments.length) {
    molecularFeatureParagraphs.push(`<p>MLST 分型结果支持：${escapeHtml(mlstFragments.join("；"))}</p>`);
  }
  if (serotypeKnowledgeFragments.length) {
    molecularFeatureParagraphs.push(`<p>分型/家系知识库结果支持：${escapeHtml(serotypeKnowledgeFragments.join("；"))}</p>`);
  } else if (serotypeFragments.length) {
    molecularFeatureParagraphs.push(`<p>血清型/血清群结果支持：${escapeHtml(serotypeFragments.join("；"))}</p>`);
  }
  if (tbAmrFragments.length) {
    molecularFeatureParagraphs.push(`<p>结核耐药突变结果支持：${escapeHtml(tbAmrFragments.join("；"))}</p>`);
  }
  if (neisseriaAmrFragments.length) {
    molecularFeatureParagraphs.push(`<p>脑膜炎奈瑟菌耐药突变证据链提示：${escapeHtml(neisseriaAmrFragments.join("；"))}。上述结果更适合解释为基于位点和单倍型背景的分子证据，建议结合完整序列背景、系统发育关系及药敏结果进行综合判读。</p>`);
  }
  return `
    <article class="scene-report-card research-scene-card">
      <header class="scene-report-head">
        <span class="scene-report-kicker">Research Interpretation</span>
        <h3>科研判读摘要</h3>
      </header>
      <div class="scene-report-sections">
        <section>
          <h4>主导研究对象</h4>
          <p>当前样本的主导物种线索为 <strong>${escapeHtml(meta.speciesName)}</strong>（NCBI TaxID: ${escapeHtml(meta.taxid)}；学名：${escapeHtml(meta.scientificName)}）。知识库风险层级评估为 <strong>${escapeHtml(meta.riskLevel)}</strong>。</p>
          <p>${escapeHtml(meta.significance)}</p>
        </section>
        <section>
          <h4>分子特征提示</h4>
          ${molecularAttentionCards.length ? `
            <div class="research-attention-stack">
              ${molecularAttentionCards.map((card) => `
                <article class="research-attention-card tone-${escapeHtml(card.tone)}">
                  <span class="research-attention-kicker">${escapeHtml(card.title)}</span>
                  <p>${escapeHtml(card.body)}</p>
                </article>
              `).join("")}
            </div>
          ` : ""}
          ${molecularFeatureParagraphs.join("")}
        </section>
        <section>
          <h4>证据边界</h4>
          <p>${escapeHtml(meta.confidenceHint)}</p>
          ${intraspeciesMarkup}
          ${supportSpeciesMarkup}
        </section>
        <section>
          <h4>科研提示</h4>
          <p>${escapeHtml(meta.nextSteps.join(" "))}</p>
        </section>
      </div>
    </article>
  `;
}

function renderSceneSpecificSection(data) {
  const body = document.getElementById("clinical-report-body");
  const section = document.getElementById("section-clinical-report");
  if (!body || !section) return;
  if (currentReportScenario === "clinical") {
    body.innerHTML = buildClinicalScene(data);
    section.classList.remove("hidden");
  } else if (currentReportScenario === "cdc") {
    body.innerHTML = buildCdcScene(data);
    section.classList.remove("hidden");
  } else {
    body.innerHTML = "";
    section.classList.add("hidden");
  }
}

function updateReportScenarioLayout(data) {
  const scenarioSupported = supportsReportScenario(data);
  const isVirusReport = isVirusScenarioReport(data);
  const titleNode = document.querySelector(".report-title-block h1");
  const subtitleNode = document.querySelector(".report-subtitle");
  const navClinical = document.getElementById("nav-group-clinical-wrapper");
  const summaryStrip = document.querySelector(".report-summary-strip");
  const introCard = document.querySelector(".report-intro-card");
  const shell = document.querySelector(".report-shell");
  if (!scenarioSupported) {
    currentReportScenario = "research";
  }
  if (shell) {
    shell.dataset.reportScene = currentReportScenario;
  }
  if (titleNode) {
    titleNode.textContent = currentReportScenario === "clinical"
      ? (isVirusReport ? "病毒检测医院报告" : "单菌测序临床辅助诊断报告")
      : currentReportScenario === "cdc"
        ? (isVirusReport ? "病毒疾控监测报告" : "单菌测序疾控监测报告")
        : "病原微生物分析结果";
  }
  if (subtitleNode) {
    subtitleNode.textContent = currentReportScenario === "clinical"
      ? (isVirusReport ? "面向临床阅读，突出病毒检出、分型证据、临床意义和复核建议。" : "面向临床医生阅读，突出病原体意义、耐药风险、毒力提示和原则性诊疗建议。")
      : currentReportScenario === "cdc"
        ? (isVirusReport ? "面向监测与流调场景，突出病毒分型、传播风险、公共卫生提示和处置建议。" : "面向监测与流调场景，突出病原体识别、分型特征、传播风险和公共卫生提示。")
        : "本报告按分析流程整理样本质控、组装、物种、耐药毒力及分型结果，供结果判读与存档使用。";
  }
  if (navClinical) {
    navClinical.classList.toggle("hidden", currentReportScenario === "research");
    const navLabel = navClinical.querySelector('a[href="#section-clinical-report"]');
    if (navLabel) {
      navLabel.textContent = currentReportScenario === "cdc"
        ? (isVirusReport ? "病毒疾控报告" : "疾控监测报告")
        : (isVirusReport ? "病毒医院报告" : "临床辅助诊断");
    }
  }
  const clinicalHeading = document.querySelector("#section-clinical-report .section-heading h2");
  const clinicalCopy = document.querySelector("#section-clinical-report .section-heading p:last-child");
  if (clinicalHeading) {
    clinicalHeading.textContent = currentReportScenario === "cdc"
      ? (isVirusReport ? "病毒疾控监测与风险评估报告" : "疾控监测与风险评估报告")
      : (isVirusReport ? "病毒检测医院报告" : "临床辅助诊断报告");
  }
  if (clinicalCopy) {
    clinicalCopy.textContent = currentReportScenario === "cdc"
      ? (isVirusReport ? "面向疾控与流行病学专业人员，强调病毒分型、传播风险、监测意义和防控建议。" : "面向疾控与流行病学专业人员，强调病原体公共卫生意义、传播风险和监测防控建议。")
      : (isVirusReport ? "面向临床医生的正式阅读视图，突出病毒检出、分型证据、临床意义和复核建议。" : "面向临床医生的正式阅读视图，突出病原体意义、耐药风险、毒力提示和原则性诊疗建议。");
  }
  const contentSections = Array.from(document.querySelectorAll(".report-content > .report-section"));
  const clinicalSection = document.getElementById("section-clinical-report");
  contentSections.forEach((section) => {
    if (!section?.id) return;
    if (currentReportScenario === "clinical") {
      section.classList.toggle("hidden", section !== clinicalSection);
      return;
    }
    if (currentReportScenario === "cdc") {
      section.classList.toggle("hidden", section !== clinicalSection);
      return;
    }
    section.classList.toggle("hidden", section === clinicalSection);
  });
  if (summaryStrip) {
    summaryStrip.classList.toggle("hidden", currentReportScenario !== "research");
  }
  if (introCard) {
    introCard.classList.toggle("hidden", currentReportScenario === "clinical" || currentReportScenario === "cdc");
  }
  renderSceneSpecificSection(data);
  const navGroups = Array.from(document.querySelectorAll(".report-nav > .report-nav-group"));
  if (currentReportScenario === "clinical" || currentReportScenario === "cdc") {
    navGroups.forEach((group) => {
      group.classList.toggle("hidden", group !== navClinical);
    });
  } else {
    const navPairs = [
      ["#section-overview", 'a[href="#section-overview"]'],
      ["#section-raw-qc", '[data-nav-section="section-raw-qc"]'],
      ["#section-species", '[data-nav-section="section-species"]'],
      ["#section-assembly", '[data-nav-section="section-assembly"]'],
      ["#section-rv", '[data-nav-section="section-rv"]'],
      ["#section-mlst", 'a[href="#section-mlst"]'],
      ["#section-serotype", 'a[href="#section-serotype"]'],
      ["#section-priority-serotype", 'a[href="#section-priority-serotype"]'],
      ["#section-mge", '[data-nav-section="section-mge"]'],
    ];
    navGroups.forEach((group) => {
      if (group === navClinical) {
        group.classList.add("hidden");
      } else {
        group.classList.remove("hidden");
      }
    });
    navPairs.forEach(([sectionSelector, navSelector]) => {
      const sectionNode = document.querySelector(sectionSelector);
      const navNode = document.querySelector(navSelector)?.closest(".report-nav-group");
      if (!sectionNode || !navNode) return;
      navNode.classList.toggle("hidden", sectionNode.classList.contains("hidden"));
    });
  }
}

function isMobileCompactReportKind(reportKind) {
  return window.matchMedia("(max-width: 720px)").matches && MOBILE_COMPACT_REPORT_KINDS.has(String(reportKind || "").trim());
}

function setMobileSectionOpenState(section, open) {
  if (!section) return;
  section.classList.toggle("is-mobile-open", open);
  const toggle = section.querySelector(":scope > .section-heading .report-mobile-toggle");
  if (toggle) {
    toggle.textContent = open ? "收起" : "展开";
    toggle.setAttribute("aria-expanded", open ? "true" : "false");
  }
}

function ensureMobileSectionVisible(target) {
  const section = target?.closest?.(".report-section.is-mobile-collapsible");
  if (!section) return;
  setMobileSectionOpenState(section, true);
}

function applyMobileReportCompaction(reportKind) {
  const shell = document.querySelector(".report-shell");
  if (!shell) return;
  const enabled = isMobileCompactReportKind(reportKind) && currentReportScenario === "research";
  shell.dataset.mobileCompact = enabled ? "virus" : "default";
  const sections = Array.from(document.querySelectorAll(".report-content > .report-section"));
  sections.forEach((section) => {
    section.classList.remove("is-mobile-collapsible", "is-mobile-open");
    const toggle = section.querySelector(":scope > .section-heading .report-mobile-toggle");
    if (toggle) toggle.remove();
    const wrappedBody = section.querySelector(":scope > .report-section-body-wrap");
    if (wrappedBody) {
      while (wrappedBody.firstChild) {
        section.appendChild(wrappedBody.firstChild);
      }
      wrappedBody.remove();
    }
  });
  if (!enabled) return;
  const defaultOpen = new Set(["section-overview"]);
  sections.forEach((section) => {
    if (!section?.id || section.classList.contains("hidden")) return;
    if (section.id === "section-clinical-report") return;
    const heading = section.querySelector(":scope > .section-heading");
    if (!heading) return;
    const bodyWrap = document.createElement("div");
    bodyWrap.className = "report-section-body-wrap";
    let sibling = heading.nextElementSibling;
    while (sibling) {
      const next = sibling.nextElementSibling;
      bodyWrap.appendChild(sibling);
      sibling = next;
    }
    section.appendChild(bodyWrap);
    const toggle = document.createElement("button");
    toggle.type = "button";
    toggle.className = "report-mobile-toggle";
    toggle.addEventListener("click", () => {
      setMobileSectionOpenState(section, !section.classList.contains("is-mobile-open"));
    });
    heading.appendChild(toggle);
    section.classList.add("is-mobile-collapsible");
    setMobileSectionOpenState(section, defaultOpen.has(section.id));
  });
  if (window.location.hash) {
    const target = document.querySelector(window.location.hash);
    if (target) ensureMobileSectionVisible(target);
  }
}

function buildExecutiveSummary(data) {
  const metrics = Array.isArray(data?.overview_metrics) ? data.overview_metrics : [];
  const sections = data?.sections || {};
  const isMetaMethod = String(data?.task?.method || "").trim() === "meta";

  if (isMetaMethod) {
    const assemblyMetric = getMetricByKey(metrics, "meta_assembly");
    const binningMetric = getMetricByKey(metrics, "meta_binning_quality");
    const speciesMetric = getMetricByKey(metrics, "meta_species_mge");
    const sequencingMetric = getMetricByKey(metrics, "meta_sequencing");
    const binningSummary = sections?.binning_results?.quality?.summary || {};
    const hq = Number(binningSummary.hq_bins || 0);
    const mq = Number(binningSummary.mq_bins || 0);
    const lq = Number(binningSummary.lq_bins || 0);
    const dominantSpecies = speciesMetric?.items?.[0]?.display || "--";
    const dominantRatio = Number(sections?.species_identification?.abundance?.ranks?.find((rank) => rank.rank === "种")?.segments?.[0]?.ratio || 0);
    const mgeTotal = speciesMetric?.items?.[1]?.display || "--";
    const sequencingBody = `总测序数据量 ${sequencingMetric?.items?.[0]?.display || "--"}，读取记录 ${sequencingMetric?.items?.[1]?.display || "--"}。`;
    const assemblyBody = `当前组装得到 ${assemblyMetric?.items?.[0]?.display || "--"} 条 Contig，总长度 ${assemblyMetric?.items?.[1]?.display || "--"}。`;
    const binningBody = `高质量 ${hq} 个，中质量 ${mq} 个，低质量 ${lq} 个；平均完整性 ${binningMetric?.items?.[0]?.display || "--"}，平均污染率 ${binningMetric?.items?.[1]?.display || "--"}。`;
    const speciesBody = dominantSpecies && dominantSpecies !== "--"
      ? `当前优势物种为 ${dominantSpecies}，占比约 ${dominantRatio.toFixed(2)}%；移动元件共检出 ${mgeTotal} 条。`
      : `当前未形成稳定的优势物种判读，移动元件共检出 ${mgeTotal} 条。`;
    return [
      { icon: "Q", label: "Sequencing", title: "测序数据概览", body: sequencingBody, meta: "测序数据量 + 读取记录", state: "neutral", target: "#section-overview" },
      { icon: "A", label: "Assembly", title: "组装结果已汇总", body: assemblyBody, meta: "基因组组装", state: "neutral", target: "#section-assembly" },
      { icon: "B", label: "Binning", title: "Binning完整性概览", body: binningBody, meta: "Binning结果", state: hq > 0 ? "success" : mq > 0 ? "attention" : "warning", target: "#section-binning" },
      { icon: "T", label: "Taxonomy", title: dominantSpecies || "优势物种未定", body: speciesBody, meta: "优势物种 + 移动元件", state: dominantSpecies && dominantSpecies !== "--" ? "neutral" : "attention", target: "#section-species" },
    ];
  }

  const qMetric = getMetricByKey(metrics, "q_metrics");
  const q20 = parsePercentDisplay(qMetric?.items?.[0]?.display);
  const q30 = parsePercentDisplay(qMetric?.items?.[1]?.display);
  const qualityState = getQualityCardState(q20, q30);
  const hasQualityMetrics = q20 !== null && q30 !== null;
  const qualityTitle = !hasQualityMetrics
    ? "测序质量未生成"
    : qualityState === "danger"
      ? "测序质量需重点关注"
      : "测序质量总体合格";
  const qualityBody = !hasQualityMetrics
    ? "当前样本未形成可判读的 Q20/Q30 结果，可能因流程中断或质控结果未生成。"
    : qualityState === "danger"
      ? `Q20 ${qMetric?.items?.[0]?.display || "--"}、Q30 ${qMetric?.items?.[1]?.display || "--"}，建议优先复核原始质控与 fastp 结果。`
      : `Q20 ${qMetric?.items?.[0]?.display || "--"}、Q30 ${qMetric?.items?.[1]?.display || "--"}，原始测序质量达到当前报告判读阈值。`;

  const serotypeSection = sections?.serotype || {};
  const isInfluenzaMode = String(serotypeSection?.mode || "").trim() === "influenza_typing";
  const isMonkeypoxMode = String(serotypeSection?.mode || "").trim() === "monkeypox_nextclade";
  if (isInfluenzaMode) {
    const summaryCards = Array.isArray(serotypeSection?.summary_cards) ? serotypeSection.summary_cards : [];
    const influenzaType = String(
      serotypeSection?.influenza_type
      || summaryCards.find((item) => item?.label === "流感类型")?.value
      || "--",
    ).trim() || "--";
    const haSubtype = String(
      serotypeSection?.ha_subtype
      || summaryCards.find((item) => item?.label === "HA 亚型")?.value
      || "--",
    ).trim() || "--";
    const naSubtype = String(
      serotypeSection?.na_subtype
      || summaryCards.find((item) => item?.label === "NA 亚型")?.value
      || "--",
    ).trim() || "--";
    const segmentManifest = serotypeSection?.segment_manifest && typeof serotypeSection.segment_manifest === "object"
      ? serotypeSection.segment_manifest
      : { columns: [], rows: [] };
    const segmentColumns = Array.isArray(segmentManifest?.columns) ? segmentManifest.columns : [];
    const segmentRows = Array.isArray(segmentManifest?.rows) ? segmentManifest.rows : [];
    const segmentNameIndex = segmentColumns.findIndex((value) => {
      const text = String(value || "").trim();
      return text.includes("segment_group") || text.includes("片段") || text.toLowerCase().includes("segment");
    });
    const subtypeIndex = segmentColumns.findIndex((value) => {
      const text = String(value || "").trim();
      return text.includes("subtype") || text.includes("亚型");
    });
    const segmentSummaries = dedupeList(
      segmentRows.map((row) => {
        const segmentName = String(segmentNameIndex >= 0 ? row?.[segmentNameIndex] : "").trim();
        if (!segmentName) return "";
        const segmentSubtype = String(subtypeIndex >= 0 ? row?.[subtypeIndex] : "").trim();
        return segmentSubtype && segmentSubtype !== "-" ? `${segmentName}(${segmentSubtype})` : segmentName;
      }).filter(Boolean),
      16,
    );
    const assemblyTitle = segmentSummaries.length ? `${segmentSummaries.length} 个 segment 已纳入分析` : "未形成稳定的 segment 结果";
    const assemblyBody = segmentSummaries.length
      ? `当前样本的 segment 结果为：${segmentSummaries.join("、")}。`
      : "当前未形成可判读的流感 segment 列表，建议优先检查分型与组装步骤是否完整。";
    const speciesTitle = influenzaType === "Influenza B virus"
      ? "乙流"
      : influenzaType === "Influenza A virus"
        ? "甲流"
        : (influenzaType || "流感类型未定");
    const speciesBody = influenzaType === "Influenza A virus"
      ? `当前判定为甲型流感，HA 分型为 ${haSubtype}，NA 分型为 ${naSubtype}。`
      : influenzaType === "Influenza B virus"
        ? "当前判定为乙型流感。"
        : "当前未形成稳定的甲/乙流判断。";
    return [
      { icon: "Q", label: "Sequencing", title: qualityTitle, body: qualityBody, meta: "原始数据质控 + fastp", state: qualityState, target: "#section-raw-qc" },
      { icon: "A", label: "Assembly", title: assemblyTitle, body: assemblyBody, meta: "组装情况", state: segmentSummaries.length ? "neutral" : "attention", target: "#section-serotype" },
      { icon: "T", label: "Taxonomy", title: speciesTitle, body: speciesBody, meta: "物种预估", state: influenzaType !== "--" ? "neutral" : "attention", target: "#section-serotype" },
      { icon: "!", label: "Priority", title: "建议先看流感分型", body: "当前结果更适合按流感类型、HA/NA 分型、segment 组成和变异注释的顺序阅读。", meta: "优先阅读建议", state: "neutral", target: "#section-serotype" },
    ];
  }
  if (isMonkeypoxMode) {
    const summaryCards = Array.isArray(serotypeSection?.summary_cards) ? serotypeSection.summary_cards : [];
    const clade = String(
      serotypeSection?.predicted_clade
      || summaryCards.find((item) => item?.label === "Clade")?.value
      || "--",
    ).trim() || "--";
    const lineage = String(
      serotypeSection?.predicted_lineage
      || summaryCards.find((item) => item?.label === "Lineage")?.value
      || "--",
    ).trim() || "--";
    const outbreak = String(
      serotypeSection?.predicted_outbreak
      || summaryCards.find((item) => item?.label === "Outbreak")?.value
      || "--",
    ).trim() || "--";
    const qualityMetrics = Array.isArray(serotypeSection?.quality_metrics) ? serotypeSection.quality_metrics : [];
    const qcLabel = String(
      qualityMetrics.find((item) => String(item?.label || "").toLowerCase().includes("qc"))?.value || "--",
    ).trim() || "--";
    const mutationTable = serotypeSection?.mutation_table && typeof serotypeSection.mutation_table === "object"
      ? serotypeSection.mutation_table
      : { rows: [], high_quality_variants: "--", low_quality_variants: "--" };
    const totalVariants = Array.isArray(mutationTable?.rows) ? mutationTable.rows.length : 0;
    const highVariants = String(mutationTable?.high_quality_variants ?? "--");
    const lowVariants = String(mutationTable?.low_quality_variants ?? "--");
    const coverageMetric = getMetricByKey(metrics, "monkeypox_assembly_coverage");
    const speciesMetric = getMetricByKey(metrics, "monkeypox_species_estimation");
    const coveragePieces = Array.isArray(coverageMetric?.items)
      ? coverageMetric.items.map((item) => `${item?.label || "--"} ${item?.display || "--"}`).filter(Boolean)
      : [];
    const assemblyTitle = clade !== "--"
      ? `猴痘谱系定位为 ${clade}${lineage !== "--" ? ` / ${lineage}` : ""}`
      : "猴痘谱系结果待确认";
    const assemblyBody = coveragePieces.length
      ? `当前覆盖度表现为 ${coveragePieces.join("，")}。`
      : "当前未形成可判读的猴痘覆盖度摘要。";
    const speciesTitle = clade !== "--"
      ? `${clade}${lineage !== "--" ? ` / ${lineage}` : ""}`
      : (speciesMetric?.items?.[0]?.display || "猴痘分型未定");
    const speciesBody = outbreak !== "--"
      ? `Nextclade 结果提示该样本归属于 ${speciesTitle}，并关联 outbreak ${outbreak}。`
      : `Nextclade 结果提示该样本归属于 ${speciesTitle}。`;
    const focusBody = totalVariants
      ? `当前共整理 ${totalVariants} 个突变位点，其中高质量 ${highVariants} 个、低质量 ${lowVariants} 个，建议结合变异表和 IGV 查看关键位点。`
      : "建议优先查看猴痘分型总表、质量指标和 IGV 比对结果。";
    return [
      { icon: "Q", label: "Sequencing", title: qualityTitle, body: qualityBody, meta: "原始数据质控 + fastp", state: qualityState, target: "#section-raw-qc" },
      { icon: "A", label: "Assembly", title: assemblyTitle, body: assemblyBody, meta: "覆盖度 + 组装情况", state: coveragePieces.length ? "neutral" : "attention", target: "#section-assembly" },
      { icon: "T", label: "Typing", title: speciesTitle, body: speciesBody, meta: qcLabel !== "--" ? `Clade / Lineage + QC ${qcLabel}` : "Clade / Lineage", state: clade !== "--" ? "neutral" : "attention", target: "#section-serotype" },
      { icon: "!", label: "Priority", title: "建议先看猴痘分型", body: focusBody, meta: "变异位点 + IGV", state: totalVariants ? "neutral" : "attention", target: "#section-serotype" },
    ];
  }
  if (["rsv_nextclade", "hmpv_nextclade", "denv_nextclade", "zikav_nextclade", "chikv_nextclade", "ebola_nextclade", "hpiv_typing", "hiv_resistance", "hadv_typing", "norovirus_typing", "enterovirus_typing", "hepatovirus_typing", "bandavirus_typing", "orthohantavirus_typing", "orthoebolavirus_typing", "astroviridae_typing", "rhinovirus_typing", "seasonal_hcov_typing", "rotavirus_typing"].includes(String(serotypeSection?.mode || "").trim())) {
    const meta = extractRespiratoryNextcladeResearchInterpretation(data);
    const coverageMetric = getMetricByKey(
      metrics,
      String(serotypeSection?.mode || "").trim() === "hmpv_nextclade"
        ? "hmpv_assembly_coverage"
      : ((String(serotypeSection?.mode || "").trim() === "denv_nextclade" || String(serotypeSection?.mode || "").trim() === "zikav_nextclade" || String(serotypeSection?.mode || "").trim() === "chikv_nextclade")
          ? "denv_assembly_coverage"
          : (String(serotypeSection?.mode || "").trim() === "hpiv_typing"
            ? "hpiv_assembly_coverage"
            : (String(serotypeSection?.mode || "").trim() === "hiv_resistance"
              ? "hiv_assembly_coverage"
            : (String(serotypeSection?.mode || "").trim() === "hadv_typing"
              ? "hadv_assembly_coverage"
                : (String(serotypeSection?.mode || "").trim() === "norovirus_typing"
                ? "norovirus_assembly_coverage"
                : (String(serotypeSection?.mode || "").trim() === "enterovirus_typing"
                  ? "rhinovirus_assembly_coverage"
                  : (String(serotypeSection?.mode || "").trim() === "hepatovirus_typing"
                    ? "hepatovirus_assembly_coverage"
                  : (String(serotypeSection?.mode || "").trim() === "bandavirus_typing"
                    ? "bandavirus_assembly_coverage"
                  : (String(serotypeSection?.mode || "").trim() === "orthohantavirus_typing"
                    ? "orthohantavirus_assembly_coverage"
                  : (String(serotypeSection?.mode || "").trim() === "astroviridae_typing"
                    ? "astroviridae_assembly_coverage"
                  : (String(serotypeSection?.mode || "").trim() === "rhinovirus_typing"
                  ? "rhinovirus_assembly_coverage"
                  : (String(serotypeSection?.mode || "").trim() === "seasonal_hcov_typing"
                    ? "seasonal_hcov_assembly_coverage"
                    : (String(serotypeSection?.mode || "").trim() === "rotavirus_typing" ? "rotavirus_assembly_coverage" : "rsv_assembly_coverage"))))))))))))),
    );
    const coveragePieces = Array.isArray(coverageMetric?.items)
      ? coverageMetric.items.map((item) => `${item?.label || "--"} ${item?.display || "--"}`).filter(Boolean)
      : [];
    const assemblyTitle = (meta.isHpiv || meta.isHadv || meta.isNorovirus || meta.isEnterovirus || meta.isHepatovirus || meta.isBandavirus || meta.isOrthohantavirus || meta.isAstroviridae || meta.isRhinovirus || meta.isSeasonalHcov || meta.isRotavirus)
      ? (meta.clade !== "--" ? `${meta.virusShort} 分型定位为 ${meta.clade}` : `${meta.virusShort} 分型结果待确认`)
      : (meta.clade !== "--"
        ? `${meta.virusShort} 谱系定位为 ${meta.clade}${meta.lineage !== "--" ? ` / ${meta.lineage}` : ""}`
        : `${meta.virusShort} 谱系结果待确认`);
    const assemblyBody = coveragePieces.length
      ? `当前覆盖度表现为 ${coveragePieces.join("，")}。`
      : "当前未形成可判读的覆盖度摘要。";
    const speciesTitle = meta.clade !== "--"
      ? `${meta.clade}${!meta.isHpiv && !meta.isRhinovirus && !meta.isSeasonalHcov && meta.lineage !== "--" ? ` / ${meta.lineage}` : ""}`
      : `${meta.virusShort} 分型未定`;
    const speciesBody = (meta.isHpiv || meta.isHiv || meta.isHadv || meta.isHepatovirus || meta.isBandavirus || meta.isOrthohantavirus || meta.isRhinovirus || meta.isSeasonalHcov)
      ? `最优参考比对结果提示该样本更接近 ${speciesTitle}。`
      : `Nextclade 结果提示该样本归属于 ${speciesTitle}。`;
    const focusBody = meta.mutationCount
      ? `当前共整理 ${meta.mutationCount} 个突变位点，其中高质量 ${meta.highMutationCount} 个、低质量 ${meta.lowMutationCount} 个，建议结合变异表和 IGV 查看关键位点。`
      : `建议优先查看 ${meta.virusShort} 分型总表、质量指标和 IGV 比对结果。`;
    return [
      { icon: "Q", label: "Sequencing", title: qualityTitle, body: qualityBody, meta: "原始数据质控 + fastp", state: qualityState, target: "#section-raw-qc" },
      { icon: "A", label: "Assembly", title: assemblyTitle, body: assemblyBody, meta: "覆盖度 + 组装情况", state: coveragePieces.length ? "neutral" : "attention", target: "#section-assembly" },
      { icon: "T", label: "Typing", title: speciesTitle, body: speciesBody, meta: (meta.isHpiv || meta.isHiv || meta.isHadv || meta.isHepatovirus || meta.isRhinovirus || meta.isSeasonalHcov) ? "参考分型 + 覆盖度筛选" : (meta.qualityLabel !== "--" ? `Clade / Lineage + QC ${meta.qualityLabel}` : "Clade / Lineage"), state: meta.clade !== "--" ? "neutral" : "attention", target: "#section-serotype" },
      { icon: "!", label: "Priority", title: `建议先看${meta.virusShort}分型`, body: focusBody, meta: "变异位点 + IGV", state: meta.mutationCount ? "neutral" : "attention", target: "#section-serotype" },
    ];
  }

  const assemblyMetric = getMetricByKey(metrics, "assembly_profile");
  const contigCount = parseOptionalNumber(assemblyMetric?.contig_count);
  const plasmidCount = parseOptionalNumber(assemblyMetric?.plasmid_count);
  const totalLength = parseOptionalNumber(assemblyMetric?.total_length);
  const checkmMetric = getMetricByKey(metrics, "checkm_metrics");
  const completeness = parsePercentDisplay(checkmMetric?.items?.[0]?.display);
  const contamination = parsePercentDisplay(checkmMetric?.items?.[1]?.display);
  const assemblyAssessment = assessAssemblyQuality({ totalLength, contigCount, completeness, contamination });
  const assemblyState = assemblyAssessment.state;
  const assemblyTitle = assemblyAssessment.title;
  const assemblyBody = !assemblyAssessment.hasMetrics
    ? assemblyAssessment.body
    : `当前共检出 ${contigCount !== null ? contigCount : "--"} 条 Contig、${plasmidCount !== null ? plasmidCount : "--"} 条质粒；总长度 ${totalLength !== null ? formatBases(totalLength) : "--"}；完整性 ${checkmMetric?.items?.[0]?.display || "--"}，污染率 ${checkmMetric?.items?.[1]?.display || "--"}。${assemblyAssessment.note ? ` 关注点：${assemblyAssessment.note}。` : ""}`;

  const speciesMetric = getMetricByKey(metrics, "species_estimation");
  const speciesName = speciesMetric?.items?.[0]?.display || "--";
  const mlstSpecies = speciesMetric?.items?.[1]?.display || "--";
  const speciesTitle = speciesName && speciesName !== "--" ? speciesName : "未形成稳定物种预估";
  const speciesBody = mlstSpecies && mlstSpecies !== "--"
    ? `综合 CheckM 与 MLST 结果，当前样本更接近 ${mlstSpecies} 对应的分型背景。`
    : "当前报告未提供稳定的 MLST 物种辅助信息，建议以分类与分型结果综合判断。";

  const resistanceRows = sections?.resistance_virulence?.resistance_elements?.rows?.length || 0;
  const virulenceRows = sections?.resistance_virulence?.virulence_elements?.rows?.length || 0;
  const priorityRows = countPrioritySerotypeHits(sections?.priority_serotype || {});
  let focusState = "neutral";
  let focusTitle = "建议按模块顺序判读";
  let focusBody = "当前未检出需优先跳读的高风险模块，可依次查看组装、分型与耐药毒力结果。";
  if (priorityRows > 0) {
    focusState = "danger";
    focusTitle = "建议优先查看关注毒力血清型";
    focusBody = `当前已检出 ${priorityRows} 条重点关注结果，建议先结合血清型与毒力信息进行判读。`;
  } else if (resistanceRows > 0 || virulenceRows > 0) {
    focusState = resistanceRows > 0 && virulenceRows > 0 ? "warning" : "neutral";
    focusTitle = "建议重点查看耐药与毒力分析";
    focusBody = `当前耐药元件 ${resistanceRows} 条、毒力元件 ${virulenceRows} 条，建议结合关系图与汇总表共同判读。`;
  }

  return [
    { icon: "Q", label: "Sequencing", title: qualityTitle, body: qualityBody, meta: "原始数据质控 + fastp", state: qualityState, target: "#section-raw-qc" },
    { icon: "A", label: "Assembly", title: assemblyTitle, body: assemblyBody, meta: "组装情况 + CheckM", state: assemblyState, target: "#section-assembly" },
    { icon: "T", label: "Taxonomy", title: speciesTitle, body: speciesBody, meta: "物种预估", state: "neutral", target: "#section-species" },
    { icon: "!", label: "Priority", title: focusTitle, body: focusBody, meta: "优先阅读建议", state: focusState, target: priorityRows > 0 ? "#section-priority-serotype" : "#section-rv" },
  ];
}

function renderExecutiveSummary(data) {
  const container = document.getElementById("report-summary-grid");
  const scenePanel = document.getElementById("report-scene-panel");
  const scenarioSupported = supportsReportScenario(data);
  document.querySelectorAll("[data-report-scene]").forEach((button) => {
    const active = button.dataset.reportScene === currentReportScenario;
    button.classList.toggle("active", active);
    button.disabled = !scenarioSupported && button.dataset.reportScene !== "research";
  });
  if (!container || !scenePanel) return;
  const cards = buildExecutiveSummary(data);
  container.innerHTML = cards.map((card) => buildSummaryCard(card)).join("");
  container.classList.toggle("hidden", currentReportScenario !== "research");
  if (currentReportScenario === "research") {
    scenePanel.classList.remove("hidden");
    scenePanel.innerHTML = buildResearchScene(data);
  } else {
    scenePanel.classList.add("hidden");
    scenePanel.innerHTML = "";
  }
  bindSummaryCardJumps(container);
}

function bindReportSceneSwitcher(data) {
  const buttons = Array.from(document.querySelectorAll("[data-report-scene]"));
  if (!buttons.length) return;
  const scenarioSupported = supportsReportScenario(data);
  if (!scenarioSupported) {
    currentReportScenario = "research";
  }
  buttons.forEach((button) => {
    button.onclick = () => {
      const scene = String(button.dataset.reportScene || "research");
      if (!scenarioSupported && scene !== "research") return;
      currentReportScenario = scene;
      renderExecutiveSummary(data);
      updateReportScenarioLayout(data);
      applyMobileReportCompaction(document.querySelector(".report-shell")?.dataset.reportKind || "default");
    };
  });
  renderExecutiveSummary(data);
  updateReportScenarioLayout(data);
  applyMobileReportCompaction(document.querySelector(".report-shell")?.dataset.reportKind || "default");
}

function buildMetricCards(metrics) {
  const totalMetric = metrics.find((metric) => metric.key === "total_bases");
  const qMetric = metrics.find((metric) => metric.key === "q_metrics");
  const mergedMetrics = metrics.filter((metric) => metric.key !== "q_metrics").map((metric) => {
    if (metric.key === "total_bases" && totalMetric && qMetric) {
      return { ...metric, type: "quality_overview", items: qMetric.items || [] };
    }
    return metric;
  });

  return mergedMetrics.map((metric) => {
    if (metric.type === "quality_overview") {
      const q20 = parsePercentDisplay(metric.items?.[0]?.display);
      const q30 = parsePercentDisplay(metric.items?.[1]?.display);
      const cardState = getQualityCardState(q20, q30);
      return `
        <article class="metric-card quality-overview-card metric-state-${cardState}">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <strong>${escapeHtml(metric.display ?? "--")}</strong>
          <div class="paired-metric-grid compact-paired-grid">
            ${(metric.items || []).map((item) => {
              const value = parsePercentDisplay(item.display);
              const itemState = getQMetricState(item.label, value);
              return `
              <div class="paired-metric-item inline-metric-item metric-state-${itemState}">
                <span>${escapeHtml(item.label)}</span>
                <strong>${escapeHtml(item.display ?? "--")}</strong>
                <div class="metric-progress-track" aria-hidden="true">
                  <div class="metric-progress-fill metric-progress-${itemState}" style="width:${clampPercent(value)}%"></div>
                </div>
              </div>
            `;}).join("")}
          </div>
          <div class="metric-status-note">${cardState === "danger" ? "测序质量差" : "测序质量合格"}</div>
        </article>
      `;
    }
    if (metric.type === "paired") {
      const items = (metric.items || []).map((item) => {
        let state = "neutral";
        const value = parsePercentDisplay(item.display);
        if (metric.key === "checkm_metrics") {
          state = item.label.includes("完整") ? getCompletenessState(value) : getContaminationState(value);
        }
        return { ...item, state, numericValue: value };
      });
      return `
        <article class="metric-card paired-metric-card">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <div class="paired-metric-grid">
            ${items.map((item) => `
              <div class="paired-metric-item metric-state-${item.state}">
                <span>${escapeHtml(item.label)}</span>
                <strong>${escapeHtml(item.display ?? "--")}</strong>
                ${metric.key === "checkm_metrics" ? `
                  <div class="metric-progress-track" aria-hidden="true">
                    <div class="metric-progress-fill metric-progress-${item.state}" style="width:${clampPercent(item.numericValue)}%"></div>
                  </div>
                ` : ""}
              </div>
            `).join("")}
          </div>
        </article>
      `;
    }
    if (metric.type === "influenza_segments") {
      const segments = Array.isArray(metric.segments) ? metric.segments.filter((item) => String(item || "").trim()) : [];
      const segmentCount = Number(metric.segment_count || segments.length || 0);
      return `
        <article class="metric-card paired-metric-card">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <div class="paired-metric-grid">
            <div class="paired-metric-item metric-state-neutral">
              <span>Segment 数</span>
              <strong>${escapeHtml(String(segmentCount || "--"))}</strong>
            </div>
            <div class="paired-metric-item metric-state-neutral">
              <span>片段信息</span>
              <strong>${escapeHtml(segments.length ? segments.join("、") : "--")}</strong>
            </div>
          </div>
        </article>
      `;
    }
    if (metric.type === "influenza_species_estimation") {
      const influenzaType = String(metric.influenza_type || "--").trim() || "--";
      const haSubtype = String(metric.ha_subtype || "--").trim() || "--";
      const naSubtype = String(metric.na_subtype || "--").trim() || "--";
      const typeDisplay = influenzaType === "Influenza A virus"
        ? "甲流"
        : influenzaType === "Influenza B virus"
          ? "乙流"
          : influenzaType;
      const subtypeDisplay = influenzaType === "Influenza A virus"
        ? `HA=${haSubtype} / NA=${naSubtype}`
        : influenzaType === "Influenza B virus"
          ? "乙流"
          : "--";
      return `
        <article class="metric-card paired-metric-card">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <div class="paired-metric-grid">
            <div class="paired-metric-item metric-state-neutral">
              <span>流感类型</span>
              <strong>${escapeHtml(typeDisplay)}</strong>
            </div>
            <div class="paired-metric-item metric-state-neutral">
              <span>分型结果</span>
              <strong>${escapeHtml(subtypeDisplay)}</strong>
            </div>
          </div>
        </article>
      `;
    }
    if (metric.type === "assembly_profile") {
      const contig = parseOptionalNumber(metric.contig_count);
      const plasmid = parseOptionalNumber(metric.plasmid_count);
      const totalCount = parseOptionalNumber(metric.total_count);
      const totalLength = parseOptionalNumber(metric.total_length);
      const hasAssemblyMetrics = contig !== null || plasmid !== null || totalCount !== null;
      const total = Math.max(totalCount ?? 0, (contig ?? 0) + (plasmid ?? 0), 1);
      const checkmMetric = metrics.find((item) => item.key === "checkm_metrics");
      const completeness = parsePercentDisplay(checkmMetric?.items?.[0]?.display);
      const contamination = parsePercentDisplay(checkmMetric?.items?.[1]?.display);
      const assemblyAssessment = assessAssemblyQuality({ totalLength, contigCount: contig, completeness, contamination });
      const contigState = assemblyAssessment.state;
      const contigNote = assemblyAssessment.note;
      const contigRatio = Math.max(0, Math.min(1, (contig ?? 0) / total));
      const plasmidRatio = Math.max(0, Math.min(1, (plasmid ?? 0) / total));
      const contigDash = `${(contigRatio * 251.2).toFixed(2)} 251.2`;
      const plasmidDash = `${(plasmidRatio * 251.2).toFixed(2)} 251.2`;
      return `
        <article class="metric-card assembly-metric-card metric-state-${contigState}">
          <span class="metric-label">${escapeHtml(metric.label)}</span>
          <div class="assembly-metric-layout">
            <div class="donut-chart" aria-label="组装情况环形图">
              <svg viewBox="0 0 120 120" role="img">
                <circle class="donut-track" cx="60" cy="60" r="40"></circle>
                <circle class="donut-segment donut-segment-contig" cx="60" cy="60" r="40" stroke-dasharray="${contigDash}"></circle>
                <circle class="donut-segment donut-segment-plasmid" cx="60" cy="60" r="40" stroke-dasharray="${plasmidDash}" stroke-dashoffset="-${(contigRatio * 251.2).toFixed(2)}"></circle>
              </svg>
              <div class="donut-total">${escapeHtml(String(totalCount ?? "--"))}</div>
            </div>
            <div class="assembly-metric-values">
              <div class="paired-metric-item">
                <span>总长度</span>
                <strong>${escapeHtml(formatBases(metric.total_length))}</strong>
              </div>
              <div class="paired-metric-item">
                <span>Contig</span>
                <strong>${escapeHtml(String(contig ?? "--"))}</strong>
              </div>
              <div class="paired-metric-item">
                <span>质粒</span>
                <strong>${escapeHtml(String(plasmid ?? "--"))}</strong>
              </div>
            </div>
          </div>
          <div class="metric-status-note">${contigNote}</div>
        </article>
      `;
    }
    return `
      <article class="metric-card">
        <span class="metric-label">${escapeHtml(metric.label)}</span>
        <strong>${escapeHtml(metric.display ?? "--")}</strong>
      </article>
    `;
  }).join("");
}


function compareTableValues(left, right, direction) {
  const leftNumber = Number(left);
  const rightNumber = Number(right);
  const bothNumeric = Number.isFinite(leftNumber) && Number.isFinite(rightNumber);
  if (bothNumeric) {
    return direction === "asc" ? leftNumber - rightNumber : rightNumber - leftNumber;
  }
  const leftText = String(left ?? "");
  const rightText = String(right ?? "");
  return direction === "asc"
    ? leftText.localeCompare(rightText, "zh-CN", { numeric: true })
    : rightText.localeCompare(leftText, "zh-CN", { numeric: true });
}

function extractNumericValue(value) {
  if (typeof value === "number" && Number.isFinite(value)) return value;
  const text = String(value ?? "").trim();
  if (!text) return null;
  const normalized = text.replaceAll(",", "");
  const match = normalized.match(/-?\d+(?:\.\d+)?/);
  if (!match) return null;
  const numeric = Number(match[0]);
  return Number.isFinite(numeric) ? numeric : null;
}

function inferColumnKinds(columns, rows) {
  return columns.map((column, index) => {
    let numericCount = 0;
    let valueCount = 0;
    rows.forEach((row) => {
      const value = Array.isArray(row) ? row[index] : row[column];
      const text = String(value ?? "").trim();
      if (!text || text === "-") return;
      valueCount += 1;
      if (extractNumericValue(text) !== null) numericCount += 1;
    });
    return valueCount > 0 && numericCount / valueCount >= 0.8 ? "number" : "text";
  });
}

function matchNumericFilter(value, rawFilter) {
  const numeric = extractNumericValue(value);
  if (numeric === null) return false;
  const filter = String(rawFilter ?? "").trim().replace(/\s+/g, "");
  if (!filter) return true;

  const rangeMatch = filter.match(/^(-?\d+(?:\.\d+)?)\s*(?:-|~)\s*(-?\d+(?:\.\d+)?)$/);
  if (rangeMatch) {
    const min = Number(rangeMatch[1]);
    const max = Number(rangeMatch[2]);
    if (!Number.isFinite(min) || !Number.isFinite(max)) return false;
    return numeric >= Math.min(min, max) && numeric <= Math.max(min, max);
  }

  const compareMatch = filter.match(/^(<=|>=|=|<|>)(-?\d+(?:\.\d+)?)$/);
  if (compareMatch) {
    const operator = compareMatch[1];
    const target = Number(compareMatch[2]);
    if (!Number.isFinite(target)) return false;
    if (operator === "<") return numeric < target;
    if (operator === "<=") return numeric <= target;
    if (operator === ">") return numeric > target;
    if (operator === ">=") return numeric >= target;
    return numeric === target;
  }

  const exact = Number(filter);
  if (Number.isFinite(exact)) {
    return numeric === exact;
  }
  return String(value ?? "").toLowerCase().includes(String(rawFilter ?? "").trim().toLowerCase());
}

function formatMutationFrequencyDisplay(column, value) {
  const columnText = String(column || "").trim();
  const raw = String(value ?? "").trim();
  if (!columnText.includes("突变频率") || !raw) return null;
  const pairMatch = raw.match(/^(.+?)\s*\/\s*(\d*\.?\d+)$/);
  if (pairMatch) {
    const depth = String(pairMatch[1] || "").trim();
    const frequency = Number(pairMatch[2]);
    if (Number.isFinite(frequency)) {
      return `${depth} / ${(frequency * 100).toFixed(2).replace(/\.00$/, "")}%`;
    }
  }
  const numeric = Number(raw);
  if (Number.isFinite(numeric) && numeric >= 0 && numeric <= 1) {
    return `${(numeric * 100).toFixed(2).replace(/\.00$/, "")}%`;
  }
  return null;
}

function formatTwoDecimalDisplay(value) {
  if (typeof value === "number" && Number.isFinite(value)) {
    return Number.isInteger(value) ? String(value) : value.toFixed(2);
  }
  const raw = String(value ?? "").trim();
  if (!raw) return raw;
  const percentMatch = raw.match(/^(-?\d+\.\d+)%$/);
  if (percentMatch) {
    return `${Number(percentMatch[1]).toFixed(2)}%`;
  }
  const decimalMatch = raw.match(/^-?\d+\.\d+$/);
  if (decimalMatch) {
    return Number(raw).toFixed(2);
  }
  return raw;
}

function formatReportTableDisplayValue(column, value) {
  return formatMutationFrequencyDisplay(column, value) || formatTwoDecimalDisplay(value) || String(value ?? "-");
}

function renderTableCellContent(value, column = "") {
  const text = formatReportTableDisplayValue(column, value);
  if (text.length <= 30) {
    return `<span>${escapeHtml(text)}</span>`;
  }
  return `<span title="${escapeHtml(text)}">${escapeHtml(`${text.slice(0, 30)}...`)}</span>`;
}

const tableExportStore = new WeakMap();
let currentReportData = null;

function slugifyFilename(value) {
  const normalized = String(value ?? "").trim().replace(/[^\w\u4e00-\u9fff-]+/g, "_");
  return normalized.replace(/^_+|_+$/g, "") || "export";
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function renderTableExportToolbar() {
  return `
    <div class="table-export-toolbar">
      <div class="table-export-menu">
        <button class="table-export-button table-export-trigger" type="button">导出表格</button>
        <div class="table-export-dropdown">
          <button class="table-export-option" type="button" data-table-export-format="xlsx">XLSX</button>
          <button class="table-export-option" type="button" data-table-export-format="tsv">TSV</button>
          <button class="table-export-option" type="button" data-table-export-format="csv">CSV</button>
        </div>
      </div>
    </div>
  `;
}

function renderChartExportToolbar() {
  return `
    <div class="chart-export-toolbar">
      <div class="chart-export-menu">
        <button class="chart-export-button chart-export-trigger" type="button">下载图表</button>
        <div class="chart-export-dropdown">
          <button class="chart-export-option" type="button" data-chart-export-format="svg">SVG</button>
          <button class="chart-export-option" type="button" data-chart-export-format="png">PNG</button>
          <button class="chart-export-option" type="button" data-chart-export-format="pdf">PDF</button>
        </div>
      </div>
    </div>
  `;
}

const SARS_COV_2_GENE_RANGES = [
  { gene: "ORF1a", start: 266, end: 13468 },
  { gene: "ORF1b", start: 13468, end: 21555 },
  { gene: "S", start: 21563, end: 25384 },
  { gene: "ORF3a", start: 25393, end: 26220 },
  { gene: "E", start: 26245, end: 26472 },
  { gene: "M", start: 26523, end: 27191 },
  { gene: "ORF6", start: 27202, end: 27387 },
  { gene: "ORF7a", start: 27394, end: 27759 },
  { gene: "ORF7b", start: 27756, end: 27887 },
  { gene: "ORF8", start: 27894, end: 28259 },
  { gene: "N", start: 28274, end: 29533 },
  { gene: "ORF9b", start: 28284, end: 28577 },
];

const SARS_COV_2_GENE_COLORS = {
  "5'UTR": "#b5bcc7",
  "3'UTR": "#b5bcc7",
  ORF1ab: "#7f93ab",
  ORF1a: "#7f93ab",
  ORF1b: "#8f9d83",
  S: "#8c7f9f",
  ORF3a: "#9d8577",
  E: "#7f9f9a",
  M: "#8d8a74",
  ORF6: "#9a8695",
  ORF7a: "#7e99a8",
  ORF7b: "#9a927b",
  ORF8: "#8f8476",
  N: "#7f8f79",
  ORF9b: "#9b8a85",
  ORF10: "#8e8796",
};

const EBOLA_GENE_COLORS = {
  "5'UTR": "#9aa3ad",
  "3'UTR": "#9aa3ad",
  NP: "#0072b2",
  VP35: "#d55e00",
  VP40: "#009e73",
  GP: "#cc79a7",
  VP30: "#b07d2b",
  VP24: "#e69f00",
  L: "#1b9e9e",
};

const ORTHOHANTAVIRUS_GENE_COLORS = {
  N: "#0072b2",
  GPC: "#009e73",
  Gn: "#56b4e9",
  Gc: "#cc79a7",
  L: "#d55e00",
};

const GENE_TEXT_COLORS = {
  "5'UTR": "#263241",
  "3'UTR": "#263241",
  NP: "#151a20",
  VP35: "#151a20",
  VP40: "#151a20",
  L: "#151a20",
  VP24: "#263241",
  VP30: "#263241",
  GP: "#263241",
  N: "#151a20",
  GPC: "#151a20",
  Gn: "#151a20",
  Gc: "#151a20",
};

function getSarsCov2GeneColor(label, featureType = "") {
  const key = String(label || "").trim();
  if (key && EBOLA_GENE_COLORS[key]) return EBOLA_GENE_COLORS[key];
  if (key && ORTHOHANTAVIRUS_GENE_COLORS[key]) return ORTHOHANTAVIRUS_GENE_COLORS[key];
  if (key && SARS_COV_2_GENE_COLORS[key]) return SARS_COV_2_GENE_COLORS[key];
  if (String(featureType || "").trim() === "five_prime_UTR" || String(featureType || "").trim() === "three_prime_UTR") {
    return SARS_COV_2_GENE_COLORS["5'UTR"];
  }
  return "#8a96a6";
}

function getGenomeFeatureTextColor(label) {
  const key = String(label || "").trim();
  return GENE_TEXT_COLORS[key] || "#151a20";
}

function splitCsvField(value) {
  return String(value || "")
    .split(",")
    .map((item) => item.trim())
    .filter(Boolean);
}

function resolveSarsCov2GenesForRange(start, end = start) {
  const genes = SARS_COV_2_GENE_RANGES
    .filter((item) => start <= item.end && end >= item.start)
    .map((item) => item.gene);
  return genes.length ? genes : ["Intergenic"];
}

function translateNextcladeMutationType(type) {
  const normalized = String(type || "").trim().toLowerCase();
  if (normalized === "substitution") return "替换";
  if (normalized === "deletion") return "缺失";
  if (normalized === "insertion") return "插入";
  if (normalized === "frame shift" || normalized === "frameshift") return "移码";
  return type || "--";
}

function hideAllNextcladeTooltips(root = document) {
  root.querySelectorAll("[data-nextclade-plot-tooltip]").forEach((item) => {
    if (item instanceof HTMLElement) item.hidden = true;
  });
}

function parseNextcladeNucleotideMutations(row) {
  const substitutions = splitCsvField(row?.substitutions).map((token) => {
    const match = token.match(/^[A-Z](\d+)[A-Z]$/i);
    const position = match ? Number(match[1]) : null;
    const genes = position ? resolveSarsCov2GenesForRange(position) : ["Intergenic"];
    return {
      type: "Substitution",
      mutation: token,
      position_label: position != null ? String(position) : "--",
      start: position,
      end: position,
      genes,
    };
  });
  const deletions = splitCsvField(row?.deletions).map((token) => {
    const match = token.match(/^(\d+)-(\d+)$/);
    const start = match ? Number(match[1]) : null;
    const end = match ? Number(match[2]) : null;
    const genes = start != null && end != null ? resolveSarsCov2GenesForRange(start, end) : ["Intergenic"];
    return {
      type: "Deletion",
      mutation: token,
      position_label: start != null && end != null ? `${start}-${end}` : token,
      start,
      end,
      genes,
    };
  });
  const insertions = splitCsvField(row?.insertions).map((token) => {
    const [positionToken = "", inserted = ""] = token.split(":");
    const position = Number(positionToken);
    const genes = Number.isFinite(position) ? resolveSarsCov2GenesForRange(position) : ["Intergenic"];
    return {
      type: "Insertion",
      mutation: token,
      position_label: Number.isFinite(position) ? String(position) : "--",
      start: Number.isFinite(position) ? position : null,
      end: Number.isFinite(position) ? position : null,
      genes,
      inserted,
    };
  });
  return [...substitutions, ...deletions, ...insertions];
}

function parseNextcladeAaMutations(row) {
  const parseTokens = (tokens, type) => splitCsvField(tokens).map((token) => {
    const [gene = "Unknown", change = token] = token.split(":");
    const positionMatch = change.match(/(\d+)/);
    return {
      type,
      mutation: change,
      gene,
      position: positionMatch ? Number(positionMatch[1]) : null,
      label: `${gene}:${change}`,
    };
  });
  return [
    ...parseTokens(row?.aaSubstitutions, "Substitution"),
    ...parseTokens(row?.aaDeletions, "Deletion"),
    ...parseTokens(row?.aaInsertions, "Insertion"),
  ];
}

function parseNextcladeFrameShifts(row) {
  return splitCsvField(row?.frameShifts).map((token) => {
    const match = token.match(/(\d+)(?:-(\d+))?/);
    const start = match ? Number(match[1]) : null;
    const end = match && match[2] ? Number(match[2]) : start;
    return {
      type: "Frame shift",
      mutation: token,
      position_label: start != null && end != null && start !== end ? `${start}-${end}` : (start != null ? String(start) : token),
      start,
      end,
      genes: start != null ? resolveSarsCov2GenesForRange(start, end || start) : ["Intergenic"],
    };
  });
}

function initializeNextcladeGeneSummary() {
  const NEXTCLADE_KNOWLEDGE_CATEGORY_COLUMNS = [
    { key: "transmission", label: "传播/感染力" },
    { key: "antibody", label: "抗体逃逸" },
    { key: "drug", label: "药物耐药" },
    { key: "tcell", label: "T 细胞表位" },
    { key: "pathogenicity", label: "致病性/分子机制" },
    { key: "structure", label: "结构稳定性" },
  ];
  const featureNodes = Array.from(document.querySelectorAll("[data-nextclade-feature-summary]"));
  featureNodes.forEach((node) => {
    const payload = node.getAttribute("data-nextclade-feature-summary");
    if (!payload) return;
    let parsed = null;
    try {
      parsed = JSON.parse(decodeURIComponent(payload));
    } catch (_error) {
      return;
    }
    const select = node.querySelector("[data-nextclade-feature-select]");
    const geneTable = node.querySelector("[data-nextclade-feature-gene-table]");
    const overview = node.querySelector("[data-nextclade-feature-overview]");
    const plot = node.querySelector("[data-nextclade-feature-plot]");
    const linkedTable = node.querySelector("[data-nextclade-linked-mutation-table]");
    if (!(select instanceof HTMLSelectElement) || !geneTable || !overview || !plot || !linkedTable) return;

    const getNtAaLinks = (ntItem, candidateAaRows, scopeGene = "") => {
      const start = Number(ntItem?.start);
      const end = Number(ntItem?.end || start);
      if (!Number.isFinite(start)) return [];
      const links = [];
      const genes = Array.isArray(ntItem?.genes) ? ntItem.genes : [];
      const scopedGenes = scopeGene && scopeGene !== "__genome__"
        ? genes.filter((gene) => gene === scopeGene)
        : genes;
      scopedGenes.forEach((gene) => {
        const geneRange = SARS_COV_2_GENE_RANGES.find((item) => item.gene === gene);
        if (!geneRange) return;
        const overlapStart = Math.max(start, geneRange.start);
        const overlapEnd = Math.min(Number.isFinite(end) ? end : start, geneRange.end);
        if (overlapStart > overlapEnd) return;
        const aaStart = Math.max(1, Math.floor((overlapStart - geneRange.start) / 3) + 1);
        const aaEnd = Math.max(aaStart, Math.floor((overlapEnd - geneRange.start) / 3) + 1);
        const matched = candidateAaRows.filter((aaItem) => (
          aaItem.gene === gene
          && Number.isFinite(Number(aaItem.position))
          && Number(aaItem.position) >= aaStart
          && Number(aaItem.position) <= aaEnd
        ));
        if (matched.length) {
          matched.forEach((aaItem) => links.push({
            gene,
            aaPosition: aaItem.position,
            aaMutation: aaItem.label,
          }));
        } else {
          links.push({
            gene,
            aaPosition: aaStart === aaEnd ? aaStart : `${aaStart}-${aaEnd}`,
            aaMutation: "--",
          });
        }
      });
      return links;
    };

    const summarizeKnowledgeItems = (items) => {
      if (!Array.isArray(items) || !items.length) return "--";
      const parts = [];
      items.forEach((item) => {
        const snippets = [
          item?.effect_zh || item?.effect || "",
          item?.detail_zh || item?.detail || "",
        ].map((value) => String(value || "").trim()).filter(Boolean);
        const text = snippets.join(" | ");
        if (text && !parts.includes(text)) parts.push(text);
      });
      return parts.length ? parts.join("；") : "--";
    };

    const getKnowledgeSummary = (aaMutations) => {
      const knowledgeMap = parsed?.knowledge?.aa_matches && typeof parsed.knowledge.aa_matches === "object"
        ? parsed.knowledge.aa_matches
        : {};
      const labels = Array.isArray(aaMutations) ? aaMutations.filter(Boolean) : [];
      const matchedRecords = [];
      const seen = new Set();
      labels.forEach((label) => {
        const items = Array.isArray(knowledgeMap[label]) ? knowledgeMap[label] : [];
        items.forEach((item) => {
          const key = [
            item?.section_key || "",
            item?.gene || "",
            item?.mutation || "",
            item?.effect || "",
            item?.detail || "",
          ].join("::");
          if (seen.has(key)) return;
          seen.add(key);
          matchedRecords.push(item);
        });
      });
      if (!matchedRecords.length) {
        return {
          categories: "--",
          effects: "--",
          details: "--",
          pmids: "--",
          categoriesMap: Object.fromEntries(NEXTCLADE_KNOWLEDGE_CATEGORY_COLUMNS.map((item) => [item.key, "--"])),
        };
      }
      const categories = [...new Set(matchedRecords.map((item) => item?.section_label).filter(Boolean))];
      const effects = [...new Set(matchedRecords.map((item) => item?.effect_zh || item?.effect).filter(Boolean))];
      const details = [...new Set(matchedRecords.map((item) => item?.detail_zh || item?.detail).filter(Boolean))];
      const pmids = [...new Set(
        matchedRecords.flatMap((item) => String(item?.pmid || "").split(/[;,]/).map((token) => token.trim()).filter(Boolean)),
      )];
      const categoriesMap = Object.fromEntries(NEXTCLADE_KNOWLEDGE_CATEGORY_COLUMNS.map((item) => {
        const scoped = matchedRecords.filter((record) => String(record?.section_key || "") === item.key);
        return [item.key, summarizeKnowledgeItems(scoped)];
      }));
      return {
        categories: categories.length ? categories.join("；") : "--",
        effects: effects.length ? effects.join("；") : "--",
        details: details.length ? details.join("；") : "--",
        pmids: pmids.length ? pmids.join("；") : "--",
        categoriesMap,
      };
    };

    const getKnowledgeHitCountForGene = (gene) => {
      if (!gene) return 0;
      const knowledgeMap = parsed?.knowledge?.aa_matches && typeof parsed.knowledge.aa_matches === "object"
        ? parsed.knowledge.aa_matches
        : {};
      const aaRows = Array.isArray(parsed?.aminoacid) ? parsed.aminoacid : [];
      const geneLabels = aaRows
        .filter((item) => item?.gene === gene && item?.label)
        .map((item) => String(item.label));
      return geneLabels.filter((label) => Array.isArray(knowledgeMap[label]) && knowledgeMap[label].length > 0).length;
    };

    const renderSelection = () => {
      hideAllNextcladeTooltips(node);
      const selected = select.value || "__genome__";
      const selectedLabel = select.options[select.selectedIndex]?.textContent || "Nucleotide sequence";
      const nucRows = selected === "__genome__"
        ? parsed.nucleotide
        : parsed.nucleotide.filter((item) => Array.isArray(item.genes) && item.genes.includes(selected));
      const aaRows = selected === "__genome__"
        ? parsed.aminoacid
        : parsed.aminoacid.filter((item) => item.gene === selected);
      const nucCounts = {
        total: nucRows.length,
        substitutions: nucRows.filter((item) => item.type === "Substitution").length,
        deletions: nucRows.filter((item) => item.type === "Deletion").length,
        insertions: nucRows.filter((item) => item.type === "Insertion").length,
      };
      const aaCounts = {
        total: aaRows.length,
        substitutions: aaRows.filter((item) => item.type === "Substitution").length,
        deletions: aaRows.filter((item) => item.type === "Deletion").length,
        insertions: aaRows.filter((item) => item.type === "Insertion").length,
      };
      const genes = parsed.genes || [];
      const selectedGeneMeta = genes.find((item) => item.gene === selected);
      const selectedRange = SARS_COV_2_GENE_RANGES.find((item) => item.gene === selected);
      const frameShiftRows = selected === "__genome__"
        ? (parsed.frameShifts || [])
        : (parsed.frameShifts || []).filter((item) => Array.isArray(item.genes) && item.genes.includes(selected));
      const geneSummaryMarkup = `
        <div class="table-frame nextclade-feature-summary-frame">
          <table class="report-table nextclade-feature-summary-table">
            <thead>
              <tr>
                <th>基因</th>
                <th>核苷酸突变</th>
                <th>氨基酸突变</th>
                <th>移码突变</th>
                <th>知识注释位点</th>
                <th>合计</th>
              </tr>
            </thead>
            <tbody>
              ${genes.map((item) => `
                <tr class="${item.gene === selected ? "is-active" : ""}" data-nextclade-feature-row="${escapeHtml(item.gene)}">
                  <td>${escapeHtml(item.gene)}</td>
                  <td>${escapeHtml(String(item.nucleotide || 0))}</td>
                  <td>${escapeHtml(String(item.aminoacid || 0))}</td>
                  <td>${escapeHtml(String(item.frameshift || 0))}</td>
                  <td>${escapeHtml(String(item.knowledgeHits || 0))}</td>
                  <td>${escapeHtml(String(item.total || 0))}</td>
                </tr>
              `).join("")}
            </tbody>
          </table>
        </div>
      `;
      geneTable.innerHTML = geneSummaryMarkup;
      overview.innerHTML = `
        <div class="nextclade-feature-headline">
          <div>
            <span class="section-chip">Genetic Feature</span>
            <h3>${escapeHtml(selectedLabel)}</h3>
            <p>${selected === "__genome__" ? "查看全基因组核苷酸突变与全部氨基酸突变。" : `查看 ${escapeHtml(selected)} 区域的核苷酸和氨基酸突变摘要。`}</p>
          </div>
          <div class="nextclade-feature-stats">
            <article class="nextclade-feature-stat">
              <span>Nucleotide</span>
              <strong>${escapeHtml(String(nucCounts.total))}</strong>
              <small>${escapeHtml(`${nucCounts.substitutions} 替换 / ${nucCounts.deletions} 缺失 / ${nucCounts.insertions} 插入 / ${frameShiftRows.length} 移码`)}</small>
            </article>
            <article class="nextclade-feature-stat">
              <span>Amino acid</span>
              <strong>${escapeHtml(String(aaCounts.total))}</strong>
              <small>${escapeHtml(`${aaCounts.substitutions} 替换 / ${aaCounts.deletions} 缺失 / ${aaCounts.insertions} 插入`)}</small>
            </article>
          </div>
        </div>
      `;
      const plotInnerWidth = 1000;
      const leftPad = 140;
      const geneStart = selectedRange?.start || 1;
      const geneEnd = selectedRange?.end || Math.max(...(parsed.nucleotide || []).map((item) => item.end || item.start || 1), 1);
      const geneLength = Math.max(1, geneEnd - geneStart + 1);
      const ntSpan = Math.max(1, geneEnd - geneStart);
      const mapNt = (value) => leftPad + (((value - geneStart) / ntSpan) * plotInnerWidth);
      const tickCount = 6;
      const tickValues = Array.from({ length: tickCount }, (_, index) => {
        if (index === tickCount - 1) return geneEnd;
        return Math.round(geneStart + ((geneLength - 1) * index) / (tickCount - 1));
      });
      const guideMarkup = tickValues.map((value, index) => {
        const x = mapNt(value);
        return `
          <g class="nextclade-plot-guide">
            <line x1="${x}" y1="46" x2="${x}" y2="168" class="nextclade-plot-guide-line"></line>
            <text x="${x}" y="34" text-anchor="${index === 0 ? "start" : (index === tickCount - 1 ? "end" : "middle")}" class="nextclade-plot-tick">${escapeHtml(String(value))}</text>
          </g>
        `;
      }).join("");
      const tooltipData = (parts) => encodeURIComponent(JSON.stringify(parts));
      const renderNtMarkers = nucRows.map((item, index) => {
        const start = Number(item.start || geneStart);
        const end = Number(item.end || start);
        const x = mapNt(start);
        const width = Math.max(6, end > start ? (mapNt(end) - x) : 6);
        const klass = item.type === "Deletion" ? "deletion" : (item.type === "Insertion" ? "insertion" : "substitution");
        const links = getNtAaLinks(item, aaRows, selected);
        const aaMutations = links.map((link) => link.aaMutation).filter((value) => value && value !== "--");
        const aaPositions = links.map((link) => `${link.gene}:${link.aaPosition}`).filter(Boolean);
        const payload = tooltipData({
          track: "Linked mutation",
          type: item.type,
          mutation: item.mutation,
          position: item.position_label,
          gene: selected === "__genome__" ? (item.genes || []).join(", ") : selected,
          aminoacid: aaMutations.length ? [...new Set(aaMutations)].join(", ") : "--",
          aaPosition: aaPositions.length ? [...new Set(aaPositions)].join(", ") : "--",
        });
        const aaBridge = links.length
          ? `<line x1="${x + Math.min(width / 2, 8)}" y1="88" x2="${x + Math.min(width / 2, 8)}" y2="132" class="nextclade-plot-link-line"></line><circle cx="${x + Math.min(width / 2, 8)}" cy="138" r="${aaMutations.length ? 5 : 3}" class="nextclade-plot-aa-dot ${aaMutations.length ? "has-aa" : ""}"></circle>`
          : "";
        if (item.type === "Insertion") {
          return `<g class="nextclade-plot-marker is-${klass}" data-nextclade-tooltip="${payload}" tabindex="0"><title>${escapeHtml(`${item.type}: ${item.mutation}`)}</title><g transform="translate(${x},74)"><path d="M0 0 L8 12 L-8 12 Z"></path></g>${aaBridge}</g>`;
        }
        return `<g class="nextclade-plot-marker is-${klass}" data-nextclade-tooltip="${payload}" tabindex="0"><title>${escapeHtml(`${item.type}: ${item.mutation}`)}</title><rect x="${x}" y="68" width="${width}" height="18" rx="4"></rect>${aaBridge}</g>`;
      }).join("");
      const renderFsMarkers = frameShiftRows.map((item) => {
        const start = Number(item.start || geneStart);
        const x = mapNt(start);
        const payload = tooltipData({
          track: "Frame shift",
          type: item.type,
          mutation: item.mutation,
          position: item.position_label,
          gene: (item.genes || []).join(", "),
        });
        return `<g class="nextclade-plot-marker is-frameshift" data-nextclade-tooltip="${payload}" tabindex="0"><title>${escapeHtml(`Frame shift: ${item.mutation}`)}</title><rect x="${x - 4}" y="112" width="8" height="28" rx="3"></rect></g>`;
      }).join("");
      plot.innerHTML = `
        <div class="nextclade-feature-plot-shell">
          <div class="nextclade-feature-plot-head">
            <div>
              <span class="section-chip">Interactive plot</span>
              <h3>${escapeHtml(selected === "__genome__" ? "全基因组 mutation track" : `${selected} mutation track`)}</h3>
            </div>
            <div class="nextclade-feature-legend">
              <span class="is-substitution">核苷酸替换</span>
              <span class="is-deletion">缺失</span>
              <span class="is-insertion">插入</span>
              <span class="is-aa">关联氨基酸</span>
              <span class="is-frameshift">移码</span>
            </div>
          </div>
          <div class="nextclade-feature-plot-frame">
            <svg viewBox="0 0 1180 180" class="nextclade-feature-plot" role="img" aria-label="${escapeHtml(selectedLabel)} mutation track">
              <rect x="20" y="0" width="${leftPad - 32}" height="180" class="nextclade-plot-label-rail"></rect>
              <line x1="${leftPad - 16}" y1="20" x2="${leftPad - 16}" y2="166" class="nextclade-plot-label-divider"></line>
              <rect x="${leftPad}" y="58" width="${plotInnerWidth}" height="88" rx="10" class="nextclade-plot-band nextclade-plot-band-linked"></rect>
              <text x="34" y="82" class="nextclade-plot-label">NT to AA</text>
              <text x="34" y="124" class="nextclade-plot-label nextclade-plot-label-small">Frame shift</text>
              <text x="${leftPad}" y="48" class="nextclade-plot-track-title">${escapeHtml(selected === "__genome__" ? "Genome linked track" : `${selected} linked track`)}</text>
              <line x1="${leftPad}" y1="78" x2="${leftPad + plotInnerWidth}" y2="78" class="nextclade-plot-axis"></line>
              <line x1="${leftPad}" y1="138" x2="${leftPad + plotInnerWidth}" y2="138" class="nextclade-plot-axis nextclade-plot-axis-aa"></line>
              ${guideMarkup}
              ${renderNtMarkers}
              ${renderFsMarkers}
            </svg>
            <div class="nextclade-plot-tooltip" data-nextclade-plot-tooltip hidden></div>
          </div>
        </div>
      `;
      const tooltip = plot.querySelector("[data-nextclade-plot-tooltip]");
      const plotFrame = plot.querySelector(".nextclade-feature-plot-frame");
      if (tooltip instanceof HTMLElement && plotFrame instanceof HTMLElement) {
        const hideTooltip = () => {
          tooltip.hidden = true;
        };
        const isInsideTrack = (event) => {
          const rect = plotFrame.getBoundingClientRect();
          return (
            event.clientX >= rect.left
            && event.clientX <= rect.right
            && event.clientY >= rect.top
            && event.clientY <= rect.bottom
          );
        };
        const hideTooltipOutsideTrack = (event) => {
          if (!isInsideTrack(event)) hideTooltip();
        };
        const hideTooltipWhenPointerTargetsOutsideTrack = (event) => {
          const target = event.target instanceof Element ? event.target : null;
          if (!target?.closest(".nextclade-feature-plot-frame")) hideTooltip();
        };
        const showTooltip = (event, marker) => {
          const raw = marker.getAttribute("data-nextclade-tooltip");
          if (!raw) return;
          let data = null;
          try {
            data = JSON.parse(decodeURIComponent(raw));
          } catch (_error) {
            return;
          }
          tooltip.innerHTML = `
            <span class="nextclade-tooltip-kicker">${escapeHtml(data.track || "Track")}</span>
            <strong>${escapeHtml(data.mutation || "--")}</strong>
            <span><b>类型</b>${escapeHtml(translateNextcladeMutationType(data.type))}</span>
            <span><b>核苷酸位置</b>${escapeHtml(data.position || "--")}</span>
            <span><b>基因</b>${escapeHtml(data.gene || "--")}</span>
            <span><b>氨基酸位置</b>${escapeHtml(data.aaPosition || "--")}</span>
            <span><b>关联氨基酸突变</b>${escapeHtml(data.aminoacid || "--")}</span>
          `;
          tooltip.hidden = false;
          const tooltipRect = tooltip.getBoundingClientRect();
          const margin = 14;
          const preferredX = event.clientX + 14;
          const preferredY = event.clientY + 14;
          const x = Math.min(preferredX, window.innerWidth - tooltipRect.width - margin);
          const y = Math.min(preferredY, window.innerHeight - tooltipRect.height - margin);
          tooltip.style.left = `${x}px`;
          tooltip.style.top = `${y}px`;
        };
        plotFrame.addEventListener("pointerleave", hideTooltip);
        plotFrame.addEventListener("mouseleave", hideTooltip);
        plotFrame.addEventListener("pointerout", (event) => {
          if (!plotFrame.contains(event.relatedTarget)) hideTooltip();
        });
        plotFrame.addEventListener("scroll", hideTooltip, { passive: true });
        plot.addEventListener("mouseleave", hideTooltip);
        document.addEventListener("pointermove", hideTooltipOutsideTrack, { capture: true, passive: true });
        document.addEventListener("pointerover", hideTooltipWhenPointerTargetsOutsideTrack, { capture: true, passive: true });
        document.addEventListener("wheel", hideTooltip, { capture: true, passive: true });
        document.addEventListener("scroll", hideTooltip, { capture: true, passive: true });
        plot.querySelectorAll("[data-nextclade-tooltip]").forEach((marker) => {
          marker.addEventListener("pointerenter", (event) => showTooltip(event, marker));
          marker.addEventListener("pointermove", (event) => showTooltip(event, marker));
          marker.addEventListener("focus", () => {
            const rect = marker.getBoundingClientRect();
            showTooltip({ clientX: rect.left, clientY: rect.top }, marker);
          });
          marker.addEventListener("blur", hideTooltip);
        });
      }
      const linkedColumns = [
        "核苷酸突变",
        "核苷酸类型",
        "核苷酸位置",
        "基因",
        "氨基酸位置",
        "关联氨基酸突变",
        ...NEXTCLADE_KNOWLEDGE_CATEGORY_COLUMNS.map((item) => item.label),
        "PMID",
      ];
      const linkedRows = nucRows.map((item) => {
          const links = getNtAaLinks(item, aaRows, selected);
          const displayGenes = selected === "__genome__"
            ? (item.genes || [])
            : (item.genes || []).filter((gene) => gene === selected);
          const aaPositions = links.map((link) => `${link.gene}:${link.aaPosition}`).filter(Boolean);
          const aaMutations = links.map((link) => link.aaMutation).filter((value) => value && value !== "--");
          const knowledge = getKnowledgeSummary(aaMutations);
          return [
            item.mutation,
            translateNextcladeMutationType(item.type),
            item.position_label,
            displayGenes.length ? displayGenes.join(", ") : "--",
            aaPositions.length ? [...new Set(aaPositions)].join(", ") : "--",
            aaMutations.length ? [...new Set(aaMutations)].join(", ") : "--",
            ...NEXTCLADE_KNOWLEDGE_CATEGORY_COLUMNS.map((category) => knowledge.categoriesMap?.[category.key] || "--"),
            knowledge.pmids,
          ];
        });
      linkedTable.dataset.exportTitle = `${selectedLabel}_突变位点与知识库关联表`;
      renderInteractiveContigTable(linkedTable, linkedColumns, linkedRows, "nextclade-linked-mutation-table");
    };

    select.addEventListener("change", renderSelection);
    window.addEventListener("scroll", () => hideAllNextcladeTooltips(node), { passive: true });
    geneTable.addEventListener("click", (event) => {
      const row = event.target instanceof Element ? event.target.closest("[data-nextclade-feature-row]") : null;
      const gene = row?.getAttribute("data-nextclade-feature-row");
      if (!gene || select.value === gene) return;
      select.value = gene;
      renderSelection();
    });
    renderSelection();
  });
}

function getChartExportFilename(container) {
  const title = container.querySelector('.mini-chart-title')?.textContent?.trim()
    || container.closest('.result-card')?.querySelector('.card-title-stack h3, .card-head h3')?.textContent?.trim()
    || 'chart';
  return slugifyFilename(title);
}

function serializeChartSvg(svg) {
  const clone = svg.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  if (!clone.getAttribute('viewBox')) {
    const width = svg.viewBox?.baseVal?.width || svg.clientWidth || 1200;
    const height = svg.viewBox?.baseVal?.height || svg.clientHeight || 480;
    clone.setAttribute('viewBox', `0 0 ${width} ${height}`);
  }
  return new XMLSerializer().serializeToString(clone);
}

async function requestChartExport(container, format) {
  const svg = container.querySelector('svg');
  if (!svg) return;
  const filename = getChartExportFilename(container);
  const svgText = serializeChartSvg(svg);
  const svgBlob = new Blob([svgText], { type: 'image/svg+xml;charset=utf-8' });
  if (format === 'svg') {
    downloadBlob(svgBlob, `${filename}.svg`);
    return;
  }
  if (format === 'png') {
    const url = URL.createObjectURL(svgBlob);
    try {
      const img = new Image();
      img.decoding = 'async';
      await new Promise((resolve, reject) => {
        img.onload = resolve;
        img.onerror = reject;
        img.src = url;
      });
      const viewBox = svg.viewBox?.baseVal;
      const width = Math.max(1, Math.round(viewBox?.width || svg.clientWidth || img.width || 1200));
      const height = Math.max(1, Math.round(viewBox?.height || svg.clientHeight || img.height || 480));
      const canvas = document.createElement('canvas');
      canvas.width = width * 2;
      canvas.height = height * 2;
      const ctx = canvas.getContext('2d');
      if (!ctx) throw new Error('图表导出失败');
      ctx.scale(2, 2);
      ctx.fillStyle = '#fbf8f1';
      ctx.fillRect(0, 0, width, height);
      ctx.drawImage(img, 0, 0, width, height);
      const pngBlob = await new Promise((resolve) => canvas.toBlob(resolve, 'image/png'));
      if (!pngBlob) throw new Error('PNG 导出失败');
      downloadBlob(pngBlob, `${filename}.png`);
    } finally {
      URL.revokeObjectURL(url);
    }
    return;
  }
  const title = container.querySelector('.mini-chart-title')?.textContent?.trim() || '图表';
  const printWindow = window.open('', '_blank', 'noopener,noreferrer');
  if (!printWindow) return;
  printWindow.document.write(`<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>${escapeHtml(title)}</title><style>
    @page { size: A4 landscape; margin: 12mm; }
    body { margin: 0; font-family: "Source Sans 3", "PingFang SC", sans-serif; background: #fffaf2; color: #243447; }
    .sheet { padding: 16px 18px; }
    h1 { margin: 0 0 10px; font-size: 20px; }
    p { margin: 0 0 16px; color: #526170; }
    .frame { border: 1px solid rgba(36,52,71,0.12); border-radius: 16px; padding: 12px; background: #fffdf8; }
    svg { width: 100%; height: auto; display: block; }
  </style></head><body><div class="sheet"><h1>${escapeHtml(title)}</h1><p>图表导出副本</p><div class="frame">${svgText}</div></div><script>window.addEventListener('load',function(){setTimeout(function(){window.print();},300);});</script></body></html>`);
  printWindow.document.close();
}

function bindChartExportButtons(root = document) {
  root.querySelectorAll('.mini-chart-card').forEach((card) => {
    if (card.dataset.chartEnhanced === 'true') return;
    if (!card.querySelector('svg')) return;
    const title = card.querySelector('.mini-chart-title');
    if (!title) return;
    const bar = document.createElement('div');
    bar.className = 'chart-titlebar';
    title.replaceWith(bar);
    bar.appendChild(title);
    bar.insertAdjacentHTML('beforeend', renderChartExportToolbar());
    card.dataset.chartEnhanced = 'true';
  });
  root.querySelectorAll('[data-chart-export-format]').forEach((button) => {
    if (button.dataset.bound === 'true') return;
    button.dataset.bound = 'true';
    button.addEventListener('click', async (event) => {
      event.stopPropagation();
      try {
        await requestChartExport(button.closest('.mini-chart-card'), button.dataset.chartExportFormat || 'svg');
      } catch (error) {
        console.error(error);
      }
    });
  });
}

function renderInteractiveTableSummary(columns, columnKinds, state, totalRows, filteredCount) {
  const activeFilters = columns
    .map((column, index) => ({ column, value: String(state.filters[index] || "").trim(), kind: columnKinds[index], index }))
    .filter((item) => item.value);
  const sortLabel = state.sortIndex >= 0
    ? `${columns[state.sortIndex]} · ${state.sortDirection === "asc" ? "升序" : "降序"}`
    : "未排序";
  return `
    <div class="table-statebar">
      <div class="table-state-toprow">
        <div class="table-state-meta">
          <span class="table-state-pill">显示 ${filteredCount} / ${totalRows} 行</span>
          <span class="table-state-pill">${escapeHtml(sortLabel)}</span>
          ${activeFilters.length ? `<span class="table-state-pill">已筛选 ${activeFilters.length} 列</span>` : ""}
        </div>
        <div class="table-view-switch" role="tablist" aria-label="表格视图切换">
          <button class="table-view-button ${state.viewMode === "key" ? "active" : ""}" type="button" data-table-view="key">关键列</button>
          <button class="table-view-button ${state.viewMode === "all" ? "active" : ""}" type="button" data-table-view="all">全部列</button>
        </div>
      </div>
      <div class="table-filter-chip-row">
        ${activeFilters.length ? activeFilters.map((item) => `
          <button class="table-filter-chip" type="button" data-clear-filter-index="${item.index}">
            <span>${escapeHtml(item.column)}</span>
            <strong>${escapeHtml(item.value)}</strong>
          </button>
        `).join("") : `<span class="table-filter-placeholder">当前未设置筛选条件</span>`}
      </div>
    </div>
  `;
}

function getInteractiveTableSpec(tableId, columns) {
  const aliasSets = {
    "patho-mlst-table": ["样本名称", "MLST方案", "ST", "adk", "fumC", "glyA", "tyrB", "icd", "pepA", "pgm"],
    "patho-mutate-table": ["染色体", "变异位点位置", "参考位点碱基", "突变样本数", "替代碱基统计", "发生突变的样本"],
    "binning-quality-table": ["Bin名称", "质量等级", "完整性", "污染率", "基因组大小", "Contig总数", "Contig N50"],
    "binning-taxonomy-table": ["Bin名称", "门", "属", "种", "参考基因组", "ANI", "AF", "分类方法"],
    "assembly-species-table": ["序列名称", "基因组/质粒", "平均深度", "taxid", "物种名称", "NCBI TaxID", "NCBI学名", "NCBI分类等级", "属", "种", "NCBI属", "NCBI种", "科", "NCBI科", "目", "NCBI目", "界", "门"],
    "contig-annotation-table": ["seq_name", "length", "cov.", "circ.", "repeat", "mult.", "alt_group"],
    "rv-summary-table": ["Bin名称", "序列名称", "物种名称", "平均深度", "基因组/质粒", "毒力基因", "耐药基因", "血清型"],
    "virulence-table": ["Contig名称", "物种名称", "基因名称", "VF分类", "覆盖度%", "一致性%", "平均深度"],
    "resistance-table": ["Contig名称", "物种名称", "基因名称", "耐药药物", "覆盖度%", "一致性%", "平均深度"],
    "mge-resistance-table": ["样本名称", "基因名称", "所在序列", "关联元件类型", "转移风险等级", "位于预测 MGE 边界", "最近核心模块", "最近距离(bp)"],
    "mge-virulence-table": ["样本名称", "基因名称", "所在序列", "关联元件类型", "转移风险等级", "位于预测 MGE 边界", "最近核心模块", "最近距离(bp)"],
    "nextclade-linked-mutation-table": ["核苷酸突变", "核苷酸类型", "核苷酸位置", "基因", "氨基酸位置", "关联氨基酸突变"],
    "nextclade-variant-annotation-table": ["基因名", "位置", "参考碱基", "突变碱基", "测序深度 / 突变频率", "HGVS.c", "HGVS.p"],
    "rsv-typing-mutation-table": ["基因名", "位置", "参考碱基", "突变碱基", "测序深度 / 突变频率", "HGVS.c", "HGVS.p"],
    "rsv-typing-nmdc-table": ["亚型", "基因", "突变位点", "基因组位置", "我们的突变结果", "质量分层", "抗体亲和", "氨基酸突变"],
  };
  const preferred = aliasSets[tableId] || [];
  const normalized = columns.map((column) => String(column || "").toLowerCase());
  const keyIndexes = preferred
    .map((alias) => normalized.findIndex((column) => column.includes(String(alias).toLowerCase())))
    .filter((index, position, list) => index >= 0 && list.indexOf(index) === position);
  return {
    keyIndexes: keyIndexes.length ? keyIndexes : columns.slice(0, Math.min(columns.length, 6)).map((_, index) => index),
    stickyFirstColumn: !["nextclade-linked-mutation-table", "nextclade-variant-annotation-table"].includes(tableId),
  };
}

function getTableCellTone(column, value) {
  const text = String(value ?? "-").trim();
  const displayText = formatReportTableDisplayValue(column, value);
  const numeric = extractNumericValue(text);
  if (!text || text === "-") return { tone: "", render: `<span>${escapeHtml(text || "-")}</span>` };
  if (column === "ST") {
    return { tone: "tag", render: `<span class="table-inline-tag">${escapeHtml(`ST ${text}`)}</span>` };
  }
  if (["基因组/质粒", "是否成环", "VF分类", "耐药药物", "物种名称", "属", "种"].includes(column)) {
    return { tone: "tag", render: `<span class="table-inline-tag">${escapeHtml(text)}</span>` };
  }
  if (column.includes("完整性")) {
    const tone = getCompletenessState(numeric);
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(displayText)}</span>` };
  }
  if (column.includes("污染率")) {
    const tone = getContaminationState(numeric);
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(displayText)}</span>` };
  }
  if (column.includes("覆盖度") && column.includes("%")) {
    const tone = numeric >= 95 ? "success" : numeric >= 80 ? "warning" : "danger";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(displayText)}</span>` };
  }
  if (column.includes("一致性")) {
    const tone = numeric >= 95 ? "success" : numeric >= 85 ? "warning" : "danger";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(displayText)}</span>` };
  }
  if (column.includes("平均深度")) {
    const tone = numeric >= 30 ? "success" : numeric >= 10 ? "warning" : "danger";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone}">${escapeHtml(displayText)}</span>` };
  }
  if (column === "质量分层") {
    const tone = text.includes("高质量") ? "success" : text.includes("低质量") ? "warning" : "";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone || "neutral"}">${escapeHtml(text)}</span>` };
  }
  if (column === "影响等级") {
    const toneMap = { HIGH: "danger", MODERATE: "warning", LOW: "success", MODIFIER: "neutral" };
    const tone = toneMap[text] || "";
    return { tone, render: `<span class="table-inline-chip table-inline-chip-${tone || "neutral"}">${escapeHtml(text)}</span>` };
  }
  if (column === "变异类型") {
    return { tone: "tag", render: `<span class="table-inline-tag">${escapeHtml(text)}</span>` };
  }
  return { tone: "", render: renderTableCellContent(value, column) };
}

async function requestTableExport(container, format) {
  const shell = document.querySelector(".report-shell");
  const payload = tableExportStore.get(container);
  if (!shell || !payload) return;
  const response = await fetch(shell.dataset.exportEndpoint, {
    method: "POST",
    credentials: "same-origin",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      title: payload.title,
      filename: slugifyFilename(payload.title),
      columns: payload.columns,
      rows: payload.rows,
      format,
    }),
  });
  if (!response.ok) {
    const data = await response.json().catch(() => ({}));
    throw new Error(data.error || "表格导出失败");
  }
  const blob = await response.blob();
  downloadBlob(blob, `${slugifyFilename(payload.title)}.${format}`);
}

function bindTableExportButtons(container, title, columns, rows) {
  tableExportStore.set(container, { title, columns, rows });
  container.querySelectorAll("[data-table-export-format]").forEach((button) => {
    button.addEventListener("click", async () => {
      try {
        await requestTableExport(container, button.dataset.tableExportFormat || "csv");
      } catch (error) {
        console.error(error);
      }
    });
  });
}

async function fetchReportCssText() {
  const stylesheet = document.querySelector('link[href*="report.css"]');
  if (!stylesheet) return "";
  const response = await fetch(stylesheet.href, { credentials: "same-origin" });
  return response.ok ? response.text() : "";
}

async function fetchReportJsText() {
  const runtimeScript = document.querySelector('script[src*="report_runtime.js"]');
  if (!runtimeScript?.src) return "";
  const response = await fetch(runtimeScript.src, { credentials: "same-origin" });
  return response.ok ? response.text() : "";
}

async function buildEmbeddedInteractiveReportData(data) {
  const cloned = JSON.parse(JSON.stringify(data || {}));
  const serotype = cloned?.sections?.serotype;
  if (String(serotype?.mode || "").trim() !== "hiv_resistance") {
    return cloned;
  }
  const taskId = String(cloned?.task?.id || "").trim();
  const bootscanAssets = serotype?.bootscan_assets && typeof serotype.bootscan_assets === "object"
    ? serotype.bootscan_assets
    : {};
  const embedded = {};
  await Promise.all([
    ["overall_csv_text", "overall_csv"],
    ["pure_csv_text", "pure_csv"],
  ].map(async ([targetKey, sourceKey]) => {
    const assetName = String(bootscanAssets?.[sourceKey] || "").trim();
    if (!taskId || !assetName) return;
    try {
      embedded[targetKey] = await fetchReportAsset(taskId, assetName, "text");
    } catch (_error) {
      // 导出时尽量嵌入可离线重绘的数据；若单个资产失败，保留线上路径作为兜底。
    }
  }));
  if (Object.keys(embedded).length) {
    serotype.bootscan_embedded = embedded;
  }
  return cloned;
}

function isVirusFocusedReport(data) {
  return isSarsCov2NextcladeReport(data)
    || isHmpvNextcladeReport(data)
    || isDenvNextcladeReport(data)
    || isZikavNextcladeReport(data)
    || isChikvNextcladeReport(data)
    || isEbolaNextcladeReport(data)
    || isHpivTypingReport(data)
    || isHadvTypingReport(data)
    || isNorovirusTypingReport(data)
    || isEnterovirusTypingReport(data)
    || isHivTypingReport(data)
    || isHepatovirusTypingReport(data)
    || isBandavirusTypingReport(data)
    || isOrthohantavirusTypingReport(data)
    || isAstroviridaeTypingReport(data)
    || isRhinovirusTypingReport(data)
    || isSeasonalHcovTypingReport(data)
    || isRotavirusTypingReport(data)
    || isMonkeypoxNextcladeReport(data)
    || isInfluenzaTypingReport(data);
}

function pruneBacteriaOnlySectionsForVirus(root) {
  if (!root) return;
  [
    "section-contig-annotation",
    "section-cgview",
    "section-checkm",
    "section-gene-annotation",
    "section-rv",
    "section-rv-summary",
    "section-virulence",
    "section-resistance",
    "section-resistance-mutation",
    "section-mlst",
    "section-priority-serotype",
    "section-mge",
    "section-mge-resistance",
    "section-mge-virulence",
  ].forEach((id) => {
    root.querySelector?.(`#${id}`)?.remove();
  });
  [
    '#section-contig-annotation',
    '#section-cgview',
    '#section-checkm',
    '#section-gene-annotation',
    '#section-rv',
    '#section-rv-summary',
    '#section-virulence',
    '#section-resistance',
    '#section-resistance-mutation',
    '#section-mlst',
    '#section-priority-serotype',
    '#section-mge',
    '#section-mge-resistance',
    '#section-mge-virulence',
  ].forEach((href) => {
    root.querySelectorAll?.(`.report-nav-link[href="${href}"]`).forEach((node) => node.remove());
  });
}

function escapeInlineScriptText(value) {
  return String(value ?? "")
    .replace(/<\/script/gi, "<\\/script")
    .replace(/<!--/g, "<\\!--")
    .replace(/\u2028/g, "\\u2028")
    .replace(/\u2029/g, "\\u2029");
}

async function buildReportExportHtml(printMode = false, interactiveMode = false) {
  const clone = document.documentElement.cloneNode(true);
  const shell = clone.querySelector(".report-shell");
  const taskName = document.querySelector(".report-shell")?.dataset.taskName || document.getElementById("report-task-name")?.textContent || "analysis_report";
  const taskStatus = document.querySelector(".report-status")?.textContent?.trim() || "未知";
  const exportedAt = new Intl.DateTimeFormat("zh-CN", {
    year: "numeric",
    month: "2-digit",
    day: "2-digit",
    hour: "2-digit",
    minute: "2-digit",
    second: "2-digit",
    hour12: false,
  }).format(new Date());
  clone.body?.classList.add("report-export-document");
  if (printMode) clone.body?.classList.add("report-print-export");
  if (interactiveMode) clone.body?.classList.add("report-interactive-export");
  clone.querySelectorAll("script").forEach((node) => node.remove());
  clone.querySelectorAll(".report-export-toolbar, .table-export-toolbar, .report-back").forEach((node) => node.remove());
  if (!interactiveMode) {
    clone.querySelectorAll(".table-filter-input, .table-statebar, .table-view-switch").forEach((node) => node.remove());
    clone.querySelectorAll(".table-sort-button").forEach((button) => {
      const replacement = clone.ownerDocument.createElement("span");
      replacement.textContent = button.querySelector("span")?.textContent || button.textContent || "";
      button.replaceWith(replacement);
    });
    clone.querySelectorAll(".mlst-hostgene-button").forEach((button) => {
      const replacement = clone.ownerDocument.createElement("span");
      replacement.textContent = button.textContent || "";
      button.replaceWith(replacement);
    });
    clone.querySelectorAll("td span[title]").forEach((span) => {
      span.textContent = span.getAttribute("title") || span.textContent || "";
      span.removeAttribute("title");
    });
  }
  clone.querySelectorAll('link[href*="report.css"]').forEach((node) => node.remove());
  const cssText = await fetchReportCssText();
  if (cssText) {
    const style = clone.ownerDocument.createElement("style");
    style.textContent = cssText;
    clone.querySelector("head")?.appendChild(style);
  }
  if (shell) {
    const masthead = clone.ownerDocument.createElement("section");
    masthead.className = "report-document-masthead";
    masthead.innerHTML = `
      <div class="report-document-headline">
        <span class="report-document-label">${printMode ? "打印版报告" : interactiveMode ? "交互导出版" : "归档导出版"}</span>
        <h2>病原微生物分析结果归档副本</h2>
        <p>本副本用于结果存档、离线查看与纸质输出。表格与图形已整理为正式报告阅读顺序。</p>
      </div>
      <dl class="report-document-meta">
        <div><dt>任务名称</dt><dd>${escapeHtml(taskName)}</dd></div>
        <div><dt>任务状态</dt><dd>${escapeHtml(taskStatus)}</dd></div>
        <div><dt>导出时间</dt><dd>${escapeHtml(exportedAt)}</dd></div>
        <div><dt>导出形式</dt><dd>${printMode ? "PDF打印" : interactiveMode ? "HTML交互文档" : "Word归档文档"}</dd></div>
      </dl>
    `;
    shell.insertBefore(masthead, shell.querySelector(".report-layout"));
  }
  if (currentReportData && isVirusFocusedReport(currentReportData)) {
    pruneBacteriaOnlySectionsForVirus(clone);
  }
  if (interactiveMode && currentReportData) {
    if (shell) {
      shell.dataset.reportEndpoint = "";
      shell.dataset.exportEndpoint = "";
    }
    const embeddedReportData = await buildEmbeddedInteractiveReportData(currentReportData);
    const dataScript = clone.ownerDocument.createElement("script");
    dataScript.textContent = `window.__EMBEDDED_REPORT_DATA__ = ${escapeInlineScriptText(JSON.stringify(embeddedReportData))};`;
    clone.querySelector("body")?.appendChild(dataScript);
    const runtimeText = await fetchReportJsText();
    if (runtimeText) {
      const script = clone.ownerDocument.createElement("script");
      script.textContent = escapeInlineScriptText(runtimeText);
      clone.querySelector("body")?.appendChild(script);
    }
  }
  if (printMode) {
    const script = clone.ownerDocument.createElement("script");
    script.textContent = escapeInlineScriptText("window.addEventListener('load', function(){ setTimeout(function(){ window.print(); }, 300); });");
    clone.querySelector("body")?.appendChild(script);
  }
  return `<!doctype html>
${clone.outerHTML}`;
}

async function exportReportPage(format) {
  const shell = document.querySelector(".report-shell");
  const baseName = slugifyFilename(shell?.dataset.taskName || "analysis_report");
  if (format === "pdf") {
    const html = await buildReportExportHtml(true, false);
    const blob = new Blob([html], { type: "text/html;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    window.open(url, "_blank");
    window.setTimeout(() => URL.revokeObjectURL(url), 60000);
    return;
  }
  const html = await buildReportExportHtml(false, format === "html");
  if (format === "word") {
    downloadBlob(new Blob([html], { type: "application/msword;charset=utf-8" }), `${baseName}.doc`);
    return;
  }
  downloadBlob(new Blob([html], { type: "text/html;charset=utf-8" }), `${baseName}.html`);
}

function normalizeInteractiveTableRows(columns, rows) {
  const safeColumns = Array.isArray(columns) ? columns : [];
  const safeRows = Array.isArray(rows) ? rows : [];
  return safeRows.map((row) => {
    if (Array.isArray(row)) {
      if (row.length === safeColumns.length) return row;
      const normalized = new Array(safeColumns.length).fill("");
      safeColumns.forEach((_, index) => {
        normalized[index] = index < row.length ? row[index] : "";
      });
      return normalized;
    }
    if (row && typeof row === "object") {
      return safeColumns.map((column) => row[column] ?? "");
    }
    return safeColumns.map(() => "");
  });
}

const REPORT_TABLE_HEADER_LABELS = {
  sample: "样本",
  samples: "样本",
  sample_name: "样本名称",
  broad_type: "大分型",
  af_group: "A/F 分组",
  cj_group: "CJ 分组",
  typed_label: "最优亚型",
  segment_groups: "片段分组",
  reassortment_flag: "重配标记",
  selection_summary_path: "结果摘要路径",
  selection_path: "结果路径",
  summary_path: "摘要路径",
  file_path: "文件路径",
  broadtype: "大分型",
  msa_method: "比对方法",
  tree_method: "建树方法",
  tree: "系统树",
  ref_name: "参考名称",
  ref_accession: "参考登录号",
  ref_length: "参考长度",
  query_length: "查询长度",
  genome_length: "基因组长度",
  contig_count: "Contig 数",
  contig_length: "Contig 长度",
  gc_content: "GC 含量",
  n50: "N50",
  taxid: "TaxID",
  species: "物种",
  genus: "属",
  family: "科",
  order: "目",
  class: "纲",
  phylum: "门",
  lineage: "谱系",
  clade: "分支",
  genotype: "基因型",
  subtype: "亚型",
  serotype: "血清型",
  reads: "读段数",
  read_count: "读段数",
  support_reads: "支持读段数",
  coverage: "覆盖度",
  coverage_percent: "覆盖度(%)",
  depth: "深度",
  identity: "一致性",
  abundance: "丰度",
  frequency: "频率",
  score: "评分",
  status: "状态",
  result: "结果",
  source: "来源",
  source_info: "来源信息",
  accession: "登录号",
  reference: "参考",
  reference_id: "参考 ID",
  gene: "基因",
  genes: "基因",
  locus: "位点",
  position: "位置",
  mutation: "突变",
  mutations: "突变",
  host: "宿主",
  country: "国家",
  region: "地区",
  province: "省份",
  city: "城市",
  date: "日期",
  path: "路径",
  bootstrap: "自举值",
  distance: "距离",
  cluster: "聚类",
  host_gene: "宿主基因",
  host_gene_id: "宿主基因 ID",
  "host gene": "宿主基因",
  "host gene id": "宿主基因 ID",
  "host gene 展示": "宿主基因展示",
};

const REPORT_TABLE_HEADER_TOKENS = {
  sample: "样本",
  samples: "样本",
  broad: "大",
  type: "类型",
  af: "A/F",
  group: "分组",
  groups: "分组",
  segment: "片段",
  reassortment: "重配",
  flag: "标记",
  selection: "选择",
  summary: "摘要",
  path: "路径",
  read: "读段",
  reads: "读段",
  count: "数",
  support: "支持",
  coverage: "覆盖度",
  depth: "深度",
  identity: "一致性",
  abundance: "丰度",
  frequency: "频率",
  score: "评分",
  species: "物种",
  genus: "属",
  family: "科",
  order: "目",
  class: "纲",
  phylum: "门",
  lineage: "谱系",
  clade: "分支",
  genotype: "基因型",
  subtype: "亚型",
  serotype: "血清型",
  host: "宿主",
  gene: "基因",
  genes: "基因",
  locus: "位点",
  position: "位置",
  mutation: "突变",
  mutations: "突变",
  source: "来源",
  info: "信息",
  status: "状态",
  result: "结果",
  reference: "参考",
  ref: "参考",
  alt: "替代",
  alternate: "替代",
  accession: "登录号",
  bootstrap: "自举",
  distance: "距离",
  cluster: "聚类",
  country: "国家",
  region: "地区",
  province: "省份",
  city: "城市",
  date: "日期",
};

function humanizeReportColumnLabel(column) {
  const raw = String(column || "").trim();
  if (!raw) return "";
  const exact = REPORT_TABLE_HEADER_LABELS[raw] || REPORT_TABLE_HEADER_LABELS[raw.toLowerCase()];
  if (exact) return exact;
  if (/[\u4e00-\u9fa5]/.test(raw)) return raw;
  const normalized = raw
    .replace(/([a-z0-9])([A-Z])/g, "$1_$2")
    .replace(/[()/.-]+/g, " ")
    .replace(/\s+/g, "_")
    .trim()
    .toLowerCase();
  const mapped = normalized
    .split("_")
    .filter(Boolean)
    .map((token) => REPORT_TABLE_HEADER_TOKENS[token] || token.toUpperCase())
    .join("");
  return mapped || raw;
}

function renderInteractiveContigTable(container, columns, rows, tableId = "") {
  const normalizedColumns = Array.isArray(columns) ? columns : [];
  const normalizedRows = normalizeInteractiveTableRows(normalizedColumns, rows);
  const columnKinds = inferColumnKinds(normalizedColumns, normalizedRows);
  const spec = getInteractiveTableSpec(tableId, normalizedColumns);
  const isIgvLinkedTable = [
    "influenza-typing-mutation-table",
    "nextclade-variant-annotation-table",
    "monkeypox-typing-mutation-table",
    "rsv-typing-mutation-table",
  ].includes(tableId);
  const locusColumnIndex = isIgvLinkedTable ? normalizedColumns.indexOf("染色体") : -1;
  const posColumnIndex = isIgvLinkedTable ? normalizedColumns.indexOf("位置") : -1;
  const state = {
    sortIndex: -1,
    sortDirection: "asc",
    filters: normalizedColumns.map(() => ""),
    scrollLeft: 0,
    scrollTop: 0,
    viewMode: spec.keyIndexes?.length ? "key" : "all",
  };

  const applyState = (focusIndex = null, caretPosition = null) => {
    const previousFrame = container.querySelector(".table-frame");
    if (previousFrame) {
      state.scrollLeft = previousFrame.scrollLeft;
      state.scrollTop = previousFrame.scrollTop;
    }

    const filteredRows = normalizedRows.filter((row) => normalizedColumns.every((column, index) => {
      const keyword = state.filters[index]?.trim();
      if (!keyword) return true;
      const value = row[index];
      if (columnKinds[index] === "number") {
        return matchNumericFilter(value, keyword);
      }
      return String(value ?? "").toLowerCase().includes(keyword.toLowerCase());
    }));

    if (state.sortIndex >= 0) {
      filteredRows.sort((leftRow, rightRow) => {
        const left = leftRow[state.sortIndex];
        const right = rightRow[state.sortIndex];
        return compareTableValues(left, right, state.sortDirection);
      });
    }

    const visibleIndexes = state.viewMode === "key" && spec.keyIndexes?.length
      ? spec.keyIndexes
      : normalizedColumns.map((_, index) => index);
    const exportTitle = container.dataset.exportTitle || "结果表";
    container.innerHTML = `
      ${renderTableExportToolbar()}
      ${renderInteractiveTableSummary(normalizedColumns, columnKinds, state, normalizedRows.length, filteredRows.length)}
      <div class="table-scroll-assist" aria-label="表格横向滚动辅助">
        <button type="button" data-table-scroll="left">左移</button>
        <span data-table-scroll-status>横向位置 0%</span>
        <button type="button" data-table-scroll="right">右移</button>
      </div>
      <div class="table-frame interactive-table-frame tall-table-frame ${spec.stickyFirstColumn === false ? "table-no-sticky-first" : ""}">
        <table class="report-table report-table-interactive ${spec.stickyFirstColumn === false ? "table-no-sticky-first" : ""}">
          <thead>
            <tr>
              ${visibleIndexes.map((index) => {
                const column = normalizedColumns[index];
                return `
                <th>
                  <div class="table-head-stack">
                    <button class="table-sort-button" type="button" data-sort-index="${index}">
                      <span>${escapeHtml(humanizeReportColumnLabel(column))}</span>
                      <span class="table-sort-indicator">${state.sortIndex === index ? (state.sortDirection === "asc" ? "▲" : "▼") : "↕"}</span>
                    </button>
                  </div>
                  <div class="table-filter-wrap">
                    <input class="table-filter-input" type="text" placeholder="${columnKinds[index] === "number" ? "如 >100 或 50-200" : "筛选"}" value="${escapeHtml(state.filters[index] || "")}" data-filter-index="${index}">
                  </div>
                </th>`;
              }).join("")}
            </tr>
          </thead>
          <tbody>
            ${filteredRows.map((row) => {
              const locus = locusColumnIndex >= 0 && posColumnIndex >= 0
                ? `${String(row[locusColumnIndex] ?? "").trim()}:${String(row[posColumnIndex] ?? "").trim()}`
                : "";
              const rowAttr = locus && !locus.startsWith(":") && !locus.endsWith(":")
                ? ` data-igv-locus="${escapeHtml(locus)}" class="igv-locus-row" title="点击跳转到 IGV：${escapeHtml(locus)}"`
                : "";
              return `<tr${rowAttr}>${visibleIndexes.map((index) => {
              const column = normalizedColumns[index];
              const value = row[index];
              const cell = getTableCellTone(column, value);
              const cellAttr = isIgvLinkedTable && locus
                ? ` data-igv-locus="${escapeHtml(locus)}" title="点击跳转到 IGV：${escapeHtml(locus)}" onclick="window.__influenzaIgvSelectLocus && window.__influenzaIgvSelectLocus(this.dataset.igvLocus)" onpointerdown="window.__influenzaIgvSelectLocus && window.__influenzaIgvSelectLocus(this.dataset.igvLocus)" ontouchstart="window.__influenzaIgvSelectLocus && window.__influenzaIgvSelectLocus(this.dataset.igvLocus)"`
                : "";
              return `<td class="${cell.tone ? `table-cell-${cell.tone}` : ""}"${cellAttr}>${cell.render}</td>`;
            }).join("")}</tr>`;
            }).join("")}
          </tbody>
        </table>
      </div>
    `;
    bindTableExportButtons(container, exportTitle, normalizedColumns, filteredRows);

    container.querySelectorAll('.table-sort-button').forEach((button) => {
      button.addEventListener('click', () => {
        const index = Number(button.dataset.sortIndex);
        if (state.sortIndex === index) {
          state.sortDirection = state.sortDirection === 'asc' ? 'desc' : 'asc';
        } else {
          state.sortIndex = index;
          state.sortDirection = 'asc';
        }
        applyState();
      });
    });

    container.querySelectorAll('.table-filter-input').forEach((input) => {
      input.addEventListener('compositionstart', () => {
        input.dataset.imeComposing = '1';
      });
      input.addEventListener('compositionend', () => {
        input.dataset.imeComposing = '0';
        const index = Number(input.dataset.filterIndex);
        state.filters[index] = input.value;
        applyState(index, input.selectionStart ?? input.value.length);
      });
      input.addEventListener('input', (event) => {
        if (event.isComposing || input.dataset.imeComposing === '1') return;
        const index = Number(input.dataset.filterIndex);
        state.filters[index] = input.value;
        applyState(index, input.selectionStart ?? input.value.length);
      });
    });

    container.querySelectorAll('[data-clear-filter-index]').forEach((button) => {
      button.addEventListener('click', () => {
        const index = Number(button.dataset.clearFilterIndex);
        state.filters[index] = '';
        state.viewMode = state.viewMode || 'key';
        applyState(index, 0);
      });
    });

    container.querySelectorAll('[data-table-view]').forEach((button) => {
      button.addEventListener('click', () => {
        state.viewMode = button.dataset.tableView === 'all' ? 'all' : 'key';
        applyState();
      });
    });

    const nextFrame = container.querySelector(".table-frame");
    if (nextFrame) {
      const scrollStatus = container.querySelector("[data-table-scroll-status]");
      const syncScrollAssist = () => {
        if (!(scrollStatus instanceof HTMLElement)) return;
        const maxLeft = Math.max(0, nextFrame.scrollWidth - nextFrame.clientWidth);
        const percent = maxLeft ? Math.round((nextFrame.scrollLeft / maxLeft) * 100) : 0;
        scrollStatus.textContent = maxLeft ? `横向位置 ${percent}%` : "无需横向滚动";
      };
      nextFrame.scrollLeft = state.scrollLeft;
      nextFrame.scrollTop = state.scrollTop;
      nextFrame.addEventListener("scroll", () => {
        state.scrollLeft = nextFrame.scrollLeft;
        state.scrollTop = nextFrame.scrollTop;
        syncScrollAssist();
      }, { passive: true });
      container.querySelectorAll("[data-table-scroll]").forEach((button) => {
        button.addEventListener("click", () => {
          const direction = button.dataset.tableScroll === "left" ? -1 : 1;
          nextFrame.scrollBy({ left: direction * Math.max(260, nextFrame.clientWidth * 0.72), behavior: "smooth" });
        });
      });
      syncScrollAssist();
      enableInteractiveTableDragScroll(nextFrame);
    }

    if (focusIndex !== null) {
      const target = container.querySelector(`.table-filter-input[data-filter-index="${focusIndex}"]`);
      if (target) {
        target.focus();
        const nextPos = Math.min(caretPosition ?? target.value.length, target.value.length);
        target.setSelectionRange(nextPos, nextPos);
      }
    }
  };

  applyState();
}

function enableInteractiveTableDragScroll(frame) {
  if (!frame || frame.dataset.dragScrollBound === "true") return;
  frame.dataset.dragScrollBound = "true";
  let active = false;
  let startX = 0;
  let startY = 0;
  let startLeft = 0;
  let startTop = 0;

  frame.addEventListener("pointerdown", (event) => {
    if (event.pointerType === "mouse" && event.button !== 0) return;
    const target = event.target;
    if (target instanceof HTMLElement) {
      if (target.closest("input, button, a, select, textarea")) return;
    }
    active = true;
    startX = event.clientX;
    startY = event.clientY;
    startLeft = frame.scrollLeft;
    startTop = frame.scrollTop;
    frame.classList.add("is-dragging");
    try {
      frame.setPointerCapture(event.pointerId);
    } catch (_error) {
      // ignore unsupported capture states
    }
  });

  frame.addEventListener("pointermove", (event) => {
    if (!active) return;
    frame.scrollLeft = startLeft - (event.clientX - startX);
    frame.scrollTop = startTop - (event.clientY - startY);
  });

  const release = (event) => {
    if (!active) return;
    active = false;
    frame.classList.remove("is-dragging");
    try {
      frame.releasePointerCapture(event.pointerId);
    } catch (_error) {
      // ignore unsupported capture states
    }
  };

  frame.addEventListener("pointerup", release);
  frame.addEventListener("pointercancel", release);
  frame.addEventListener("lostpointercapture", () => {
    active = false;
    frame.classList.remove("is-dragging");
  });
}

function applyTableTone(container, containerId) {
  const tones = {
    "binning-quality-table": "quality",
    "binning-taxonomy-table": "taxonomy",
    "assembly-species-table": "taxonomy",
    "contig-annotation-table": "assembly",
    "assembly-summary-table": "assembly",
    "checkm-table": "quality",
    "gene-annotation-summary-table": "gene",
    "rv-summary-table": "risk",
    "virulence-table": "virulence",
    "resistance-table": "resistance",
    "mlst-table": "typing",
    "serotype-table": "typing",
    "priority-serotype-table": "risk",
    "mge-resistance-table": "resistance",
    "mge-virulence-table": "virulence",
    "taxonomy-source-table": "taxonomy",
    "taxonomy-summary-table": "taxonomy",
  };
  const tone = tones[containerId] || "neutral";
  [
    "table-tone-neutral",
    "table-tone-taxonomy",
    "table-tone-assembly",
    "table-tone-quality",
    "table-tone-gene",
    "table-tone-risk",
    "table-tone-virulence",
    "table-tone-resistance",
    "table-tone-typing",
  ].forEach((className) => container.classList.remove(className));
  container.classList.add(`table-tone-${tone}`);
}

function renderBinningMetricCards(containerId, items = []) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (!Array.isArray(items) || !items.length) {
    container.innerHTML = "";
    return;
  }
  container.innerHTML = items.map((item) => `
    <article class="result-card">
      <div class="card-head">
        <h3>${escapeHtml(item.title || item.label || "摘要")}</h3>
        ${item.tag ? `<span class="card-tag">${escapeHtml(item.tag)}</span>` : ""}
      </div>
      <div class="qc-summary-feature">
        <strong>${escapeHtml(String(item.value ?? "--"))}</strong>
        <p>${escapeHtml(item.label || "")}</p>
      </div>
      ${item.note ? `<p class="empty-copy">${escapeHtml(item.note)}</p>` : ""}
    </article>
  `).join("");
}

function renderBinningChartCards(containerId, charts = []) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (!Array.isArray(charts) || !charts.length) {
    container.innerHTML = "";
    return;
  }
  container.innerHTML = charts.map((chart) => `
    <article class="result-card">
      ${renderBarSvg(Array.isArray(chart.points) ? chart.points : [], {
        label: chart.label || "统计图",
        width: 980,
        height: 360,
        padX: 54,
        padTop: 24,
        padBottom: 96,
        xLabel: chart.x_label || "分类",
        yLabel: chart.y_label || "数量",
        xValues: Array.isArray(chart.x_values) ? chart.x_values : [],
      })}
    </article>
  `).join("");
}

function buildTableCard(containerId, title, columns, rows) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.dataset.exportTitle = title;
  applyTableTone(container, containerId);
  if (!Array.isArray(rows) || rows.length === 0) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>${escapeHtml(title)} 暂无数据</strong>
        <p class="empty-copy">本节未检出可展示结果。</p>
      </div>
    `;
    return;
  }
  if (['binning-quality-table', 'binning-taxonomy-table', 'assembly-species-table', 'contig-annotation-table', 'rv-summary-table', 'virulence-table', 'resistance-table', 'mge-resistance-table', 'mge-virulence-table', 'patho-mutate-table', 'community-demux-table', 'community-denoise-table', 'community-taxonomy-table', 'community-alpha-table', 'community-alpha-pairwise-table', 'community-differential-table', 'community-rf-table', 'community-beta-pcoa-table', 'community-beta-nmds-table', 'community-beta-stats-table', 'community-beta-distance-table', 'community-network-module-table', 'community-network-node-table', 'community-network-edge-table', 'community-outputs-table', 'resistance-neisseria-amr-table'].includes(containerId)) {
    renderInteractiveContigTable(container, columns, rows, containerId);
    return;
  }
  container.innerHTML = `
    ${renderTableExportToolbar()}
    <div class="table-frame">
      <table class="report-table">
        <thead><tr>${columns.map((column) => `<th>${escapeHtml(humanizeReportColumnLabel(column))}</th>`).join("")}</tr></thead>
        <tbody>
          ${rows.map((row) => `<tr>${columns.map((column, index) => {
            const value = Array.isArray(row) ? row[index] : row[column];
            return `<td>${renderTableCellContent(value, column)}</td>`;
          }).join("")}</tr>`).join("")}
        </tbody>
      </table>
    </div>
  `;
  bindTableExportButtons(container, title, columns, rows.map((row) => (
    Array.isArray(row) ? row : columns.map((column) => row[column] ?? "")
  )));
}

function renderSerotypeSection(section) {
  const container = document.getElementById("serotype-table");
  if (!container) return;
  const mode = section?.mode || "generic";
  if (mode === "sars_cov_2_nextclade") {
    const sectionNode = document.getElementById("section-serotype");
    const headingTitle = sectionNode?.querySelector("h2");
    const headingCopy = sectionNode?.querySelector(".section-heading p:last-child");
    if (headingTitle) headingTitle.textContent = "Nextclade 分型";
    if (headingCopy) headingCopy.textContent = "展示 SARS-CoV-2 的 Nextclade clade、Pango 谱系、质控状态与关键突变摘要。";
    const fastpStatus = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "").trim();
    const coveragePoints = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
      ? currentReportData.sections.assembly.coverage.points
      : [];
    const keepVisible = new Set();
    if (fastpStatus === "ready") {
      keepVisible.add("section-raw-qc");
      keepVisible.add("section-fastp");
    }
    if (coveragePoints.length) {
      keepVisible.add("section-assembly");
      keepVisible.add("section-assembly-summary");
    }
    [
      "section-raw-qc",
      "section-fastp",
      "section-species",
      "section-species-identification",
      "section-assembly-species-identification",
      "section-taxonomy-abundance",
      "section-assembly",
      "section-assembly-summary",
      "section-contig-annotation",
      "section-cgview",
      "section-checkm",
      "section-gene-annotation",
      "section-rv",
      "section-rv-summary",
      "section-virulence",
      "section-resistance",
      "section-resistance-mutation",
      "section-mlst",
      "section-priority-serotype",
      "section-mge",
      "section-mge-resistance",
      "section-mge-virulence",
    ].forEach((id) => {
      const element = document.getElementById(id);
      if (!element) return;
      element.classList.toggle("hidden", !keepVisible.has(id));
    });
    const columns = Array.isArray(section?.columns) ? section.columns : [];
    const rows = Array.isArray(section?.rows) ? section.rows : [];
    const normalizeRow = (row) => {
      if (Array.isArray(row)) {
        return Object.fromEntries(columns.map((column, index) => [column, row[index]]));
      }
      return row && typeof row === "object" ? row : {};
    };
    const getRowValue = (row, key) => {
      if (Array.isArray(row)) {
        const index = columns.indexOf(key);
        return index >= 0 ? row[index] : "";
      }
      return row?.[key];
    };
    const summaryRow = normalizeRow(rows[0] || {});
    const coverageValueRaw = Number(getRowValue(summaryRow, "coverage"));
    const coverageValue = Number.isFinite(coverageValueRaw) ? `${(coverageValueRaw * 100).toFixed(1)}%` : "--";
    const qcStatusRaw = String(getRowValue(summaryRow, "qc.overallStatus") || "").trim();
    const qcStatusLabel = qcStatusRaw ? qcStatusRaw.toUpperCase() : "--";
    const summaryColumns = [
      { key: "index", label: "#" },
      { key: "seqName", label: "Sequence name" },
      { key: "qc.overallStatus", label: "QC", value: `<span class="nextclade-qc-pill is-${escapeHtml(qcStatusRaw || "unknown")}">${escapeHtml(qcStatusLabel)}</span>` },
      { key: "clade", label: "Clade" },
      { key: "Nextclade_pango", label: "Pango lineage" },
      { key: "clade_who", label: "WHO name" },
      { key: "totalSubstitutions", label: "Mut." },
      { key: "totalNonACGTNs", label: "non-ACGTN" },
      { key: "totalMissing", label: "Ns" },
      { key: "coverage", label: "Cov.", value: escapeHtml(coverageValue) },
      { key: "totalDeletions", label: "Gaps" },
      { key: "totalInsertions", label: "Ins." },
      { key: "totalFrameShifts", label: "FS" },
      { key: "qc.overallScore", label: "Score" },
    ];
    const summaryHeaderMarkup = summaryColumns
      .map((item) => `<th>${escapeHtml(item.label)}</th>`)
      .join("");
    const summaryValueMarkup = summaryColumns
      .map((item) => {
        const value = item.value != null ? item.value : renderTableCellContent(getRowValue(summaryRow, item.key));
        return `<td>${value}</td>`;
      })
      .join("");
    const nucleotideMutations = parseNextcladeNucleotideMutations(summaryRow);
    const aminoacidMutations = parseNextcladeAaMutations(summaryRow);
    const knowledgeMap = section?.mutation_knowledge?.aa_matches && typeof section.mutation_knowledge.aa_matches === "object"
      ? section.mutation_knowledge.aa_matches
      : {};
    const geneOptions = [
      { value: "__genome__", label: "Nucleotide sequence" },
      ...SARS_COV_2_GENE_RANGES.map((item) => ({ value: item.gene, label: item.gene })),
    ];
    const frameShifts = parseNextcladeFrameShifts(summaryRow);
    const variantAnnotation = section?.variant_annotation && typeof section.variant_annotation === "object"
      ? section.variant_annotation
      : { status: "empty", columns: [], rows: [] };
    const resistanceAnnotation = section?.resistance_annotation && typeof section.resistance_annotation === "object"
      ? section.resistance_annotation
      : { status: "empty", columns: [], rows: [] };
    const igvView = section?.igv && typeof section.igv === "object"
      ? section.igv
      : { status: "empty" };
    const geneSummaryRows = SARS_COV_2_GENE_RANGES.map((item) => {
      const nucleotideCount = nucleotideMutations.filter((entry) => Array.isArray(entry.genes) && entry.genes.includes(item.gene)).length;
      const aminoacidCount = aminoacidMutations.filter((entry) => entry.gene === item.gene).length;
      const frameShiftCount = frameShifts.filter((entry) => Array.isArray(entry.genes) && entry.genes.includes(item.gene)).length;
      const knowledgeHits = aminoacidMutations.filter((entry) => (
        entry.gene === item.gene
        && entry.label
        && Array.isArray(knowledgeMap[entry.label])
        && knowledgeMap[entry.label].length > 0
      )).length;
      return {
        gene: item.gene,
        nucleotide: nucleotideCount,
        aminoacid: aminoacidCount,
        frameshift: frameShiftCount,
        knowledgeHits,
        total: nucleotideCount + aminoacidCount + frameShiftCount,
      };
    });
    const featurePayload = encodeURIComponent(JSON.stringify({
      nucleotide: nucleotideMutations,
      aminoacid: aminoacidMutations,
      frameShifts,
      genes: geneSummaryRows,
      knowledge: section?.mutation_knowledge || {},
    }));
    applyTableTone(container, "serotype-table");
    const compactFacts = [
      { label: "Clade", value: section?.predicted_clade || getRowValue(summaryRow, "clade") || "--" },
      { label: "Pango", value: section?.pango_lineage || getRowValue(summaryRow, "Nextclade_pango") || "--" },
      { label: "QC", value: qcStatusRaw || "--" },
      { label: "Coverage", value: coverageValue },
    ];
    const compactFactMarkup = compactFacts.map((item) => `
      <span class="nextclade-compact-fact">
        <b>${escapeHtml(item.label)}</b>
        <strong>${escapeHtml(String(item.value ?? "--"))}</strong>
      </span>
    `).join("");
    const variantSummaryCards = [
      { label: "总变异位点", value: String(variantAnnotation?.total_variants ?? "--") },
      { label: "高质量突变", value: String(variantAnnotation?.high_quality_variants ?? "--") },
      { label: "低质量突变", value: String(variantAnnotation?.low_quality_variants ?? "--") },
      { label: "注释VCF", value: String(variantAnnotation?.source_vcf ? "snps.raw.ann.vcf" : "--") },
    ];
    const variantSummaryMarkup = variantSummaryCards.map((item) => `
      <article class="mini-stat-card">
        <span>${escapeHtml(item.label)}</span>
        <strong>${escapeHtml(item.value)}</strong>
      </article>
    `).join("");
    container.dataset.exportTitle = "新冠 Nextclade 分型";
    container.innerHTML = `
      <div class="serotype-special-layout serotype-nextclade-layout">
        <section id="nextclade-summary" class="nextclade-summary-table-card" data-report-nav-anchor>
          <div class="table-frame nextclade-summary-frame">
            <table class="report-table nextclade-summary-table">
              <thead><tr>${summaryHeaderMarkup}</tr></thead>
              <tbody><tr>${summaryValueMarkup}</tr></tbody>
            </table>
          </div>
        </section>
        <section id="nextclade-assignment" class="nextclade-compact-assignment" data-report-nav-anchor>
          <span class="section-chip">Assignment</span>
          <div class="nextclade-compact-facts">${compactFactMarkup}</div>
        </section>
        ${String(variantAnnotation?.status || "") === "ready" && Array.isArray(variantAnnotation?.rows) && variantAnnotation.rows.length ? `
          <section id="nextclade-variant-annotation" class="result-card nextclade-variant-annotation-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">SnpEff</span>
                <h3>变异注释表</h3>
              </div>
              <span class="card-tag">snps.raw.mutation_table</span>
            </div>
            <p class="nextclade-variant-annotation-copy">读取 freebayes 的原始突变位点，基于新冠参考注释执行 snpEff 注释，并按 QUAL、深度和 MAF 划分高低质量。</p>
            <div id="nextclade-variant-annotation-summary" class="mini-stat-grid">${variantSummaryMarkup}</div>
            <div class="nextclade-variant-tabs" id="nextclade-variant-tabs" role="tablist" aria-label="新冠变异质量分层切换">
              <button type="button" class="report-tab-button active" data-nextclade-variant-tab="high">高质量突变</button>
              <button type="button" class="report-tab-button" data-nextclade-variant-tab="low">低质量突变</button>
            </div>
            <div id="nextclade-variant-annotation-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${String(igvView?.status || "") === "ready" ? `
          <section id="nextclade-igv" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">IGV</span>
                <h3>IGV 比对结果</h3>
              </div>
              <span class="card-tag">${escapeHtml(String(igvView?.viewer_label || "参考比对视图"))}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">${escapeHtml(String(igvView?.note || "点击上方变异位点后，IGV 会自动跳转到对应位置。"))}</p>
            <div class="report-igv-shell">
              <div id="nextclade-igv-lazy" class="empty-box">
                <p>IGV 改为按需加载，避免页面初始卡顿。点击下方按钮或上方突变位点后会自动开始加载。</p>
                <button type="button" id="nextclade-igv-load" class="table-export-button">加载 IGV</button>
              </div>
              <iframe id="nextclade-igv-frame" class="report-igv-frame" title="新冠 IGV 比对结果" hidden loading="lazy"></iframe>
            </div>
          </section>
        ` : ""}
        ${String(section?.phylogeny_tree?.status || "") === "ready" ? `
          <section id="nextclade-phylogeny" class="nextclade-phylogeny-section" data-report-nav-anchor>
            <div id="nextclade-phylogeny-tree"></div>
          </section>
        ` : ""}
        <section id="nextclade-gene-workspace" class="result-card nextclade-feature-workspace" data-nextclade-feature-summary="${featurePayload}" data-report-nav-anchor>
          <div class="nextclade-feature-toolbar">
            <h3>按基因查看突变</h3>
            <div class="nextclade-feature-select-wrap">
              <select id="nextclade-feature-select" class="nextclade-feature-select" data-nextclade-feature-select>
                ${geneOptions.map((item) => `<option value="${escapeHtml(item.value)}"${item.value === "S" ? " selected" : ""}>${escapeHtml(item.label)}</option>`).join("")}
              </select>
            </div>
          </div>
          <div data-nextclade-feature-gene-table></div>
          <div data-nextclade-feature-overview></div>
          <div id="nextclade-linked-track" data-nextclade-feature-plot data-report-nav-anchor></div>
          <div id="nextclade-linked-table" class="nextclade-feature-detail-workspace" data-nextclade-detail-workspace data-report-nav-anchor>
            <div class="nextclade-detail-head">
              <div>
                <span class="section-chip">Mutation details</span>
                <h3>突变位点与知识库关联表</h3>
              </div>
              <p>按核苷酸位点映射 CDS 氨基酸位置，并关联本地 NGDC 新冠突变表型知识库注释。</p>
            </div>
            <div data-nextclade-linked-mutation-table></div>
          </div>
        </section>
      </div>
    `;
    if (String(section?.phylogeny_tree?.status || "") === "ready") {
      renderITOLTreeCard("nextclade-phylogeny-tree", {
        ...section.phylogeny_tree,
        task_id: currentReportData?.task?.id || "",
        label: section.phylogeny_tree?.label || "Nextclade 系统发育树",
      });
    }
    if (String(variantAnnotation?.status || "") === "ready" && Array.isArray(variantAnnotation?.rows) && variantAnnotation.rows.length) {
      const annotationTable = document.getElementById("nextclade-variant-annotation-table");
      if (annotationTable) {
        annotationTable.dataset.exportTitle = "新冠_snpEff_变异注释表";
        const columns = Array.isArray(variantAnnotation?.columns) ? variantAnnotation.columns : [];
        const rows = Array.isArray(variantAnnotation?.rows) ? variantAnnotation.rows : [];
        const qualityIndex = columns.indexOf("质量分层");
        const highRows = qualityIndex >= 0 ? rows.filter((row) => String(row[qualityIndex] ?? "").trim() === "高质量突变") : rows;
        const lowRows = qualityIndex >= 0 ? rows.filter((row) => String(row[qualityIndex] ?? "").trim() === "低质量突变") : [];
        const renderVariantTab = (tabKey = "high") => {
          const activeRows = tabKey === "low" ? lowRows : highRows;
          annotationTable.dataset.exportTitle = tabKey === "low" ? "新冠_snpEff_低质量变异注释表" : "新冠_snpEff_高质量变异注释表";
          renderInteractiveContigTable(
            annotationTable,
            columns,
            activeRows,
            "nextclade-variant-annotation-table",
          );
          document.querySelectorAll("[data-nextclade-variant-tab]").forEach((button) => {
            button.classList.toggle("active", button.getAttribute("data-nextclade-variant-tab") === tabKey);
          });
        };
        renderVariantTab("high");
        document.querySelectorAll("[data-nextclade-variant-tab]").forEach((button) => {
          button.addEventListener("click", () => {
            renderVariantTab(String(button.getAttribute("data-nextclade-variant-tab") || "high"));
          });
        });
      }
    }
    const annotationTable = document.getElementById("nextclade-variant-annotation-table");
    const igvFrame = document.getElementById("nextclade-igv-frame");
    if (igvFrame && annotationTable && String(igvView?.status || "") === "ready") {
      const firstMutation = Array.isArray(variantAnnotation?.rows) && variantAnnotation.rows.length ? variantAnnotation.rows[0] : null;
      const chromIndex = Array.isArray(variantAnnotation?.columns) ? variantAnnotation.columns.indexOf("染色体") : -1;
      const posIndex = Array.isArray(variantAnnotation?.columns) ? variantAnnotation.columns.indexOf("位置") : -1;
      const initialLocus = firstMutation && chromIndex >= 0 && posIndex >= 0
        ? `${String(firstMutation[chromIndex] || "").trim()}:${String(firstMutation[posIndex] || "").trim()}`
        : "";
      initializeDeferredIgvEmbed({
        frameId: "nextclade-igv-frame",
        panelId: "nextclade-igv-lazy",
        buttonId: "nextclade-igv-load",
        task: currentReportData?.task || {},
        igvView,
        initialLocus,
        mutationTableNode: annotationTable,
      });
    }
    return;
  }
  if (mode === "rsv_nextclade" || mode === "hmpv_nextclade" || mode === "denv_nextclade" || mode === "zikav_nextclade" || mode === "chikv_nextclade" || mode === "ebola_nextclade" || mode === "hpiv_typing" || mode === "hadv_typing" || mode === "norovirus_typing" || mode === "enterovirus_typing" || mode === "hepatovirus_typing" || mode === "bandavirus_typing" || mode === "orthohantavirus_typing" || mode === "orthoebolavirus_typing" || mode === "astroviridae_typing" || mode === "rhinovirus_typing" || mode === "seasonal_hcov_typing" || mode === "rotavirus_typing" || mode === "hiv_resistance") {
    const isHadv = mode === "hadv_typing";
    const isHpiv = mode === "hpiv_typing";
    const isHmpv = mode === "hmpv_nextclade";
    const isDenv = mode === "denv_nextclade";
    const isZikav = mode === "zikav_nextclade";
    const isChikv = mode === "chikv_nextclade";
    const isEbola = mode === "ebola_nextclade";
    const isNorovirus = mode === "norovirus_typing";
    const isEnterovirus = mode === "enterovirus_typing";
    const isHiv = mode === "hiv_resistance";
    const isHepatovirus = mode === "hepatovirus_typing";
    const isBandavirus = mode === "bandavirus_typing";
    const isOrthohantavirus = mode === "orthohantavirus_typing";
    const isAstroviridae = mode === "astroviridae_typing";
    const isRhinovirus = mode === "rhinovirus_typing";
    const isSeasonalHcov = mode === "seasonal_hcov_typing";
    const isRotavirus = mode === "rotavirus_typing";
    const hepatovirusBroad = String(section?.predicted_group || "").trim().toUpperCase();
    const hepatovirusLabel = ({ HAV: "甲型肝炎病毒", HBV: "乙型肝炎病毒", HCV: "丙型肝炎病毒", HDV: "丁型肝炎病毒", HEV: "戊型肝炎病毒" })[hepatovirusBroad] || "肝炎病毒";
    const virusShort = isHiv ? "HIV" : (isRotavirus ? "RotaV" : (isNorovirus ? "NoV" : (isEnterovirus ? "EV" : (isHepatovirus ? (hepatovirusBroad || "HepV") : (isBandavirus ? "BandV" : (isOrthohantavirus ? "HTNV" : (isEbola ? "EBOV" : (isAstroviridae ? "AstV" : (isRhinovirus ? "HRV" : (isSeasonalHcov ? "HCoV" : (isChikv ? "CHIKV" : (isZikav ? "ZIKV" : (isDenv ? "DENV" : (isHmpv ? "HMPV" : (isHpiv ? "HPIV" : (isHadv ? "HAdV" : "RSV"))))))))))))))));
    const virusLabel = isHiv ? "HIV" : (isRotavirus ? "轮状病毒" : (isNorovirus ? "诺如病毒" : (isEnterovirus ? "肠道病毒" : (isHepatovirus ? hepatovirusLabel : (isBandavirus ? "班达病毒" : (isOrthohantavirus ? "汉坦病毒" : (isEbola ? "埃博拉病毒" : (isAstroviridae ? "星状病毒" : (isRhinovirus ? "鼻病毒" : (isSeasonalHcov ? "季节性冠状病毒" : (isChikv ? "基孔肯雅病毒" : (isZikav ? "寨卡病毒" : (isDenv ? "登革热病毒" : (isHmpv ? "人偏肺病毒" : (isHpiv ? "人副流感病毒" : (isHadv ? "人腺病毒" : "RSV"))))))))))))))));
    const typingCopy = isNorovirus
      ? "基于 CDC RdRp 与 VP1 双位点参考库进行诺如病毒分型，再从对应候选全基因组参考中选择覆盖度最优者完成组装展示、变异注释与结果判读。"
      : isEnterovirus
      ? "基于 EV-A/B/C/D 的 VP1 参考库完成肠道病毒分型，在同亚型完整基因组候选集中经 95% 去冗余和分批竞争选择最优参考，并结合 VADR 注释展示组装与变异结果。"
      : isHiv
      ? "先基于 HIV-1/HIV-2 broad 参考库完成大亚型区分；若命中 HIV-1，再按子亚型代表株覆盖度选择一致性生成参考，随后整合 REGA-like 子亚型/重组证据与 HIVDB 耐药解释。"
      : isHepatovirus
      ? "先基于肝炎病毒 broad 参考库完成大亚型判定；再进入对应大亚型参考库选择覆盖度最优的子亚型/基因型参考，并据此展示组装、突变与比对结果。"
      : isBandavirus
      ? "先用 Bandavirus reference_genomes 判定大亚型；若为 SFTSV，再结合 A_F 与 CJ 的 L/M/S 三片段分型结果，选择最优参考并汇总重组提示、分段证据与 IGV 结果。"
      : isOrthohantavirus
      ? "先用 Orthohantavirus 本地三片段参考库完成 broad 型别筛选，再查看 L/M/S 三片段最优参考与 consensus 复核结果，并结合 IGV 验证当前样本的分段支持证据。"
      : isEbola
      ? "先用 Orthoebolavirus 本地参考库完成最优参考筛选，再对 EBOV 样本调用 Ebola Nextclade 数据集，展示 clade、lineage、质控指标、突变位点、系统发育树与参考比对结果。"
      : isAstroviridae
      ? "基于 ORF2 参考库完成星状病毒分型，在同亚型完整基因组候选集中经 95% 去冗余和分批竞争选择最优参考，并结合 VADR 注释与 ORF2 系统树完成二次分类。"
      : isRhinovirus
      ? "基于 VP1 参考库进行鼻病毒分型，确定物种组后从对应全基因组候选集中选择覆盖度最优参考，并结合 VADR 注释完成组装展示、变异注释与结果判读。"
      : isSeasonalHcov
      ? "先对 229E、NL63、OC43、HKU1 四类季节性冠状病毒进行参考筛选，再基于 consensus 序列执行 VADR 注释、S 基因提取和系统发育树比对完成大类分型与子亚型判读。"
      : isRotavirus
      ? "先基于 A/B/C 三个大组完整参考的全长覆盖度完成大组判定，再对 A 组样本结合 VP4/VP7 组合分型选择最优参考株，并汇总 group coverage 与 subtype typing 结果。"
      : isHpiv
      ? "根据 HPIV1/2/3/4A/4B 参考序列覆盖度自动选择最优参考，并使用对应 GFF 注释文件完成有参组装、变异注释与结果展示。"
      : (isHadv
      ? "先对 Penton、Hexon、Fiber 三个基因片段进行分型，再根据 PHF 组合确定 HAdV 总分型，并从对应分型全基因组候选库中选择覆盖度最优参考完成结果展示。"
      : (isChikv
      ? "基于 CHIKV 固定参考基因组与对应 Nextclade 数据集，展示 clade、lineage、质控指标、突变位点、系统发育树与参考比对结果。"
      : (isZikav
      ? "基于固定 ZIKV 参考基因组与对应 Nextclade 数据集，展示 clade、lineage 与参考比对结果。"
      : (isDenv
      ? "根据自动选择的 DENV1-4 参考株传入对应注释文件，并使用对应亚型的 Nextclade 数据集展示 clade、lineage 与参考比对结果。"
      : (isHmpv
        ? "基于固定 HMPV 参考基因组与对应 Nextclade 数据集，展示 clade、lineage 与参考比对结果。"
        : "根据自动选择的 RSV A/B 参考株分别调用对应 Nextclade 数据集，并展示 clade、lineage 与参考比对结果。")))));
    const igvNote = isNorovirus
      ? "展示自动选择的 Norovirus 最优参考株比对结果。"
      : isEnterovirus
      ? "展示自动选择的 Enterovirus 最优参考株比对结果。"
      : isHiv
      ? "展示当前 HIV 样本与自动选择代表株参考的比对结果，可结合子亚型、重组和耐药结果一起判读。"
      : isHepatovirus
      ? "展示自动选择的肝炎病毒最优参考株比对结果；可结合大亚型、子亚型/基因型与知识库参考关联一起判读。"
      : isAstroviridae
      ? "展示自动选择的 Astrovirus 最优参考株比对结果。"
      : isRhinovirus
      ? "展示自动选择的 Rhinovirus 最优参考株比对结果。"
      : isOrthohantavirus
      ? "展示自动选择的 Orthohantavirus 最优参考株比对结果，可与 broad 筛选和 L/M/S 三片段证据联合判读。"
      : isEbola
      ? "展示自动选择的 Orthoebolavirus 参考序列比对结果，可与本地参考筛选、Ebola Nextclade clade / lineage 和突变位点联合判读。"
      : isSeasonalHcov
      ? "展示自动选择的季节性冠状病毒最优参考株比对结果。"
      : isRotavirus
      ? "当前 rotavirus demo 主要展示参考选择与分型结果，尚未生成可嵌入的 IGV 会话。"
      : isHpiv
      ? "展示自动选择的 HPIV 最优参考株比对结果。"
      : (isHadv
      ? "展示自动选择的 HAdV 最优参考株比对结果。"
      : (isChikv
      ? "展示 CHIKV 固定参考株比对结果，可与分型总表、突变位点表和系统发育树联合判读。"
      : (isZikav
      ? "展示 ZIKV 参考株比对结果。"
      : (isDenv
      ? "展示自动选择的 DENV 参考株比对结果。"
      : (isHmpv
        ? "展示 HMPV 参考株比对结果。"
        : "展示自动选择的 RSV 参考株比对结果。")))));
    const sectionNode = document.getElementById("section-serotype");
    const headingTitle = sectionNode?.querySelector("h2");
    const headingCopy = sectionNode?.querySelector(".section-heading p:last-child");
    const hasFastp = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "") === "ready";
    const hasSpecies = Array.isArray(currentReportData?.sections?.species_identification?.species?.rows)
      && currentReportData.sections.species_identification.species.rows.length > 0;
    const hasSubspecies = Array.isArray(currentReportData?.sections?.species_identification?.subspecies?.rows)
      && currentReportData.sections.species_identification.subspecies.rows.length > 0;
    const hasCoverage = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
      && currentReportData.sections.assembly.coverage.points.length > 0;
    const hasAssemblySummary = Array.isArray(currentReportData?.sections?.assembly?.summary?.rows)
      && currentReportData.sections.assembly.summary.rows.length > 0;
    const keepVisible = new Set(["section-raw-qc", "section-serotype"]);
    if (hasFastp) keepVisible.add("section-fastp");
    if (hasSpecies || hasSubspecies) keepVisible.add("section-species-identification");
    if (hasCoverage) keepVisible.add("section-assembly");
    if (hasAssemblySummary) keepVisible.add("section-assembly-summary");
    [
      "section-raw-qc",
      "section-fastp",
      "section-species",
      "section-species-identification",
      "section-taxonomy-abundance",
      "section-assembly",
      "section-assembly-summary",
      "section-contig-annotation",
      "section-cgview",
      "section-checkm",
      "section-gene-annotation",
      "section-rv",
      "section-rv-summary",
      "section-virulence",
      "section-resistance",
      "section-resistance-mutation",
      "section-mlst",
      "section-priority-serotype",
      "section-mge",
      "section-mge-resistance",
      "section-mge-virulence",
    ].forEach((id) => {
      const element = document.getElementById(id);
      if (!element) return;
      element.classList.toggle("hidden", !keepVisible.has(id));
    });
    if (headingTitle) headingTitle.textContent = `${virusLabel} 分型`;
    if (headingCopy) headingCopy.textContent = typingCopy;
    applyTableTone(container, "serotype-table");
    const summaryCards = Array.isArray(section?.summary_cards) ? section.summary_cards : [];
    const qualityMetrics = Array.isArray(section?.quality_metrics) ? section.quality_metrics : [];
    const notes = String(section?.notes || "").trim();
    const columns = Array.isArray(section?.columns) ? section.columns : [];
    const rows = Array.isArray(section?.rows) ? section.rows : [];
    const mutationTable = section?.mutation_table && typeof section.mutation_table === "object"
      ? section.mutation_table
      : { status: "empty", columns: [], rows: [] };
    const nmdcAnnotation = !isHmpv && !isZikav && !isChikv && !isEbola && !isNorovirus && !isEnterovirus && !isHepatovirus && !isBandavirus && !isRhinovirus && !isSeasonalHcov && !isRotavirus && section?.nmdc_annotation && typeof section.nmdc_annotation === "object"
      ? section.nmdc_annotation
      : { status: "empty", columns: [], rows: [] };
    const functionalAnnotation = isHpiv && section?.functional_annotation && typeof section.functional_annotation === "object"
      ? section.functional_annotation
      : { status: "empty", columns: [], rows: [] };
    const phfTable = isHadv && section?.phf_table && typeof section.phf_table === "object"
      ? section.phf_table
      : { status: "empty", columns: [], rows: [] };
    const phfCoverage = isHadv && section?.phf_coverage && typeof section.phf_coverage === "object"
      ? section.phf_coverage
      : { status: "empty", columns: [], rows: [] };
    const phfSnp = isHadv && section?.phf_snp && typeof section.phf_snp === "object"
      ? section.phf_snp
      : { status: "empty", columns: [], rows: [] };
    const igvView = section?.igv && typeof section.igv === "object"
      ? section.igv
      : { status: "empty" };
    const groupTyping = isRotavirus && section?.group_typing && typeof section.group_typing === "object"
      ? section.group_typing
      : { status: "empty", columns: [], rows: [] };
    const subtypeTyping = isRotavirus && section?.subtype_typing && typeof section.subtype_typing === "object"
      ? section.subtype_typing
      : { status: "empty", columns: [], rows: [] };
    const bandavirusSelection = isBandavirus && section?.bandavirus_selection && typeof section.bandavirus_selection === "object"
      ? section.bandavirus_selection
      : { status: "empty", columns: [], rows: [] };
    const orthohantavirusSelection = isOrthohantavirus && section?.orthohantavirus_selection && typeof section.orthohantavirus_selection === "object"
      ? section.orthohantavirus_selection
      : { status: "empty", columns: [], rows: [] };
    const orthoebolavirusSelection = isEbola && section?.orthoebolavirus_selection && typeof section.orthoebolavirus_selection === "object"
      ? section.orthoebolavirus_selection
      : { status: "empty", columns: [], rows: [] };
    const bandavirusAfSegments = isBandavirus && section?.af_segment_typing && typeof section.af_segment_typing === "object"
      ? section.af_segment_typing
      : { status: "empty", columns: [], rows: [] };
    const bandavirusCjSegments = isBandavirus && section?.cj_segment_typing && typeof section.cj_segment_typing === "object"
      ? section.cj_segment_typing
      : { status: "empty", columns: [], rows: [] };
    const orthohantavirusBroadTyping = isOrthohantavirus && section?.broad_typing && typeof section.broad_typing === "object"
      ? section.broad_typing
      : { status: "empty", columns: [], rows: [] };
    const orthohantavirusSegments = isOrthohantavirus && section?.segment_typing && typeof section.segment_typing === "object"
      ? section.segment_typing
      : { status: "empty", columns: [], rows: [] };
    const consensusTyping = (isBandavirus || isOrthohantavirus) && section?.consensus_typing && typeof section.consensus_typing === "object"
      ? section.consensus_typing
      : { status: "empty", columns: [], rows: [] };
    const typingGenePhylogeny = (isNorovirus || isEnterovirus || isAstroviridae || isRhinovirus || isSeasonalHcov) && section?.gene_phylogeny && typeof section.gene_phylogeny === "object"
      ? section.gene_phylogeny
      : { status: "empty", trees: [] };
    const hivResistanceTable = isHiv && section?.resistance_table && typeof section.resistance_table === "object"
      ? section.resistance_table
      : { status: "empty", columns: [], rows: [] };
    const hivResistancePayload = isHiv && section?.resistance_payload && typeof section.resistance_payload === "object"
      ? section.resistance_payload
      : {};
    const hivReferenceSelection = isHiv && section?.reference_selection && typeof section.reference_selection === "object"
      ? section.reference_selection
      : { status: "empty", columns: [], rows: [] };
    const hivBroadTyping = isHiv && section?.broad_typing && typeof section.broad_typing === "object"
      ? section.broad_typing
      : { status: "empty", columns: [], rows: [] };
    const hivSubtypeReferenceTyping = isHiv && section?.subtype_reference_typing && typeof section.subtype_reference_typing === "object"
      ? section.subtype_reference_typing
      : { status: "empty", columns: [], rows: [] };
    const hivMutationPanels = isHiv && Array.isArray(section?.mutation_panels)
      ? section.mutation_panels
      : [];
    const hivSubtypingSummary = isHiv && section?.subtyping_summary && typeof section.subtyping_summary === "object"
      ? section.subtyping_summary
      : {};
    const hivBootscanAssets = isHiv && section?.bootscan_assets && typeof section.bootscan_assets === "object"
      ? section.bootscan_assets
      : {};
    const hivBootscanEmbedded = isHiv && section?.bootscan_embedded && typeof section.bootscan_embedded === "object"
      ? section.bootscan_embedded
      : {};
    const typingPhylogenyTrees = Array.isArray(typingGenePhylogeny?.trees)
      ? typingGenePhylogeny.trees.filter((item) => String(item?.status || "") === "ready")
      : [];
    const mutationSummaryCards = [
      { label: "总变异位点", value: String(mutationTable?.total_variants ?? (Array.isArray(mutationTable?.rows) ? mutationTable.rows.length : "--")) },
      { label: "高质量突变", value: String(mutationTable?.high_quality_variants ?? "--") },
      { label: "低质量突变", value: String(mutationTable?.low_quality_variants ?? "--") },
      { label: "注释VCF", value: String(mutationTable?.source_vcf ? "snps.anno.vcf" : "--") },
    ];
    const mutationSummaryMarkup = mutationSummaryCards.map((item) => `
      <article class="mini-stat-card">
        <span>${escapeHtml(item.label)}</span>
        <strong>${escapeHtml(item.value)}</strong>
      </article>
    `).join("");
    const reportTaskId = String(currentReportData?.task?.id || "").trim();
    const hivOverallSvgSrc = isHiv && hivBootscanAssets?.overall_svg ? buildReportAssetUrl(reportTaskId, hivBootscanAssets.overall_svg) : "";
    const hivPureSvgSrc = isHiv && hivBootscanAssets?.pure_svg ? buildReportAssetUrl(reportTaskId, hivBootscanAssets.pure_svg) : "";
    const hivOverallCsvAsset = isHiv ? String(hivBootscanAssets?.overall_csv || "").trim() : "";
    const hivPureCsvAsset = isHiv ? String(hivBootscanAssets?.pure_csv || "").trim() : "";
    const summaryMarkup = summaryCards.length ? `
      <div class="mini-stat-grid">
        ${summaryCards.map((item) => `
          <article class="mini-stat-card">
            <span>${escapeHtml(item?.label || "--")}</span>
            <strong>${escapeHtml(String(item?.value ?? "--"))}</strong>
          </article>
        `).join("")}
      </div>
    ` : "";
    const qualityMarkup = qualityMetrics.length ? `
      <div class="mini-stat-grid">
        ${qualityMetrics.map((item) => `
          <article class="mini-stat-card">
            <span>${escapeHtml(item?.label || "--")}</span>
            <strong>${escapeHtml(String(item?.value ?? "--"))}</strong>
          </article>
        `).join("")}
      </div>
    ` : "";
    container.dataset.exportTitle = `${virusShort} Nextclade 分型`;
    container.innerHTML = `
      <div class="serotype-special-layout influenza-typing-layout">
        ${summaryMarkup}
        ${notes ? `
          <div class="chart-insight serotype-insight" role="note" aria-label="${virusLabel} 分型说明">
            <span class="chart-insight-label">结果说明</span>
            <p>${escapeHtml(notes)}</p>
          </div>
        ` : ""}
        <section id="rsv-typing-summary" class="result-card" data-report-nav-anchor>
          <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">${(isHpiv || isHadv || isNorovirus || isBandavirus || isSeasonalHcov || isRotavirus) ? "Reference typing" : "Nextclade summary"}</span>
              <h3>${virusLabel} 分型总表</h3>
            </div>
          </div>
          <div id="rsv-typing-summary-table" class="report-table-card report-table-card-embedded"></div>
        </section>
        ${isHadv && Array.isArray(phfTable?.rows) && phfTable.rows.length ? `
          <section id="rsv-typing-phf" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">PHF typing</span>
                <h3>PHF 三基因分型表</h3>
              </div>
            </div>
            <div id="rsv-typing-phf-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isHadv && Array.isArray(phfCoverage?.rows) && phfCoverage.rows.length ? `
          <section id="rsv-typing-phf-coverage" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">PHF coverage</span>
                <h3>PHF 三基因覆盖度图</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">基于 PHF 三基因分型阶段的 <code>samtools coverage</code> 汇总结果，展示 Penton、Hexon、Fiber 最终命中参考的覆盖度和平均深度，便于快速判断三个分型基因的命中稳定性。</p>
            <div id="rsv-typing-phf-coverage-charts"></div>
            <div id="rsv-typing-phf-coverage-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isHadv && Array.isArray(phfSnp?.rows) && phfSnp.rows.length ? `
          <section id="rsv-typing-phf-snp" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">PHF SNP</span>
                <h3>PHF 三基因差异 SNP 数</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">以当前样本 consensus 序列为基础，分别和 PHF 分型阶段命中的 Penton、Hexon、Fiber 参考基因序列进行比对，统计三者各自的差异 SNP 数。</p>
            <div id="rsv-typing-phf-snp-chart"></div>
            <div id="rsv-typing-phf-snp-table" class="report-table-card report-table-card-embedded"></div>
            ${Array.isArray(phfSnp?.detail_rows) && phfSnp.detail_rows.length ? `
              <div id="rsv-typing-phf-snp-detail-table" class="report-table-card report-table-card-embedded"></div>
            ` : ""}
          </section>
        ` : ""}
        ${qualityMarkup ? `
          <section id="rsv-typing-quality" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">${(isHpiv || isHadv || isNorovirus || isHiv || isBandavirus || isOrthohantavirus || isRhinovirus || isSeasonalHcov || isRotavirus) ? "Typing metrics" : "QC metrics"}</span>
                <h3>${(isHpiv || isHadv || isNorovirus || isHiv || isBandavirus || isOrthohantavirus || isRhinovirus || isSeasonalHcov || isRotavirus) ? "分型指标" : "Nextclade 质量指标"}</h3>
              </div>
            </div>
            ${qualityMarkup}
          </section>
        ` : ""}
        ${isHiv && Array.isArray(hivBroadTyping?.rows) && hivBroadTyping.rows.length ? `
          <section id="rsv-typing-hiv-broad" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Broad typing</span>
                <h3>HIV-1 / HIV-2 大亚型筛选</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">先使用 HIV-1 / HIV-2 代表全基因组参考库按覆盖度和支持 reads 进行 broad typing，确定是否继续进入 HIV-1 子亚型和耐药分析链路。</p>
            <div id="rsv-typing-hiv-broad-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isHiv && Array.isArray(hivSubtypeReferenceTyping?.rows) && hivSubtypeReferenceTyping.rows.length ? `
          <section id="rsv-typing-hiv-reference" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">代表参考</span>
                <h3>子亚型代表株覆盖度筛选</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">若 broad typing 命中 HIV-1，则在各子亚型代表株集合中再次比较覆盖度，选择用于一致性生成的代表株参考。</p>
            <div id="rsv-typing-hiv-reference-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isHiv && (Array.isArray(hivReferenceSelection?.rows) && hivReferenceSelection.rows.length) ? `
          <section id="rsv-typing-hiv-selection" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">选择总表</span>
                <h3>HIV 参考选择总表</h3>
              </div>
            </div>
            <div id="rsv-typing-hiv-selection-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isHiv ? `
          <section id="rsv-typing-hiv-subtyping" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">REGA 风格分型</span>
                <h3>子亚型 / 重组摘要</h3>
              </div>
            </div>
            <div class="mini-stat-grid">
              <article class="mini-stat-card">
                <span>分型结论</span>
                <strong>${escapeHtml(localizeHivAssignmentLabel(String(hivSubtypingSummary?.assignment_label || "--")))}</strong>
              </article>
              <article class="mini-stat-card">
                <span>预测大亚型</span>
                <strong>${escapeHtml(String(hivSubtypingSummary?.predicted_group || "--"))}</strong>
              </article>
              <article class="mini-stat-card">
                <span>预测子亚型</span>
                <strong>${escapeHtml(String(hivSubtypingSummary?.predicted_clade || "--"))}</strong>
              </article>
              <article class="mini-stat-card">
                <span>纯亚型树支持度</span>
                <strong>${escapeHtml(String(hivSubtypingSummary?.pure_tree_support || "--"))}</strong>
              </article>
              <article class="mini-stat-card">
                <span>全参考树支持度</span>
                <strong>${escapeHtml(String(hivSubtypingSummary?.overall_tree_support || "--"))}</strong>
              </article>
              <article class="mini-stat-card">
                <span>候选父本数</span>
                <strong>${escapeHtml(String(hivSubtypingSummary?.candidate_parent_count || "--"))}</strong>
              </article>
            </div>
            ${hivOverallCsvAsset || hivPureCsvAsset ? `
              <article class="result-card" style="margin-top:16px;">
                <div class="card-head">
                  <div class="card-title-stack">
                    <span class="section-chip">Bootscan</span>
                    <h3>Bootscan 交互图</h3>
                  </div>
                  <span class="card-tag">交互图</span>
                </div>
                <p class="nextclade-variant-annotation-copy">改用结果页内交互式重绘，默认聚焦主信号窗口；可切换整体 / 纯亚型视角，并按亚型逐条开关曲线，避免原始 SVG 在单屏里把所有信息挤成一团。</p>
                <div id="rsv-typing-hiv-bootscan"></div>
              </article>
            ` : ""}
          </section>
        ` : ""}
        ${isHiv && Array.isArray(hivResistanceTable?.rows) && hivResistanceTable.rows.length ? `
          <section id="rsv-typing-hiv-resistance" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">HIVDB resistance</span>
                <h3>药物耐药解释</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">参考 HIVDB / VISTAS 的阅读节奏重排为“等级概览 -> 分药物解释 -> 基因突变说明”。原始 TSV 仍然保留在折叠明细里，但不再让它占据第一屏。</p>
            <div id="rsv-typing-hiv-resistance-workspace"></div>
            <details class="report-detail-block" style="margin-top:18px;">
              <summary>查看原始药物总表</summary>
              <div id="rsv-typing-hiv-resistance-table" class="report-table-card report-table-card-embedded" style="margin-top:12px;"></div>
            </details>
          </section>
        ` : ""}
        ${isHiv && hivMutationPanels.length ? `
          <section id="rsv-typing-hiv-mutations" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Mutation evidence</span>
                <h3>PR / RT / IN 位点图</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">主视图改成基因坐标上的突变分布图：所有命中位置都会落在对应基因条带上，关键位点再单独挂标签。原始长文本清单保留在折叠明细里，避免第一屏再次变成字墙。</p>
            <div id="rsv-typing-hiv-mutation-map"></div>
          </section>
        ` : ""}
        ${isBandavirus && Array.isArray(bandavirusSelection?.rows) && bandavirusSelection.rows.length ? `
          <section id="rsv-typing-bandavirus-selection" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Reference typing</span>
                <h3>Bandavirus 参考筛选结果</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">展示 Bandavirus 大亚型筛选、SFTSV 的 A_F 子亚型判定，以及疑似重组/重配提示。</p>
            <div id="rsv-typing-bandavirus-selection-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isOrthohantavirus && Array.isArray(orthohantavirusSelection?.rows) && orthohantavirusSelection.rows.length ? `
          <section id="rsv-typing-orthohantavirus-selection" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Reference typing</span>
                <h3>Orthohantavirus 参考筛选结果</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">展示当前样本的 broad 型别筛选结果、最终选中的型别标签，以及是否进入后续 L/M/S 参考组合。</p>
            <div id="rsv-typing-orthohantavirus-selection-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isEbola && Array.isArray(orthoebolavirusSelection?.rows) && orthoebolavirusSelection.rows.length ? `
          <section id="rsv-typing-ebola-selection" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Reference typing</span>
                <h3>Orthoebolavirus 参考筛选结果</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">展示本地 Orthoebolavirus 参考库中覆盖度最高的候选参考，解释当前样本为何进入 EBOV / Ebola Nextclade 判读链路。</p>
            <div id="rsv-typing-ebola-selection-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isOrthohantavirus && Array.isArray(orthohantavirusBroadTyping?.rows) && orthohantavirusBroadTyping.rows.length ? `
          <section id="rsv-typing-orthohantavirus-broad" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Broad typing</span>
                <h3>Broad 分型支持摘要</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">汇总 broad 分型阶段的命中片段数、coverage_sum、depth_sum 和支持 reads，用于判断当前样本是否具备继续进行汉坦 L/M/S 细分的证据。</p>
            <div id="rsv-typing-orthohantavirus-broad-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isOrthohantavirus && Array.isArray(orthohantavirusSegments?.rows) && orthohantavirusSegments.rows.length ? `
          <section id="rsv-typing-orthohantavirus-segments" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">L/M/S evidence</span>
                <h3>L/M/S 三片段最优参考</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">按 L、M、S 三个片段分别展示命中的最优参考 accession、覆盖度、平均深度和支持 reads，作为汉坦病毒分段证据的主体。</p>
            <div id="rsv-typing-orthohantavirus-segments-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isBandavirus && Array.isArray(bandavirusAfSegments?.rows) && bandavirusAfSegments.rows.length ? `
          <section id="rsv-typing-bandavirus-af" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">A_F typing</span>
                <h3>A_F 三片段最优参考</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">按 L/M/S 三个片段分别展示 A_F 分型阶段命中的最优参考 accession、覆盖度与支持深度。</p>
            <div id="rsv-typing-bandavirus-af-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isBandavirus && Array.isArray(bandavirusCjSegments?.rows) && bandavirusCjSegments.rows.length ? `
          <section id="rsv-typing-bandavirus-cj" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">CJ typing</span>
                <h3>CJ 三片段分型结果</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">按 L/M/S 三个片段分别展示 CJ 分型命中的最优参考与 genotype，用于补充 SFTSV 子亚型解释。</p>
            <div id="rsv-typing-bandavirus-cj-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isBandavirus && Array.isArray(consensusTyping?.rows) && consensusTyping.rows.length ? `
          <section id="rsv-typing-bandavirus-consensus" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Consensus review</span>
                <h3>Consensus 分型复核</h3>
              </div>
            </div>
            <div id="rsv-typing-bandavirus-consensus-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isOrthohantavirus && Array.isArray(consensusTyping?.rows) && consensusTyping.rows.length ? `
          <section id="rsv-typing-orthohantavirus-consensus" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Consensus review</span>
                <h3>Consensus 分型复核</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">利用当前样本的 consensus 结果回看 broad 与 L/M/S 选择是否一致，便于确认最终汉坦分型的稳定性。</p>
            <div id="rsv-typing-orthohantavirus-consensus-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isRotavirus && Array.isArray(groupTyping?.rows) && groupTyping.rows.length ? `
          <section id="rsv-typing-group" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Group coverage</span>
                <h3>A/B/C 大组覆盖度比较</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">按完整参考全长覆盖率、命中片段数与支持 reads 展示当前样本对各组候选参考株的匹配情况，用于解释为何最终选择当前大组与参考株。</p>
            <div id="rsv-typing-group-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${isRotavirus && Array.isArray(subtypeTyping?.rows) && subtypeTyping.rows.length ? `
          <section id="rsv-typing-subtype" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">VP4/VP7 typing</span>
                <h3>G/P 组合分型结果</h3>
              </div>
            </div>
            <p class="nextclade-variant-annotation-copy">对 A 组样本进一步结合 VP7 的 G 分型和 VP4 的 P 分型生成组合分型结果，并回填最优参考株。</p>
            <div id="rsv-typing-subtype-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${Array.isArray(mutationTable?.rows) && mutationTable.rows.length ? `
          <section id="rsv-typing-mutations" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Variants</span>
                <h3>突变位点表</h3>
              </div>
              <span class="card-tag">${escapeHtml(`${mutationTable.rows.length} 条`)}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">${isOrthohantavirus ? "读取 <code>snps.filt1.vcf</code> 与基于对应参考 GFF 生成的 <code>snps.anno.vcf</code>，按流程过滤结果区分高低质量并展示汉坦病毒的位点注释。" : "读取 <code>snps.raw.vcf</code>、<code>snps.filt1.vcf</code> 与 <code>snps.anno.vcf</code>，按流程过滤结果区分高低质量并展示 " + virusLabel + " 的位点注释。"} </p>
            <div id="rsv-typing-mutation-summary" class="mini-stat-grid">${mutationSummaryMarkup}</div>
            <div class="nextclade-variant-tabs" id="rsv-variant-tabs" role="tablist" aria-label="${virusLabel} 变异质量分层切换">
              <button type="button" class="report-tab-button active" data-rsv-variant-tab="high">高质量突变</button>
              <button type="button" class="report-tab-button" data-rsv-variant-tab="low">低质量突变</button>
            </div>
            <div id="rsv-typing-mutation-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${String(igvView?.status || "") === "ready" ? `
          <section id="rsv-typing-igv" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">IGV</span>
                <h3>IGV 比对结果</h3>
              </div>
              <span class="card-tag">${escapeHtml(String(igvView?.viewer_label || "参考比对视图"))}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">${escapeHtml(String(igvView?.note || igvNote))}</p>
            <div class="report-igv-shell">
              <div id="rsv-igv-lazy" class="empty-box">
                <p>IGV 改为按需加载，避免页面初始卡顿。点击下方按钮或上方突变位点后会自动开始加载。</p>
                <button type="button" id="rsv-igv-load" class="table-export-button">加载 IGV</button>
              </div>
              <iframe id="rsv-igv-frame" class="report-igv-frame" title="${virusLabel} IGV 比对结果" hidden loading="lazy"></iframe>
            </div>
          </section>
        ` : ""}
        ${isHpiv && String(functionalAnnotation?.status || "") === "ready" ? `
          <section id="rsv-typing-functional-impact" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Functional impact</span>
                <h3>突变功能影响</h3>
              </div>
              <span class="card-tag">${escapeHtml(`${Array.isArray(functionalAnnotation?.rows) ? functionalAnnotation.rows.length : 0} 条`)}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">已将当前 HPIV 样本的突变位点与文献整理的功能位点库进行关联，便于在 IGV 核查后继续查看这些位点已知的功能影响、表型变化和证据来源。</p>
            <div id="rsv-typing-functional-impact-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${!isHmpv && !isZikav && !isChikv && !isEbola && !isNorovirus && !isRhinovirus && !isSeasonalHcov && String(nmdcAnnotation?.status || "") === "ready" ? `
          <section id="rsv-typing-nmdc" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">NMDC annotation</span>
                <h3>RSV 变异风险数据库注释</h3>
              </div>
              <span class="card-tag">${escapeHtml(`${Array.isArray(nmdcAnnotation?.rows) ? nmdcAnnotation.rows.length : 0} 条`)}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">已根据亚型、基因和突变位点，将当前样本的突变结果与 NMDC RSV 变异风险数据库进行关联，便于在 IGV 核查后继续查看已知位点的抗体亲和和氨基酸替换风险信息。</p>
            <div id="rsv-typing-nmdc-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${(((isNorovirus || isEnterovirus || isAstroviridae || isRhinovirus || isSeasonalHcov) && typingPhylogenyTrees.length > 0) || (!isHpiv && !isNorovirus && !isEnterovirus && !isAstroviridae && !isRhinovirus && !isSeasonalHcov && String(section?.phylogeny_tree?.status || "") === "ready")) ? `
          <section id="rsv-phylogeny" class="nextclade-phylogeny-section" data-report-nav-anchor>
            ${(isNorovirus || isEnterovirus || isAstroviridae || isRhinovirus || isSeasonalHcov) ? `
              <div id="rsv-phylogeny-vp1"></div>
              ${isNorovirus ? `<div id="rsv-phylogeny-rdrp"></div>` : ``}
            ` : `
              <div id="rsv-phylogeny-tree"></div>
            `}
          </section>
        ` : ""}
      </div>
    `;
    const summaryTable = document.getElementById("rsv-typing-summary-table");
    if (summaryTable) {
      let preferredColumns = [];
      if (isHpiv) {
        preferredColumns = [
          "sample",
          "hpiv_type",
          "coverage",
          "mean_depth",
          "covered_bases",
          "num_reads",
        ].filter((column) => columns.includes(column));
      } else if (isNorovirus) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "双位点分型",
          "RdRp分型",
          "VP1分型",
          "参考序列",
          "覆盖度",
          "平均深度",
        ].filter((column) => columns.includes(column));
      } else if (isHiv) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "大亚型",
          "子亚型",
          "重组判定",
          "代表株参考",
          "NRTI最高等级",
          "NNRTI最高等级",
          "PI最高等级",
          "INSTI最高等级",
        ].filter((column) => columns.includes(column));
      } else if (isEnterovirus) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "大亚型",
          "VP1分型",
          "参考序列",
          "覆盖度",
          "平均深度",
        ].filter((column) => columns.includes(column));
      } else if (isBandavirus) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "大亚型",
          "A_F(LMS)",
          "CJ(LMS)",
          "分型结果",
          "参考命中",
        ].filter((column) => columns.includes(column));
      } else if (isRhinovirus) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "VP1分型",
          "物种组",
          "参考序列",
          "覆盖度",
          "平均深度",
        ].filter((column) => columns.includes(column));
      } else if (isSeasonalHcov) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "大类分型",
          "S子亚型",
          "最近参考",
          "注释文件",
          "覆盖度",
          "平均深度",
        ].filter((column) => columns.includes(column));
      } else if (isRotavirus) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "大组分型",
          "G分型",
          "P分型",
          "组合分型",
          "最优参考株",
          "参考片段数",
          "全长覆盖度",
        ].filter((column) => columns.includes(column));
      } else if (isHadv) {
        preferredColumns = [
          "样本名称",
          "病毒类型",
          "HAdV分型",
          "Penton分型",
          "Hexon分型",
          "Fiber分型",
          "参考序列",
          "覆盖度",
          "平均深度",
        ].filter((column) => columns.includes(column));
      } else {
        preferredColumns = [
          "seqName",
          "clade",
          "lineage",
          "genotype",
          "qc.overallStatus",
          "qc.overallScore",
          "coverage",
          "totalSubstitutions",
          "totalAminoacidSubstitutions",
          "totalDeletions",
          "totalInsertions",
        ].filter((column) => columns.includes(column));
      }
      const indexes = preferredColumns.map((column) => columns.indexOf(column)).filter((index) => index >= 0);
      const displayColumns = preferredColumns.length ? preferredColumns : columns;
      const displayRows = preferredColumns.length
        ? rows.map((row) => indexes.map((index) => Array.isArray(row) ? row[index] : ""))
        : rows;
      let summaryExportTitle = `${virusShort}_Nextclade_分型总表`;
      if (isHpiv) summaryExportTitle = "HPIV_参考选择分型表";
      else if (isHadv) summaryExportTitle = "HAdV_参考选择分型表";
      else if (isNorovirus) summaryExportTitle = "Norovirus_双位点分型表";
      else if (isHiv) summaryExportTitle = "HIV_分型与耐药摘要表";
      else if (isEnterovirus) summaryExportTitle = "Enterovirus_VP1_分型表";
      else if (isBandavirus) summaryExportTitle = "Bandavirus_LMS_分型总表";
      else if (isEbola) summaryExportTitle = "Ebola_Nextclade_分型总表";
      else if (isAstroviridae) summaryExportTitle = "Astrovirus_ORF2_分型表";
      else if (isRhinovirus) summaryExportTitle = "Rhinovirus_VP1_分型表";
      else if (isSeasonalHcov) summaryExportTitle = "Seasonal_HCoV_S_分型表";
      else if (isRotavirus) summaryExportTitle = "Rotavirus_分型总表";
      summaryTable.dataset.exportTitle = summaryExportTitle;
      renderInteractiveContigTable(summaryTable, displayColumns, displayRows, "rsv-typing-summary-table");
    }
    if (isHiv && ((hivOverallCsvAsset || hivPureCsvAsset) || hivBootscanEmbedded?.overall_csv_text || hivBootscanEmbedded?.pure_csv_text)) {
      initializeHivBootscanExplorer({
        containerId: "rsv-typing-hiv-bootscan",
        taskId: reportTaskId,
        overallCsvAsset: hivOverallCsvAsset,
        pureCsvAsset: hivPureCsvAsset,
        embeddedData: hivBootscanEmbedded,
        summary: hivSubtypingSummary,
      });
    }
    const hivBroadTableNode = document.getElementById("rsv-typing-hiv-broad-table");
    if (hivBroadTableNode && isHiv && Array.isArray(hivBroadTyping?.rows) && hivBroadTyping.rows.length) {
      hivBroadTableNode.dataset.exportTitle = "HIV_broad_typing";
      renderInteractiveContigTable(
        hivBroadTableNode,
        Array.isArray(hivBroadTyping?.columns) ? hivBroadTyping.columns : [],
        hivBroadTyping.rows,
        "rsv-typing-hiv-broad-table",
      );
    }
    const hivReferenceTableNode = document.getElementById("rsv-typing-hiv-reference-table");
    if (hivReferenceTableNode && isHiv && Array.isArray(hivSubtypeReferenceTyping?.rows) && hivSubtypeReferenceTyping.rows.length) {
      hivReferenceTableNode.dataset.exportTitle = "HIV_subtype_reference_selection";
      renderInteractiveContigTable(
        hivReferenceTableNode,
        Array.isArray(hivSubtypeReferenceTyping?.columns) ? hivSubtypeReferenceTyping.columns : [],
        hivSubtypeReferenceTyping.rows,
        "rsv-typing-hiv-reference-table",
      );
    }
    const hivSelectionTableNode = document.getElementById("rsv-typing-hiv-selection-table");
    if (hivSelectionTableNode && isHiv && Array.isArray(hivReferenceSelection?.rows) && hivReferenceSelection.rows.length) {
      hivSelectionTableNode.dataset.exportTitle = "HIV_reference_selection_summary";
      renderInteractiveContigTable(
        hivSelectionTableNode,
        Array.isArray(hivReferenceSelection?.columns) ? hivReferenceSelection.columns : [],
        hivReferenceSelection.rows,
        "rsv-typing-hiv-selection-table",
      );
    }
    const hivResistanceTableNode = document.getElementById("rsv-typing-hiv-resistance-table");
    if (hivResistanceTableNode && isHiv && Array.isArray(hivResistanceTable?.rows) && hivResistanceTable.rows.length) {
      hivResistanceTableNode.dataset.exportTitle = "HIVDB_药物耐药解释表";
      renderHivResistanceWorkspace({
        workspaceId: "rsv-typing-hiv-resistance-workspace",
        rawTableId: "rsv-typing-hiv-resistance-table",
        table: hivResistanceTable,
        payload: hivResistancePayload,
        mutationPanels: hivMutationPanels,
      });
    }
    const hivMutationMapNode = document.getElementById("rsv-typing-hiv-mutation-map");
    if (hivMutationMapNode && isHiv && Array.isArray(hivResistanceTable?.rows) && hivResistanceTable.rows.length) {
      renderHivMutationMap({
        containerId: "rsv-typing-hiv-mutation-map",
        table: hivResistanceTable,
        payload: hivResistancePayload,
        mutationPanels: hivMutationPanels,
      });
    }
    const groupTypingNode = document.getElementById("rsv-typing-group-table");
    if (groupTypingNode && isRotavirus && Array.isArray(groupTyping?.rows) && groupTyping.rows.length) {
      groupTypingNode.dataset.exportTitle = "Rotavirus_A_B_C_大组覆盖度比较";
      renderInteractiveContigTable(
        groupTypingNode,
        Array.isArray(groupTyping?.columns) ? groupTyping.columns : [],
        Array.isArray(groupTyping?.rows) ? groupTyping.rows : [],
        "rsv-typing-group-table",
      );
    }
    const subtypeTypingNode = document.getElementById("rsv-typing-subtype-table");
    if (subtypeTypingNode && isRotavirus && Array.isArray(subtypeTyping?.rows) && subtypeTyping.rows.length) {
      subtypeTypingNode.dataset.exportTitle = "Rotavirus_VP4_VP7_G_P_组合分型";
      renderInteractiveContigTable(
        subtypeTypingNode,
        Array.isArray(subtypeTyping?.columns) ? subtypeTyping.columns : [],
        Array.isArray(subtypeTyping?.rows) ? subtypeTyping.rows : [],
        "rsv-typing-subtype-table",
      );
    }
    const phfTableNode = document.getElementById("rsv-typing-phf-table");
    if (phfTableNode && isHadv && Array.isArray(phfTable?.rows) && phfTable.rows.length) {
      phfTableNode.dataset.exportTitle = "HAdV_PHF_三基因分型表";
      renderInteractiveContigTable(
        phfTableNode,
        Array.isArray(phfTable?.columns) ? phfTable.columns : [],
        Array.isArray(phfTable?.rows) ? phfTable.rows : [],
        "rsv-typing-phf-table",
      );
    }
    const bandavirusSelectionNode = document.getElementById("rsv-typing-bandavirus-selection-table");
    if (bandavirusSelectionNode && isBandavirus && Array.isArray(bandavirusSelection?.rows) && bandavirusSelection.rows.length) {
      bandavirusSelectionNode.dataset.exportTitle = "Bandavirus_参考筛选结果";
      renderInteractiveContigTable(
        bandavirusSelectionNode,
        Array.isArray(bandavirusSelection?.columns) ? bandavirusSelection.columns : [],
        Array.isArray(bandavirusSelection?.rows) ? bandavirusSelection.rows : [],
        "rsv-typing-bandavirus-selection-table",
      );
    }
    const orthohantavirusSelectionNode = document.getElementById("rsv-typing-orthohantavirus-selection-table");
    if (orthohantavirusSelectionNode && isOrthohantavirus && Array.isArray(orthohantavirusSelection?.rows) && orthohantavirusSelection.rows.length) {
      orthohantavirusSelectionNode.dataset.exportTitle = "Orthohantavirus_参考筛选结果";
      renderInteractiveContigTable(
        orthohantavirusSelectionNode,
        Array.isArray(orthohantavirusSelection?.columns) ? orthohantavirusSelection.columns : [],
        Array.isArray(orthohantavirusSelection?.rows) ? orthohantavirusSelection.rows : [],
        "rsv-typing-orthohantavirus-selection-table",
      );
    }
    const orthoebolavirusSelectionNode = document.getElementById("rsv-typing-ebola-selection-table");
    if (orthoebolavirusSelectionNode && isEbola && Array.isArray(orthoebolavirusSelection?.rows) && orthoebolavirusSelection.rows.length) {
      orthoebolavirusSelectionNode.dataset.exportTitle = "Orthoebolavirus_参考筛选结果";
      renderInteractiveContigTable(
        orthoebolavirusSelectionNode,
        Array.isArray(orthoebolavirusSelection?.columns) ? orthoebolavirusSelection.columns : [],
        Array.isArray(orthoebolavirusSelection?.rows) ? orthoebolavirusSelection.rows : [],
        "rsv-typing-ebola-selection-table",
      );
    }
    const orthohantavirusBroadNode = document.getElementById("rsv-typing-orthohantavirus-broad-table");
    if (orthohantavirusBroadNode && isOrthohantavirus && Array.isArray(orthohantavirusBroadTyping?.rows) && orthohantavirusBroadTyping.rows.length) {
      orthohantavirusBroadNode.dataset.exportTitle = "Orthohantavirus_Broad_分型支持摘要";
      renderInteractiveContigTable(
        orthohantavirusBroadNode,
        Array.isArray(orthohantavirusBroadTyping?.columns) ? orthohantavirusBroadTyping.columns : [],
        Array.isArray(orthohantavirusBroadTyping?.rows) ? orthohantavirusBroadTyping.rows : [],
        "rsv-typing-orthohantavirus-broad-table",
      );
    }
    const orthohantavirusSegmentNode = document.getElementById("rsv-typing-orthohantavirus-segments-table");
    if (orthohantavirusSegmentNode && isOrthohantavirus && Array.isArray(orthohantavirusSegments?.rows) && orthohantavirusSegments.rows.length) {
      const rawColumns = Array.isArray(orthohantavirusSegments?.columns) ? orthohantavirusSegments.columns : [];
      const preferredColumns = [
        "segment",
        "typed_label",
        "accession",
        "coverage",
        "mean_depth",
        "covered_bases",
        "num_reads",
        "header",
      ].filter((column) => rawColumns.includes(column));
      const displayColumns = preferredColumns.length ? preferredColumns : rawColumns;
      const displayRows = preferredColumns.length
        ? orthohantavirusSegments.rows.map((row) => {
          if (!Array.isArray(row)) return displayColumns.map(() => "");
          return displayColumns.map((column) => row[rawColumns.indexOf(column)]);
        })
        : orthohantavirusSegments.rows;
      orthohantavirusSegmentNode.dataset.exportTitle = "Orthohantavirus_LMS_最优参考";
      renderInteractiveContigTable(
        orthohantavirusSegmentNode,
        displayColumns,
        displayRows,
        "rsv-typing-orthohantavirus-segments-table",
      );
    }
    const bandavirusAfNode = document.getElementById("rsv-typing-bandavirus-af-table");
    if (bandavirusAfNode && isBandavirus && Array.isArray(bandavirusAfSegments?.rows) && bandavirusAfSegments.rows.length) {
      bandavirusAfNode.dataset.exportTitle = "Bandavirus_A_F_LMS_最优参考";
      renderInteractiveContigTable(
        bandavirusAfNode,
        Array.isArray(bandavirusAfSegments?.columns) ? bandavirusAfSegments.columns : [],
        Array.isArray(bandavirusAfSegments?.rows) ? bandavirusAfSegments.rows : [],
        "rsv-typing-bandavirus-af-table",
      );
    }
    const bandavirusCjNode = document.getElementById("rsv-typing-bandavirus-cj-table");
    if (bandavirusCjNode && isBandavirus && Array.isArray(bandavirusCjSegments?.rows) && bandavirusCjSegments.rows.length) {
      bandavirusCjNode.dataset.exportTitle = "Bandavirus_CJ_LMS_分型结果";
      renderInteractiveContigTable(
        bandavirusCjNode,
        Array.isArray(bandavirusCjSegments?.columns) ? bandavirusCjSegments.columns : [],
        Array.isArray(bandavirusCjSegments?.rows) ? bandavirusCjSegments.rows : [],
        "rsv-typing-bandavirus-cj-table",
      );
    }
    const bandavirusConsensusNode = document.getElementById("rsv-typing-bandavirus-consensus-table");
    if (bandavirusConsensusNode && isBandavirus && Array.isArray(consensusTyping?.rows) && consensusTyping.rows.length) {
      bandavirusConsensusNode.dataset.exportTitle = "Bandavirus_Consensus_分型复核";
      renderInteractiveContigTable(
        bandavirusConsensusNode,
        Array.isArray(consensusTyping?.columns) ? consensusTyping.columns : [],
        Array.isArray(consensusTyping?.rows) ? consensusTyping.rows : [],
        "rsv-typing-bandavirus-consensus-table",
      );
    }
    const orthohantavirusConsensusNode = document.getElementById("rsv-typing-orthohantavirus-consensus-table");
    if (orthohantavirusConsensusNode && isOrthohantavirus && Array.isArray(consensusTyping?.rows) && consensusTyping.rows.length) {
      orthohantavirusConsensusNode.dataset.exportTitle = "Orthohantavirus_Consensus_分型复核";
      renderInteractiveContigTable(
        orthohantavirusConsensusNode,
        Array.isArray(consensusTyping?.columns) ? consensusTyping.columns : [],
        Array.isArray(consensusTyping?.rows) ? consensusTyping.rows : [],
        "rsv-typing-orthohantavirus-consensus-table",
      );
    }
    const phfCoverageChartNode = document.getElementById("rsv-typing-phf-coverage-charts");
    if (phfCoverageChartNode && isHadv && Array.isArray(phfCoverage?.rows) && phfCoverage.rows.length) {
      const xValues = Array.isArray(phfCoverage?.x_values) ? phfCoverage.x_values : [];
      const coveragePoints = Array.isArray(phfCoverage?.coverage_points) ? phfCoverage.coverage_points : [];
      const depthPoints = Array.isArray(phfCoverage?.depth_points) ? phfCoverage.depth_points : [];
      phfCoverageChartNode.innerHTML = `
        <div class="influenza-segment-coverage-grid">
          ${renderBarSvg(coveragePoints, {
            label: "PHF 三基因覆盖度",
            width: 560,
            height: 320,
            padX: 52,
            padTop: 24,
            padBottom: 84,
            xLabel: "分型基因",
            yLabel: "覆盖度(%)",
            xValues,
          })}
          ${renderBarSvg(depthPoints, {
            label: "PHF 三基因平均深度",
            width: 560,
            height: 320,
            padX: 52,
            padTop: 24,
            padBottom: 84,
            xLabel: "分型基因",
            yLabel: "平均深度",
            xValues,
          })}
        </div>
      `;
    }
    const phfCoverageTableNode = document.getElementById("rsv-typing-phf-coverage-table");
    if (phfCoverageTableNode && isHadv && Array.isArray(phfCoverage?.rows) && phfCoverage.rows.length) {
      phfCoverageTableNode.dataset.exportTitle = "HAdV_PHF_三基因覆盖度汇总表";
      renderInteractiveContigTable(
        phfCoverageTableNode,
        Array.isArray(phfCoverage?.columns) ? phfCoverage.columns : [],
        Array.isArray(phfCoverage?.rows) ? phfCoverage.rows : [],
        "rsv-typing-phf-coverage-table",
      );
    }
    const phfSnpChartNode = document.getElementById("rsv-typing-phf-snp-chart");
    if (phfSnpChartNode && isHadv && Array.isArray(phfSnp?.rows) && phfSnp.rows.length) {
      phfSnpChartNode.innerHTML = renderBarSvg(
        Array.isArray(phfSnp?.values) ? phfSnp.values : [],
        {
          label: "PHF 三基因差异 SNP 数",
          width: 1120,
          height: 360,
          padX: 52,
          padTop: 24,
          padBottom: 84,
          xLabel: "分型基因",
          yLabel: "差异SNP数",
          xValues: Array.isArray(phfSnp?.x_values) ? phfSnp.x_values : [],
        },
      );
    }
    const phfSnpTableNode = document.getElementById("rsv-typing-phf-snp-table");
    if (phfSnpTableNode && isHadv && Array.isArray(phfSnp?.rows) && phfSnp.rows.length) {
      phfSnpTableNode.dataset.exportTitle = "HAdV_PHF_三基因差异SNP数";
      renderInteractiveContigTable(
        phfSnpTableNode,
        Array.isArray(phfSnp?.columns) ? phfSnp.columns : [],
        Array.isArray(phfSnp?.rows) ? phfSnp.rows : [],
        "rsv-typing-phf-snp-table",
      );
    }
    const phfSnpDetailTableNode = document.getElementById("rsv-typing-phf-snp-detail-table");
    if (phfSnpDetailTableNode && isHadv && Array.isArray(phfSnp?.detail_rows) && phfSnp.detail_rows.length) {
      phfSnpDetailTableNode.dataset.exportTitle = "HAdV_PHF_三基因差异SNP明细";
      renderInteractiveContigTable(
        phfSnpDetailTableNode,
        Array.isArray(phfSnp?.detail_columns) ? phfSnp.detail_columns : [],
        Array.isArray(phfSnp?.detail_rows) ? phfSnp.detail_rows : [],
        "rsv-typing-phf-snp-detail-table",
      );
    }
    const mutationTableNode = document.getElementById("rsv-typing-mutation-table");
    if (mutationTableNode) {
      const mutationColumns = Array.isArray(mutationTable?.columns) ? mutationTable.columns : [];
      const mutationRows = Array.isArray(mutationTable?.rows) ? mutationTable.rows : [];
      const qualityIndex = mutationColumns.indexOf("质量分层");
      const highRows = qualityIndex >= 0 ? mutationRows.filter((row) => String(row[qualityIndex] ?? "").trim() === "高质量突变") : mutationRows;
      const lowRows = qualityIndex >= 0 ? mutationRows.filter((row) => String(row[qualityIndex] ?? "").trim() === "低质量突变") : [];
      const renderMutationTab = (tabKey = "high") => {
        const activeRows = tabKey === "low" ? lowRows : highRows;
        mutationTableNode.dataset.exportTitle = tabKey === "low" ? `${virusShort}_低质量突变位点表` : `${virusShort}_高质量突变位点表`;
        renderInteractiveContigTable(
          mutationTableNode,
          mutationColumns,
          activeRows,
          "rsv-typing-mutation-table",
        );
        document.querySelectorAll("[data-rsv-variant-tab]").forEach((button) => {
          button.classList.toggle("active", button.getAttribute("data-rsv-variant-tab") === tabKey);
        });
      };
      renderMutationTab("high");
      document.querySelectorAll("[data-rsv-variant-tab]").forEach((button) => {
        button.addEventListener("click", () => {
          renderMutationTab(String(button.getAttribute("data-rsv-variant-tab") || "high"));
        });
      });
    }
    if ((isNorovirus || isEnterovirus || isAstroviridae || isRhinovirus || isSeasonalHcov) && typingPhylogenyTrees.length) {
      typingPhylogenyTrees.forEach((treeSection) => {
        const geneKey = String(treeSection?.gene || "").trim().toLowerCase();
        const containerId = geneKey === "vp1"
          ? "rsv-phylogeny-vp1"
          : (geneKey === "orf2" ? "rsv-phylogeny-vp1" : (geneKey === "rdrp" ? "rsv-phylogeny-rdrp" : (geneKey === "s" ? "rsv-phylogeny-vp1" : "")));
        if (!containerId) return;
        renderNewickTreeCard(containerId, {
          ...treeSection,
          task_id: currentReportData?.task?.id || "",
          label: treeSection?.label || `${virusLabel} ${String(treeSection?.gene || "").trim().toUpperCase()} 系统发育树`,
        }, { rowHeight: 16, drawingWidth: 1120, labelsWidth: 240 });
      });
    } else if (!isHpiv && !isEnterovirus && !isAstroviridae && !isRhinovirus && !isSeasonalHcov && !isRotavirus && String(section?.phylogeny_tree?.status || "") === "ready") {
      renderITOLTreeCard("rsv-phylogeny-tree", {
        ...section.phylogeny_tree,
        task_id: currentReportData?.task?.id || "",
        label: section.phylogeny_tree?.label || `${virusLabel} 系统发育树`,
      });
    }
    const nmdcTableNode = document.getElementById("rsv-typing-nmdc-table");
    if (nmdcTableNode && String(nmdcAnnotation?.status || "") === "ready") {
      nmdcTableNode.dataset.exportTitle = "RSV_NMDC_变异风险数据库注释";
      renderInteractiveContigTable(
        nmdcTableNode,
        Array.isArray(nmdcAnnotation?.columns) ? nmdcAnnotation.columns : [],
        Array.isArray(nmdcAnnotation?.rows) ? nmdcAnnotation.rows : [],
        "rsv-typing-nmdc-table",
      );
    }
    const functionalImpactTableNode = document.getElementById("rsv-typing-functional-impact-table");
    if (functionalImpactTableNode && isHpiv && String(functionalAnnotation?.status || "") === "ready") {
      functionalImpactTableNode.dataset.exportTitle = "HPIV_突变功能影响注释";
      renderInteractiveContigTable(
        functionalImpactTableNode,
        Array.isArray(functionalAnnotation?.columns) ? functionalAnnotation.columns : [],
        Array.isArray(functionalAnnotation?.rows) ? functionalAnnotation.rows : [],
        "rsv-typing-functional-impact-table",
      );
    }
    const igvFrame = document.getElementById("rsv-igv-frame");
    if (igvFrame && String(igvView?.status || "") === "ready") {
      const firstMutation = Array.isArray(mutationTable?.rows) && mutationTable.rows.length ? mutationTable.rows[0] : null;
      const chromIndex = Array.isArray(mutationTable?.columns) ? mutationTable.columns.indexOf("染色体") : -1;
      const posIndex = Array.isArray(mutationTable?.columns) ? mutationTable.columns.indexOf("位置") : -1;
      const initialLocus = firstMutation && chromIndex >= 0 && posIndex >= 0
        ? `${String(firstMutation[chromIndex] || "").trim()}:${String(firstMutation[posIndex] || "").trim()}`
        : "";
      initializeDeferredIgvEmbed({
        frameId: "rsv-igv-frame",
        panelId: "rsv-igv-lazy",
        buttonId: "rsv-igv-load",
        task: currentReportData?.task || {},
        igvView,
        initialLocus,
        mutationTableNode,
      });
    }
    return;
  }
  if (mode === "monkeypox_nextclade") {
    const sectionNode = document.getElementById("section-serotype");
    const headingTitle = sectionNode?.querySelector("h2");
    const headingCopy = sectionNode?.querySelector(".section-heading p:last-child");
    const hasFastp = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "") === "ready";
    const hasSpecies = Array.isArray(currentReportData?.sections?.species_identification?.species?.rows)
      && currentReportData.sections.species_identification.species.rows.length > 0;
    const hasSubspecies = Array.isArray(currentReportData?.sections?.species_identification?.subspecies?.rows)
      && currentReportData.sections.species_identification.subspecies.rows.length > 0;
    const hasCoverage = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
      && currentReportData.sections.assembly.coverage.points.length > 0;
    const hasAssemblySummary = Array.isArray(currentReportData?.sections?.assembly?.summary?.rows)
      && currentReportData.sections.assembly.summary.rows.length > 0;
    const keepVisible = new Set(["section-raw-qc", "section-serotype"]);
    if (hasFastp) keepVisible.add("section-fastp");
    if (hasSpecies || hasSubspecies) keepVisible.add("section-species-identification");
    if (hasCoverage) keepVisible.add("section-assembly");
    if (hasAssemblySummary) keepVisible.add("section-assembly-summary");
    [
      "section-raw-qc",
      "section-fastp",
      "section-species",
      "section-species-identification",
      "section-taxonomy-abundance",
      "section-assembly",
      "section-assembly-summary",
      "section-contig-annotation",
      "section-cgview",
      "section-checkm",
      "section-gene-annotation",
      "section-rv",
      "section-rv-summary",
      "section-virulence",
      "section-resistance",
      "section-resistance-mutation",
      "section-mlst",
      "section-priority-serotype",
      "section-mge",
      "section-mge-resistance",
      "section-mge-virulence",
    ].forEach((id) => {
      const element = document.getElementById(id);
      if (!element) return;
      element.classList.toggle("hidden", !keepVisible.has(id));
    });
    if (headingTitle) headingTitle.textContent = "猴痘分型";
    if (headingCopy) headingCopy.textContent = "基于 hMPXV Nextclade 数据集展示 clade、lineage、outbreak、质控状态以及参考比对结果。";
    applyTableTone(container, "serotype-table");
    const summaryCards = Array.isArray(section?.summary_cards) ? section.summary_cards : [];
    const qualityMetrics = Array.isArray(section?.quality_metrics) ? section.quality_metrics : [];
    const notes = String(section?.notes || "").trim();
    const columns = Array.isArray(section?.columns) ? section.columns : [];
    const rows = Array.isArray(section?.rows) ? section.rows : [];
    const mutationTable = section?.mutation_table && typeof section.mutation_table === "object"
      ? section.mutation_table
      : { status: "empty", columns: [], rows: [] };
    const igvView = section?.igv && typeof section.igv === "object"
      ? section.igv
      : { status: "empty" };
    const mutationSummaryCards = [
      { label: "总变异位点", value: String(mutationTable?.total_variants ?? (Array.isArray(mutationTable?.rows) ? mutationTable.rows.length : "--")) },
      { label: "高质量突变", value: String(mutationTable?.high_quality_variants ?? "--") },
      { label: "低质量突变", value: String(mutationTable?.low_quality_variants ?? "--") },
      { label: "注释VCF", value: String(mutationTable?.source_vcf ? "snps.anno.vcf" : "--") },
    ];
    const mutationSummaryMarkup = mutationSummaryCards.map((item) => `
      <article class="mini-stat-card">
        <span>${escapeHtml(item.label)}</span>
        <strong>${escapeHtml(item.value)}</strong>
      </article>
    `).join("");
    const summaryMarkup = summaryCards.length ? `
      <div class="mini-stat-grid">
        ${summaryCards.map((item) => `
          <article class="mini-stat-card">
            <span>${escapeHtml(item?.label || "--")}</span>
            <strong>${escapeHtml(String(item?.value ?? "--"))}</strong>
          </article>
        `).join("")}
      </div>
    ` : "";
    const qualityMarkup = qualityMetrics.length ? `
      <div class="mini-stat-grid">
        ${qualityMetrics.map((item) => `
          <article class="mini-stat-card">
            <span>${escapeHtml(item?.label || "--")}</span>
            <strong>${escapeHtml(String(item?.value ?? "--"))}</strong>
          </article>
        `).join("")}
      </div>
    ` : "";
    container.dataset.exportTitle = "猴痘 Nextclade 分型";
    container.innerHTML = `
      <div class="serotype-special-layout influenza-typing-layout">
        ${summaryMarkup}
        ${notes ? `
          <div class="chart-insight serotype-insight" role="note" aria-label="猴痘分型说明">
            <span class="chart-insight-label">结果说明</span>
            <p>${escapeHtml(notes)}</p>
          </div>
        ` : ""}
        <section id="monkeypox-typing-summary" class="result-card" data-report-nav-anchor>
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">Nextclade summary</span>
              <h3>猴痘分型总表</h3>
            </div>
          </div>
          <div id="monkeypox-typing-summary-table" class="report-table-card report-table-card-embedded"></div>
        </section>
        ${qualityMarkup ? `
          <section id="monkeypox-typing-quality" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">QC metrics</span>
                <h3>Nextclade 质量指标</h3>
              </div>
            </div>
            ${qualityMarkup}
          </section>
        ` : ""}
        ${Array.isArray(mutationTable?.rows) && mutationTable.rows.length ? `
          <section id="monkeypox-typing-mutations" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Variants</span>
                <h3>突变位点表</h3>
              </div>
              <span class="card-tag">${escapeHtml(`${mutationTable.rows.length} 条`)}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">读取 <code>snps.anno.vcf</code> 中的 SnpEff 注释结果，提取猴痘样本的关键突变位点用于表格浏览和 IGV 联动。</p>
            <div id="monkeypox-typing-mutation-summary" class="mini-stat-grid">${mutationSummaryMarkup}</div>
            <div class="nextclade-variant-tabs" id="monkeypox-variant-tabs" role="tablist" aria-label="猴痘变异质量分层切换">
              <button type="button" class="report-tab-button active" data-monkeypox-variant-tab="high">高质量突变</button>
              <button type="button" class="report-tab-button" data-monkeypox-variant-tab="low">低质量突变</button>
            </div>
            <div id="monkeypox-typing-mutation-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${String(igvView?.status || "") === "ready" ? `
          <section id="monkeypox-typing-igv" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">IGV</span>
                <h3>IGV 比对结果</h3>
              </div>
              <span class="card-tag">${escapeHtml(String(igvView?.viewer_label || "参考比对视图"))}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">${escapeHtml(String(igvView?.note || "展示猴痘样本的参考比对与注释轨道。"))}</p>
            <div class="report-igv-shell">
              <div id="monkeypox-igv-lazy" class="empty-box">
                <p>IGV 改为按需加载，避免页面初始卡顿。点击下方按钮或上方突变位点后会自动开始加载。</p>
                <button type="button" id="monkeypox-igv-load" class="table-export-button">加载 IGV</button>
              </div>
              <iframe id="monkeypox-igv-frame" class="report-igv-frame" title="猴痘 IGV 比对结果" hidden loading="lazy"></iframe>
            </div>
          </section>
        ` : ""}
      </div>
    `;
    const summaryTable = document.getElementById("monkeypox-typing-summary-table");
    if (summaryTable) {
      const preferredColumns = [
        "seqName",
        "clade",
        "lineage",
        "outbreak",
        "qc.overallStatus",
        "qc.overallScore",
        "coverage",
        "totalSubstitutions",
        "totalAminoacidSubstitutions",
        "totalDeletions",
        "totalInsertions",
        "totalFrameShifts",
      ].filter((column) => columns.includes(column));
      const indexes = preferredColumns.map((column) => columns.indexOf(column)).filter((index) => index >= 0);
      const displayColumns = preferredColumns.length ? preferredColumns : columns;
      const displayRows = preferredColumns.length
        ? rows.map((row) => indexes.map((index) => Array.isArray(row) ? row[index] : ""))
        : rows;
      summaryTable.dataset.exportTitle = "猴痘_Nextclade_分型总表";
      renderInteractiveContigTable(summaryTable, displayColumns, displayRows, "monkeypox-typing-summary-table");
    }
    const mutationTableNode = document.getElementById("monkeypox-typing-mutation-table");
    if (mutationTableNode) {
      const mutationColumns = Array.isArray(mutationTable?.columns) ? mutationTable.columns : [];
      const mutationRows = Array.isArray(mutationTable?.rows) ? mutationTable.rows : [];
      const qualityIndex = mutationColumns.indexOf("质量分层");
      const highRows = qualityIndex >= 0 ? mutationRows.filter((row) => String(row[qualityIndex] ?? "").trim() === "高质量突变") : mutationRows;
      const lowRows = qualityIndex >= 0 ? mutationRows.filter((row) => String(row[qualityIndex] ?? "").trim() === "低质量突变") : [];
      const renderMutationTab = (tabKey = "high") => {
        const activeRows = tabKey === "low" ? lowRows : highRows;
        mutationTableNode.dataset.exportTitle = tabKey === "low" ? "猴痘_低质量突变位点表" : "猴痘_高质量突变位点表";
        renderInteractiveContigTable(
          mutationTableNode,
          mutationColumns,
          activeRows,
          "monkeypox-typing-mutation-table",
        );
        document.querySelectorAll("[data-monkeypox-variant-tab]").forEach((button) => {
          button.classList.toggle("active", button.getAttribute("data-monkeypox-variant-tab") === tabKey);
        });
      };
      renderMutationTab("high");
      document.querySelectorAll("[data-monkeypox-variant-tab]").forEach((button) => {
        button.addEventListener("click", () => {
          renderMutationTab(String(button.getAttribute("data-monkeypox-variant-tab") || "high"));
        });
      });
    }
    const igvFrame = document.getElementById("monkeypox-igv-frame");
    if (igvFrame && String(igvView?.status || "") === "ready") {
      const firstMutation = Array.isArray(mutationTable?.rows) && mutationTable.rows.length ? mutationTable.rows[0] : null;
      const chromIndex = Array.isArray(mutationTable?.columns) ? mutationTable.columns.indexOf("染色体") : -1;
      const posIndex = Array.isArray(mutationTable?.columns) ? mutationTable.columns.indexOf("位置") : -1;
      const initialLocus = firstMutation && chromIndex >= 0 && posIndex >= 0
        ? `${String(firstMutation[chromIndex] || "").trim()}:${String(firstMutation[posIndex] || "").trim()}`
        : "";
      initializeDeferredIgvEmbed({
        frameId: "monkeypox-igv-frame",
        panelId: "monkeypox-igv-lazy",
        buttonId: "monkeypox-igv-load",
        task: currentReportData?.task || {},
        igvView,
        initialLocus,
        mutationTableNode,
      });
    }
    return;
  }
  if (mode === "influenza_typing") {
    const sectionNode = document.getElementById("section-serotype");
    const headingTitle = sectionNode?.querySelector("h2");
    const headingCopy = sectionNode?.querySelector(".section-heading p:last-child");
    const hasFastp = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "") === "ready";
    const hasSpecies = Array.isArray(currentReportData?.sections?.species_identification?.species?.rows)
      && currentReportData.sections.species_identification.species.rows.length > 0;
    const hasSubspecies = Array.isArray(currentReportData?.sections?.species_identification?.subspecies?.rows)
      && currentReportData.sections.species_identification.subspecies.rows.length > 0;
    const hasCoverage = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
      && currentReportData.sections.assembly.coverage.points.length > 0;
    const hasAssemblySummary = Array.isArray(currentReportData?.sections?.assembly?.summary?.rows)
      && currentReportData.sections.assembly.summary.rows.length > 0;
    const keepVisible = new Set(["section-raw-qc", "section-serotype"]);
    if (hasFastp) keepVisible.add("section-fastp");
    if (hasSpecies || hasSubspecies) keepVisible.add("section-species-identification");
    if (hasCoverage) keepVisible.add("section-assembly");
    if (hasAssemblySummary) keepVisible.add("section-assembly-summary");
    [
      "section-raw-qc",
      "section-fastp",
      "section-species",
      "section-species-identification",
      "section-taxonomy-abundance",
      "section-assembly",
      "section-assembly-summary",
      "section-contig-annotation",
      "section-cgview",
      "section-checkm",
      "section-gene-annotation",
      "section-rv",
      "section-rv-summary",
      "section-virulence",
      "section-resistance",
      "section-resistance-mutation",
      "section-mlst",
      "section-priority-serotype",
      "section-mge",
      "section-mge-resistance",
      "section-mge-virulence",
    ].forEach((id) => {
      const element = document.getElementById(id);
      if (!element) return;
      element.classList.toggle("hidden", !keepVisible.has(id));
    });
    if (headingTitle) headingTitle.textContent = "流感分型";
    if (headingCopy) headingCopy.textContent = "基于 wf_flu 风格的 IRMA reference set 初筛与 HA/NA 最优亚型选择，展示甲/乙流判断、亚型组合和最终 8 segment 参考组成。";
    applyTableTone(container, "serotype-table");
    const summaryCards = Array.isArray(section?.summary_cards) ? section.summary_cards : [];
    const summaryMarkup = summaryCards.length ? `
      <div class="mini-stat-grid">
        ${summaryCards.map((item) => `
          <article class="mini-stat-card">
            <span>${escapeHtml(item?.label || "--")}</span>
            <strong>${escapeHtml(String(item?.value ?? "--"))}</strong>
          </article>
        `).join("")}
      </div>
    ` : "";
    const notes = String(section?.notes || "").trim();
    const summaryColumns = Array.isArray(section?.columns) ? section.columns : [];
    const summaryRows = Array.isArray(section?.rows) ? section.rows : [];
    const segmentManifest = section?.segment_manifest && typeof section.segment_manifest === "object"
      ? section.segment_manifest
      : { columns: [], rows: [] };
    const mutationTable = section?.mutation_table && typeof section.mutation_table === "object"
      ? section.mutation_table
      : { columns: [], rows: [] };
    const mutationSummary = section?.mutation_summary && typeof section.mutation_summary === "object"
      ? section.mutation_summary
      : {};
    const variantAnnotation = section?.variant_annotation && typeof section.variant_annotation === "object"
      ? section.variant_annotation
      : { status: "empty", columns: [], rows: [] };
    const resistanceAnnotation = section?.resistance_annotation && typeof section.resistance_annotation === "object"
      ? section.resistance_annotation
      : { status: "empty", columns: [], rows: [], total_hits: 0 };
    const igvView = section?.igv && typeof section.igv === "object"
      ? section.igv
      : { status: "empty" };
    const variantSummaryCards = [
      { label: "总变异位点", value: String(variantAnnotation?.total_variants ?? mutationSummary?.count ?? "--") },
      { label: "高质量突变", value: String(variantAnnotation?.high_quality_variants ?? "--") },
      { label: "低质量突变", value: String(variantAnnotation?.low_quality_variants ?? "--") },
      { label: "注释VCF", value: String(variantAnnotation?.source_vcf ? "snps.filt1.snpeff.vcf" : "--") },
    ];
    const variantSummaryMarkup = variantSummaryCards.map((item) => `
      <article class="mini-stat-card">
        <span>${escapeHtml(item.label)}</span>
        <strong>${escapeHtml(item.value)}</strong>
      </article>
    `).join("");
    container.dataset.exportTitle = "流感分型";
    container.innerHTML = `
      <div class="serotype-special-layout influenza-typing-layout">
        ${summaryMarkup}
        ${notes ? `
          <div class="chart-insight serotype-insight" role="note" aria-label="流感分型说明">
            <span class="chart-insight-label">结果说明</span>
            <p>${escapeHtml(notes)}</p>
          </div>
        ` : ""}
        <section id="influenza-typing-summary" class="result-card" data-report-nav-anchor>
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">Typing summary</span>
              <h3>流感分型总表</h3>
            </div>
          </div>
          <div id="influenza-typing-summary-table" class="report-table-card report-table-card-embedded"></div>
        </section>
        <section id="influenza-typing-manifest" class="result-card" data-report-nav-anchor>
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">Reference set</span>
              <h3>最终 8 Segment 参考组成</h3>
            </div>
          </div>
          <div id="influenza-typing-segment-table" class="report-table-card report-table-card-embedded"></div>
        </section>
        ${Array.isArray(mutationTable?.rows) && mutationTable.rows.length ? `
          <section id="influenza-typing-mutations" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">SnpEff</span>
                <h3>变异注释表</h3>
              </div>
              <span class="card-tag">${escapeHtml(`${mutationSummary?.count ?? mutationTable.rows.length} 条`)}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">读取流感样本的 <code>snps.filt1.vcf</code>，结合 VADR 生成的 GFF3 与 consensus FASTA 构建 snpEff 数据库，并输出可直接筛查的变异注释结果。</p>
            <div id="influenza-typing-mutation-summary" class="mini-stat-grid">${variantSummaryMarkup}</div>
            <div id="influenza-typing-mutation-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
        ${String(igvView?.status || "") === "ready" ? `
          <section id="influenza-typing-igv" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">IGV</span>
                <h3>IGV 比对结果</h3>
              </div>
              <span class="card-tag">${escapeHtml(String(igvView?.viewer_label || "参考比对视图"))}</span>
            </div>
            <p class="nextclade-variant-annotation-copy">${escapeHtml(String(igvView?.note || "点击上方位点后，IGV 会自动跳转到对应位置。"))}</p>
            <div class="report-igv-shell">
              <div id="influenza-igv-lazy" class="empty-box">
                <p>IGV 改为按需加载，避免页面初始卡顿。点击下方按钮或上方突变位点后会自动开始加载。</p>
                <button type="button" id="influenza-igv-load" class="table-export-button">加载 IGV</button>
              </div>
              <iframe id="influenza-igv-frame" class="report-igv-frame" title="流感 IGV 比对结果" hidden loading="lazy"></iframe>
            </div>
          </section>
        ` : ""}
        ${String(resistanceAnnotation?.status || "") === "ready" ? `
          <section id="influenza-typing-resistance" class="result-card" data-report-nav-anchor>
            <div class="card-head">
              <div class="card-title-stack">
                <span class="section-chip">Resistance</span>
                <h3>耐药突变注释结果</h3>
              </div>
              <span class="card-tag">${escapeHtml(String(resistanceAnnotation?.total_hits ?? (Array.isArray(resistanceAnnotation?.rows) ? resistanceAnnotation.rows.length : 0)))} 条</span>
            </div>
            <p class="nextclade-variant-annotation-copy">基于 <code>influenza_resistance_rules.tsv</code> 对高质量流感氨基酸突变进行耐药规则匹配，输出可直接判读的耐药位点结果。</p>
            <div id="influenza-typing-resistance-table" class="report-table-card report-table-card-embedded"></div>
          </section>
        ` : ""}
      </div>
    `;
    const summaryTable = document.getElementById("influenza-typing-summary-table");
    if (summaryTable) {
      summaryTable.dataset.exportTitle = "流感分型总表";
      renderInteractiveContigTable(summaryTable, summaryColumns, summaryRows, "influenza-typing-summary-table");
    }
    const segmentTable = document.getElementById("influenza-typing-segment-table");
    if (segmentTable) {
      segmentTable.dataset.exportTitle = "流感_8_segment_参考组成";
      renderInteractiveContigTable(
        segmentTable,
        Array.isArray(segmentManifest?.columns) ? segmentManifest.columns : [],
        Array.isArray(segmentManifest?.rows) ? segmentManifest.rows : [],
        "influenza-typing-segment-table",
      );
    }
    const mutationTableNode = document.getElementById("influenza-typing-mutation-table");
    if (mutationTableNode) {
      mutationTableNode.dataset.exportTitle = "流感_snpEff_变异注释表";
      renderInteractiveContigTable(
        mutationTableNode,
        Array.isArray(mutationTable?.columns) ? mutationTable.columns : [],
        Array.isArray(mutationTable?.rows) ? mutationTable.rows : [],
        "influenza-typing-mutation-table",
      );
    }
    const igvFrame = document.getElementById("influenza-igv-frame");
    if (igvFrame && String(igvView?.status || "") === "ready") {
      const firstMutation = Array.isArray(mutationTable?.rows) && mutationTable.rows.length ? mutationTable.rows[0] : null;
      const chromIndex = Array.isArray(mutationTable?.columns) ? mutationTable.columns.indexOf("染色体") : -1;
      const posIndex = Array.isArray(mutationTable?.columns) ? mutationTable.columns.indexOf("位置") : -1;
      const initialLocus = firstMutation && chromIndex >= 0 && posIndex >= 0
        ? `${String(firstMutation[chromIndex] || "").trim()}:${String(firstMutation[posIndex] || "").trim()}`
        : "";
      initializeDeferredIgvEmbed({
        frameId: "influenza-igv-frame",
        panelId: "influenza-igv-lazy",
        buttonId: "influenza-igv-load",
        task: currentReportData?.task || {},
        igvView,
        initialLocus,
        mutationTableNode,
      });
    }
    const resistanceTable = document.getElementById("influenza-typing-resistance-table");
    if (resistanceTable) {
      resistanceTable.dataset.exportTitle = "流感_耐药突变注释结果";
      renderInteractiveContigTable(
        resistanceTable,
        Array.isArray(resistanceAnnotation?.columns) ? resistanceAnnotation.columns : [],
        Array.isArray(resistanceAnnotation?.rows) ? resistanceAnnotation.rows : [],
        "influenza-typing-resistance-table",
      );
    }
    return;
  }
  if (mode === "tb_profiler") {
    applyTableTone(container, "serotype-table");
    const sectionNode = document.getElementById("section-serotype");
    const headingTitle = sectionNode?.querySelector("h2");
    const headingCopy = sectionNode?.querySelector(".section-heading p:last-child");
    if (headingTitle) headingTitle.textContent = "结核家系分析";
    if (headingCopy) headingCopy.textContent = "展示基于 tb-profiler 的结核分枝杆菌家系结果，并关联知识库中的家系背景、地区分布与关键标记提示。";
    const summaryCards = Array.isArray(section?.summary_cards) ? section.summary_cards : [];
    const notes = String(section?.notes || "").trim();
    const columns = Array.isArray(section?.columns) ? section.columns : [];
    const knowledgeSummary = section?.knowledge_summary && typeof section.knowledge_summary === "object"
      ? section.knowledge_summary
      : {};
    const knowledgeHeadline = String(knowledgeSummary?.headline || "").trim();
    const knowledgeItems = Array.isArray(knowledgeSummary?.items) ? knowledgeSummary.items : [];
    const summaryMarkup = summaryCards.length ? `
      <div class="mini-stat-grid">
        ${summaryCards.map((item) => `
          <article class="mini-stat-card">
            <span>${escapeHtml(item?.label || "--")}</span>
            <strong>${escapeHtml(String(item?.value ?? "--"))}</strong>
          </article>
        `).join("")}
      </div>
    ` : "";
    const knowledgeMarkup = (knowledgeHeadline || knowledgeItems.length) ? `
      <section class="result-card" data-report-nav-anchor>
        <div class="card-head">
          <div class="card-title-stack">
            <span class="section-chip">Knowledge Base</span>
            <h3>结核家系知识库关联</h3>
          </div>
        </div>
        ${knowledgeHeadline ? `<p class="mlst-knowledge-summary-lead">${escapeHtml(knowledgeHeadline)}</p>` : ""}
        ${knowledgeItems.length ? `
          <div class="knowledge-base-browser-grid">
            ${knowledgeItems.map((item) => `
              <article class="knowledge-base-browser-card">
                <div class="knowledge-base-browser-card-head">
                  <div>
                    <p class="knowledge-base-card-kicker">TB lineage note</p>
                    <h4>${escapeHtml(String(item?.serotype || item?.matched_on || "--"))}</h4>
                    <p class="knowledge-base-browser-subtitle">${escapeHtml(String(item?.panel || "结核分枝杆菌家系知识库"))}</p>
                  </div>
                </div>
                <div class="knowledge-base-browser-body">
                  <dl class="knowledge-base-browser-facts">
                    <div><dt>地区分布</dt><dd>${escapeHtml(Array.isArray(item?.regional) && item.regional.length ? item.regional.join("；") : "-")}</dd></div>
                    <div><dt>关键标记/提示</dt><dd>${escapeHtml(Array.isArray(item?.key_markers) && item.key_markers.length ? item.key_markers.join("；") : "-")}</dd></div>
                    <div><dt>判读提示</dt><dd>${escapeHtml(String(item?.interpretation || "-"))}</dd></div>
                  </dl>
                </div>
              </article>
            `).join("")}
          </div>
        ` : ""}
      </section>
    ` : "";
    container.dataset.exportTitle = "结核家系鉴定";
    container.innerHTML = `
      <div class="serotype-special-layout">
        ${summaryMarkup}
        ${notes ? `
          <div class="chart-insight serotype-insight" role="note" aria-label="结核家系分析说明">
            <span class="chart-insight-label">结果说明</span>
            <p>${escapeHtml(notes)}</p>
          </div>
        ` : ""}
        ${knowledgeMarkup}
      </div>
    `;
    return;
  }
  if (mode !== "bordetella_pertussis") {
    buildTableCard("serotype-table", "血清型鉴定", section?.columns || [], section?.rows || []);
    return;
  }
  const notes = Array.isArray(section?.notes) ? section.notes : [];
  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  const schemeSummary = String(section?.scheme_summary || "").trim();
  applyTableTone(container, "serotype-table");
  const noteMarkup = notes.map((note) => `
    <article class="serotype-note-card serotype-note-${escapeHtml(note?.state || "neutral")}">
      <span class="serotype-note-label">${escapeHtml(note?.label || "结果")}</span>
      <strong>${escapeHtml(note?.text || "--")}</strong>
    </article>
  `).join("");
  const summaryMarkup = schemeSummary ? `
    <div class="chart-insight serotype-insight" role="note" aria-label="抗原基因分型">
      <span class="chart-insight-label">抗原基因分型</span>
      <p>${escapeHtml(schemeSummary)}</p>
    </div>
  ` : "";
  const tableMarkup = rows.length ? `
    <div class="serotype-table-shell">
      ${renderTableExportToolbar()}
      <div class="table-frame">
        <table class="report-table">
          <thead><tr>${columns.map((column) => `<th>${escapeHtml(humanizeReportColumnLabel(column))}</th>`).join("")}</tr></thead>
          <tbody>
            ${rows.map((row) => `<tr>${columns.map((column, index) => {
              const value = Array.isArray(row) ? row[index] : row?.[column];
              return `<td>${renderTableCellContent(value, column)}</td>`;
            }).join("")}</tr>`).join("")}
          </tbody>
        </table>
      </div>
    </div>
  ` : `
    <div class="empty-table-state serotype-empty-state">
      <strong>抗原基因分型暂无表格结果</strong>
      <p class="empty-copy">当前未检出可展示的 {_scheme.tsv} 结果。</p>
    </div>
  `;
  container.dataset.exportTitle = "血清型鉴定";
  container.innerHTML = `
    <div class="serotype-special-layout">
      <div class="serotype-note-grid">${noteMarkup}</div>
      ${summaryMarkup}
      ${tableMarkup}
    </div>
  `;
  if (rows.length) {
    bindTableExportButtons(container, "血清型鉴定", columns, rows.map((row) => (
      Array.isArray(row) ? row : columns.map((column) => row[column] ?? "")
    )));
  }
}

function fillTaskMeta(task) {
  document.getElementById("report-task-name").textContent = task.name || task.id || "-";
  document.getElementById("report-task-meta").textContent = `任务编号：${task.id || "-"}`;
  document.getElementById("report-sample-title").textContent = task.sample_display_name || task.sample_name || task.name || task.id || "分析结果";
  document.getElementById("report-sample-copy").textContent = `创建时间：${formatDateTime(task.created_at)}；开始时间：${formatDateTime(task.started_at)}；结束时间：${formatDateTime(task.finished_at)}。`;
  document.getElementById("meta-owner").textContent = task.owner || "-";
  document.getElementById("meta-group").textContent = task.group || "-";
  document.getElementById("meta-asm-type").textContent = task.asm_type || "-";
  document.getElementById("meta-method").textContent = getTaskMethod(task) || "-";
  document.getElementById("meta-input").textContent = task.input_path || "-";
  document.getElementById("meta-output").textContent = task.output_dir || "-";
  renderSampleSwitcher(task);
}

function applySarsCov2ReportChrome(data) {
  if (!isSarsCov2NextcladeReport(data)) return;
  const task = data?.task || {};
  const backNode = document.querySelector(".report-back");
  const subtitleNode = document.querySelector(".report-subtitle");
  const kickerNode = document.querySelector(".report-title-block .report-kicker");
  const titleNode = document.querySelector(".report-title-block h1");
  const statusNode = document.querySelector(".report-status");
  const sampleTitleNode = document.getElementById("report-sample-title");
  const sampleCopyNode = document.getElementById("report-sample-copy");
  const sampleName = task.sample_display_name || task.sample_name || data?.sections?.serotype?.sequence_name || task.name || task.id || "SARS-CoV-2";
  if (backNode) backNode.textContent = "返回任务";
  if (kickerNode) kickerNode.textContent = "SARS-CoV-2 Report";
  if (titleNode) titleNode.textContent = "新型冠状病毒分型报告";
  if (subtitleNode) {
    subtitleNode.textContent = "围绕 SARS-CoV-2 Nextclade clade、Pango 谱系、质控指标、突变位点、系统发育树与参考比对结果组织单页展示。";
  }
  if (statusNode) {
    const rawStatus = String(task.status || "").trim();
    statusNode.textContent = rawStatus || "分析完成";
  }
  if (sampleTitleNode) sampleTitleNode.textContent = sampleName;
  if (sampleCopyNode) {
    sampleCopyNode.textContent = "该页面已收敛为 SARS-CoV-2 专用单页视图，重点展示 Nextclade 结果、变异位点、系统发育树和参考比对结果。";
  }
  if (typeof document !== "undefined") {
    document.title = `${sampleName} - 新型冠状病毒分型报告`;
  }
  buildSarsCov2ReportNav();
}

function buildSarsCov2ReportNav() {
  const nav = document.querySelector(".report-nav");
  if (!nav) return;
  const hasTree = String(currentReportData?.sections?.serotype?.phylogeny_tree?.status || "") === "ready";
  const hasVariantAnnotation = String(currentReportData?.sections?.serotype?.variant_annotation?.status || "") === "ready"
    && Array.isArray(currentReportData?.sections?.serotype?.variant_annotation?.rows)
    && currentReportData.sections.serotype.variant_annotation.rows.length > 0;
  const hasIgv = String(currentReportData?.sections?.serotype?.igv?.status || "") === "ready";
  const hasFastp = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "") === "ready";
  const hasCoverage = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
    && currentReportData.sections.assembly.coverage.points.length > 0;
  const groups = [
    {
      section: "section-raw-qc",
      title: "1. 质控与覆盖",
      id: "nav-group-ncov-qc",
      children: [
        { href: "#section-raw-qc", label: "1.1 原始数据质控" },
        ...(hasFastp ? [{ href: "#section-fastp", label: "1.2 fastp 结果可视化" }] : []),
        ...(hasCoverage ? [{ href: "#section-assembly-summary", label: hasFastp ? "1.3 测序深度" : "1.2 测序深度" }] : []),
      ],
    },
    {
      section: "section-serotype",
      title: "2. Nextclade 分型",
      id: "nav-group-ncov-serotype",
      children: [
        { href: "#nextclade-summary", label: "2.1 结果总表" },
        { href: "#nextclade-assignment", label: "2.2 分型摘要" },
        ...(hasVariantAnnotation ? [{ href: "#nextclade-variant-annotation", label: "2.3 变异注释表" }] : []),
        ...(hasIgv ? [{ href: "#nextclade-igv", label: hasVariantAnnotation ? "2.4 IGV 比对结果" : "2.3 IGV 比对结果" }] : []),
        ...(hasTree ? [{ href: "#nextclade-phylogeny", label: hasVariantAnnotation ? (hasIgv ? "2.5 系统发育树" : "2.4 系统发育树") : (hasIgv ? "2.4 系统发育树" : "2.3 系统发育树") }] : []),
      ],
    },
    {
      section: "nextclade-gene-workspace",
      title: "3. 突变分析",
      id: "nav-group-ncov-mutation",
      children: [
        { href: "#nextclade-gene-workspace", label: "3.1 按基因查看突变" },
        { href: "#nextclade-linked-track", label: "3.2 关联突变轨道" },
        { href: "#nextclade-linked-table", label: "3.3 关联突变表" },
      ],
    },
  ].filter((group) => group.children.length);
  nav.innerHTML = `
    ${groups.map((group) => `
      <div class="report-nav-group has-children" data-nav-group>
        <button class="report-nav-link report-nav-toggle" type="button" data-nav-toggle data-nav-section="${group.section}" aria-expanded="false" aria-controls="${group.id}">
          <span>${group.title}</span>
        </button>
        <div class="report-subnav" id="${group.id}" hidden>
          ${group.children.map((item) => `<a class="report-nav-link subnav-link" href="${item.href}">${item.label}</a>`).join("")}
        </div>
      </div>
    `).join("")}
  `;
}

function applyMonkeypoxReportChrome(data) {
  if (!isMonkeypoxNextcladeReport(data)) return;
  const task = data?.task || {};
  const backNode = document.querySelector(".report-back");
  const subtitleNode = document.querySelector(".report-subtitle");
  const kickerNode = document.querySelector(".report-title-block .report-kicker");
  const titleNode = document.querySelector(".report-title-block h1");
  const statusNode = document.querySelector(".report-status");
  const sampleTitleNode = document.getElementById("report-sample-title");
  const sampleCopyNode = document.getElementById("report-sample-copy");
  const sampleName = task.sample_display_name || task.sample_name || data?.sections?.serotype?.sequence_name || task.name || task.id || "hMPXV";
  if (backNode) backNode.textContent = "返回任务";
  if (kickerNode) kickerNode.textContent = "Monkeypox Report";
  if (titleNode) titleNode.textContent = "猴痘分型报告";
  if (subtitleNode) {
    subtitleNode.textContent = "围绕 hMPXV Nextclade 结果组织 clade、lineage、outbreak、质控指标和参考比对视图，适合做猴痘 demo 演示。";
  }
  if (statusNode) {
    const rawStatus = String(task.status || "").trim();
    statusNode.textContent = rawStatus || "分析完成";
  }
  if (sampleTitleNode) sampleTitleNode.textContent = sampleName;
  if (sampleCopyNode) {
    sampleCopyNode.textContent = "该页面已收敛为猴痘专用单页视图，重点展示 Nextclade 结论、质控指标以及参考比对结果。";
  }
  if (typeof document !== "undefined") {
    document.title = `${sampleName} - 猴痘分型报告`;
  }
  buildMonkeypoxReportNav();
}

function applyRsvReportChrome(data) {
  const isRsv = isRsvNextcladeReport(data);
  const isHmpv = isHmpvNextcladeReport(data);
  const isDenv = isDenvNextcladeReport(data);
  const isZikav = isZikavNextcladeReport(data);
  const isChikv = isChikvNextcladeReport(data);
  const isEbola = isEbolaNextcladeReport(data);
  const isHpiv = isHpivTypingReport(data);
  const isHadv = isHadvTypingReport(data);
  const isNorovirus = isNorovirusTypingReport(data);
  const isEnterovirus = isEnterovirusTypingReport(data);
  const isHiv = isHivTypingReport(data);
  const isHepatovirus = isHepatovirusTypingReport(data);
  const isBandavirus = isBandavirusTypingReport(data);
  const isOrthohantavirus = isOrthohantavirusTypingReport(data);
  const isAstroviridae = isAstroviridaeTypingReport(data);
  const isRhinovirus = isRhinovirusTypingReport(data);
  const isSeasonalHcov = isSeasonalHcovTypingReport(data);
  const isRotavirus = isRotavirusTypingReport(data);
  if (!isRsv && !isHmpv && !isDenv && !isZikav && !isChikv && !isEbola && !isHpiv && !isHadv && !isNorovirus && !isEnterovirus && !isHiv && !isHepatovirus && !isBandavirus && !isOrthohantavirus && !isAstroviridae && !isRhinovirus && !isSeasonalHcov && !isRotavirus) return;
  const task = data?.task || {};
  const backNode = document.querySelector(".report-back");
  const subtitleNode = document.querySelector(".report-subtitle");
  const kickerNode = document.querySelector(".report-title-block .report-kicker");
  const titleNode = document.querySelector(".report-title-block h1");
  const statusNode = document.querySelector(".report-status");
  const sampleTitleNode = document.getElementById("report-sample-title");
  const sampleCopyNode = document.getElementById("report-sample-copy");
  const hepatovirusBroad = String(currentReportData?.sections?.serotype?.predicted_group || "").trim().toUpperCase();
  const hepatovirusLabel = ({ HAV: "甲型肝炎病毒", HBV: "乙型肝炎病毒", HCV: "丙型肝炎病毒", HDV: "丁型肝炎病毒", HEV: "戊型肝炎病毒" })[hepatovirusBroad] || "肝炎病毒";
  const virusShort = isHiv ? "HIV" : (isRotavirus ? "RotaV" : (isNorovirus ? "NoV" : (isEnterovirus ? "EV" : (isHepatovirus ? (hepatovirusBroad || "HepV") : (isBandavirus ? "BandV" : (isOrthohantavirus ? "HTNV" : (isEbola ? "EBOV" : (isAstroviridae ? "AstV" : (isRhinovirus ? "HRV" : (isSeasonalHcov ? "HCoV" : (isChikv ? "CHIKV" : (isZikav ? "ZIKV" : (isDenv ? "DENV" : (isHmpv ? "HMPV" : (isHpiv ? "HPIV" : (isHadv ? "HAdV" : "RSV"))))))))))))))));
  const virusLabel = isHiv ? "HIV" : (isRotavirus ? "轮状病毒" : (isNorovirus ? "诺如病毒" : (isEnterovirus ? "肠道病毒" : (isHepatovirus ? hepatovirusLabel : (isBandavirus ? "班达病毒" : (isOrthohantavirus ? "汉坦病毒" : (isEbola ? "埃博拉病毒" : (isAstroviridae ? "星状病毒" : (isRhinovirus ? "鼻病毒" : (isSeasonalHcov ? "季节性冠状病毒" : (isChikv ? "基孔肯雅病毒" : (isZikav ? "寨卡病毒" : (isDenv ? "登革热病毒" : (isHmpv ? "人偏肺病毒" : (isHpiv ? "人副流感病毒" : (isHadv ? "人腺病毒" : "RSV"))))))))))))))));
  const sampleName = task.sample_display_name || task.sample_name || data?.sections?.serotype?.sequence_name || task.name || task.id || virusShort;
  if (backNode) backNode.textContent = "返回任务";
  if (kickerNode) kickerNode.textContent = `${virusShort} Report`;
  if (titleNode) titleNode.textContent = `${virusLabel}分型报告`;
  if (subtitleNode) {
    subtitleNode.textContent = isHpiv
      ? "围绕 HPIV1/2/3/4A/4B 自动选参考、对应注释文件、变异位点与参考比对结果组织单页展示。"
      : isNorovirus
      ? "围绕 Norovirus RdRp/VP1 双位点分型、最优参考选择、变异位点与参考比对结果组织单页展示。"
      : isEnterovirus
      ? "围绕 Enterovirus EV-A/B/C/D VP1 分型、95% 去冗余候选参考竞争、VADR 注释与参考比对结果组织单页展示。"
      : isHiv
      ? "围绕 HIV-1/HIV-2 broad 分型、HIV-1 子亚型代表株筛选、REGA-like 重组分析与 HIVDB 耐药解释组织单页展示。"
      : isHepatovirus
      ? "围绕肝炎病毒大亚型筛选、对应子亚型/基因型参考竞争、最优参考选择、突变位点与参考比对结果组织单页展示。"
      : isBandavirus
      ? "围绕 Bandavirus 大亚型筛选、SFTSV 的 A_F/CJ 三片段分型、最优参考选择、重组提示与参考比对结果组织单页展示。"
      : isOrthohantavirus
      ? "围绕 Orthohantavirus broad 分型、S 片段分型、L/M/S 三片段最优参考、知识库亚型关联与参考比对结果组织单页展示。"
      : isEbola
      ? "围绕 Orthoebolavirus 本地参考筛选、Ebola Nextclade clade / lineage、质控指标、突变位点、系统发育树与参考比对结果组织单页展示。"
      : isAstroviridae
      ? "围绕 Astrovirus ORF2 分型、Mamastrovirus/Avastrovirus 属级判定、最优参考选择、VADR 注释与 ORF2 系统树组织单页展示。"
      : isRhinovirus
      ? "围绕 Rhinovirus VP1 分型、物种组判定、最优参考选择、变异位点与参考比对结果组织单页展示。"
      : isSeasonalHcov
      ? "围绕季节性冠状病毒 229E/NL63/OC43/HKU1 参考筛选、VADR 注释、S 基因子亚型和系统发育树组织单页展示。"
      : isRotavirus
      ? "围绕轮状病毒 A/B/C 大组覆盖度比较、A 组 VP4/VP7 组合分型、最优参考株与组装覆盖度结果组织单页展示。"
      : isHadv
      ? "围绕 HAdV 的 PHF 三基因分型、总分型判定、最优参考选择、变异位点与参考比对结果组织单页展示。"
      : isChikv
      ? "围绕 CHIKV 固定参考、Nextclade 分型、质控指标、变异位点与参考比对结果组织单页展示。"
      : isZikav
      ? "围绕 ZIKV 固定参考、Nextclade 分型、质控指标、变异位点与参考比对结果组织单页展示。"
      : isDenv
      ? "围绕 DENV1-4 自动选参考、对应注释文件与对应亚型 Nextclade 数据集分型、质控指标和参考比对结果组织单页展示。"
      : isHmpv
      ? "围绕固定 HMPV 参考基因组、对应 Nextclade 数据集分型、质控指标与参考比对结果组织单页展示。"
      : "围绕 RSV A/B 自动选参考、对应 Nextclade 数据集分型、质控指标与参考比对结果组织单页展示。";
  }
  if (statusNode) {
    const rawStatus = String(task.status || "").trim();
    statusNode.textContent = rawStatus || "分析完成";
  }
  if (sampleTitleNode) sampleTitleNode.textContent = sampleName;
  if (sampleCopyNode) {
    sampleCopyNode.textContent = isHpiv
      ? "该页面已收敛为 HPIV 专用单页视图，重点展示自动选择的参考型别、变异位点和参考比对结果。"
      : isNorovirus
      ? "该页面已收敛为 Norovirus 专用单页视图，重点展示 RdRp/VP1 双位点分型、最优参考和参考比对结果。"
      : isEnterovirus
      ? "该页面已收敛为 Enterovirus 专用单页视图，重点展示 VP1 分型、大亚型、最优参考和参考比对结果。"
      : isHiv
      ? "该页面已收敛为 HIV 专用单页视图，重点展示 broad 分型、子亚型/重组、HIVDB 耐药解释和参考比对结果。"
      : isHepatovirus
      ? "该页面已收敛为肝炎病毒专用单页视图，重点展示大亚型、子亚型/基因型、最优参考和参考比对结果。"
      : isBandavirus
      ? "该页面已收敛为 Bandavirus 专用单页视图，重点展示大亚型、A_F/CJ 三片段分型、重组提示、最优参考和参考比对结果。"
      : isOrthohantavirus
      ? "该页面已收敛为汉坦病毒专用单页视图，重点展示 broad 分型、S 片段分型、L/M/S 三片段证据、知识库亚型关联和参考比对结果。"
      : isEbola
      ? "该页面已收敛为埃博拉病毒专用单页视图，重点展示本地参考分型、Ebola Nextclade 结果、突变位点和参考比对结果。"
      : isAstroviridae
      ? "该页面已收敛为 Astrovirus 专用单页视图，重点展示 ORF2 分型、病毒属/种、最优参考和参考比对结果。"
      : isRhinovirus
      ? "该页面已收敛为 Rhinovirus 专用单页视图，重点展示 VP1 分型、物种组、最优参考和参考比对结果。"
      : isSeasonalHcov
      ? "该页面已收敛为季节性冠状病毒专用单页视图，重点展示大类分型、S 子亚型、最优参考和参考比对结果。"
      : isRotavirus
      ? "该页面已收敛为轮状病毒专用单页视图，重点展示 A/B/C 大组覆盖度、G/P 组合分型、最优参考株和组装覆盖度结果。"
      : isHadv
      ? "该页面已收敛为 HAdV 专用单页视图，重点展示 PHF 三基因分型、总分型、最优参考和参考比对结果。"
      : isChikv
      ? "该页面已收敛为 CHIKV 专用单页视图，重点展示 Nextclade 结果、变异位点和参考比对结果。"
      : isZikav
      ? "该页面已收敛为 ZIKV 专用单页视图，重点展示 Nextclade 结果、变异位点和参考比对结果。"
      : isDenv
      ? "该页面已收敛为 DENV 专用单页视图，重点展示自动选择的参考型别、Nextclade 结果、变异位点和参考比对结果。"
      : isHmpv
      ? "该页面已收敛为 HMPV 专用单页视图，重点展示固定参考下的 Nextclade 结果、变异位点和参考比对结果。"
      : "该页面已收敛为 RSV 专用单页视图，重点展示自动选择的参考型别、Nextclade 结果、变异位点和参考比对结果。";
  }
  if (typeof document !== "undefined") {
    document.title = `${sampleName} - ${virusShort} 分型报告`;
  }
  buildRsvReportNav();
}

function buildRsvReportNav() {
  const nav = document.querySelector(".report-nav");
  if (!nav) return;
  const isHmpv = isHmpvNextcladeReport(currentReportData);
  const isDenv = isDenvNextcladeReport(currentReportData);
  const isZikav = isZikavNextcladeReport(currentReportData);
  const isChikv = isChikvNextcladeReport(currentReportData);
  const isEbola = isEbolaNextcladeReport(currentReportData);
  const isHpiv = isHpivTypingReport(currentReportData);
  const isHadv = isHadvTypingReport(currentReportData);
  const isNorovirus = isNorovirusTypingReport(currentReportData);
  const isEnterovirus = isEnterovirusTypingReport(currentReportData);
  const isHiv = isHivTypingReport(currentReportData);
  const isHepatovirus = isHepatovirusTypingReport(currentReportData);
  const isBandavirus = isBandavirusTypingReport(currentReportData);
  const isOrthohantavirus = isOrthohantavirusTypingReport(currentReportData);
  const isAstroviridae = isAstroviridaeTypingReport(currentReportData);
  const isRhinovirus = isRhinovirusTypingReport(currentReportData);
  const isSeasonalHcov = isSeasonalHcovTypingReport(currentReportData);
  const isRotavirus = isRotavirusTypingReport(currentReportData);
  const hepatovirusBroad = String(currentReportData?.sections?.serotype?.predicted_group || "").trim().toUpperCase();
  const virusLabel = isHiv ? "HIV" : (isRotavirus ? "RotaV" : (isNorovirus ? "NoV" : (isEnterovirus ? "EV" : (isHepatovirus ? (hepatovirusBroad || "HepV") : (isBandavirus ? "BandV" : (isOrthohantavirus ? "HTNV" : (isEbola ? "EBOV" : (isAstroviridae ? "AstV" : (isRhinovirus ? "HRV" : (isSeasonalHcov ? "HCoV" : (isChikv ? "CHIKV" : (isZikav ? "ZIKV" : (isDenv ? "DENV" : (isHmpv ? "HMPV" : (isHpiv ? "HPIV" : (isHadv ? "HAdV" : "RSV"))))))))))))))));
  const hasFastp = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "") === "ready";
  const hasSpecies = Array.isArray(currentReportData?.sections?.species_identification?.species?.rows)
    && currentReportData.sections.species_identification.species.rows.length > 0;
  const hasSubspecies = Array.isArray(currentReportData?.sections?.species_identification?.subspecies?.rows)
    && currentReportData.sections.species_identification.subspecies.rows.length > 0;
  const hasCoverage = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
    && currentReportData.sections.assembly.coverage.points.length > 0;
  const hasAssemblySummary = Array.isArray(currentReportData?.sections?.assembly?.summary?.rows)
    && currentReportData.sections.assembly.summary.rows.length > 0;
  const hasQuality = Array.isArray(currentReportData?.sections?.serotype?.quality_metrics)
    && currentReportData.sections.serotype.quality_metrics.length > 0;
  const hasPhf = isHadv && Array.isArray(currentReportData?.sections?.serotype?.phf_table?.rows)
    && currentReportData.sections.serotype.phf_table.rows.length > 0;
  const hasPhfCoverage = isHadv && Array.isArray(currentReportData?.sections?.serotype?.phf_coverage?.rows)
    && currentReportData.sections.serotype.phf_coverage.rows.length > 0;
  const hasPhfSnp = isHadv && Array.isArray(currentReportData?.sections?.serotype?.phf_snp?.rows)
    && currentReportData.sections.serotype.phf_snp.rows.length > 0;
  const hasHivBroad = isHiv && Array.isArray(currentReportData?.sections?.serotype?.broad_typing?.rows)
    && currentReportData.sections.serotype.broad_typing.rows.length > 0;
  const hasHivReference = isHiv && Array.isArray(currentReportData?.sections?.serotype?.subtype_reference_typing?.rows)
    && currentReportData.sections.serotype.subtype_reference_typing.rows.length > 0;
  const hasHivSelection = isHiv && Array.isArray(currentReportData?.sections?.serotype?.reference_selection?.rows)
    && currentReportData.sections.serotype.reference_selection.rows.length > 0;
  const hasHivResistance = isHiv && Array.isArray(currentReportData?.sections?.serotype?.resistance_table?.rows)
    && currentReportData.sections.serotype.resistance_table.rows.length > 0;
  const hasHivMutationEvidence = isHiv && Array.isArray(currentReportData?.sections?.serotype?.mutation_panels)
    && currentReportData.sections.serotype.mutation_panels.length > 0;
  const hasHivBootscan = isHiv && (
    String(currentReportData?.sections?.serotype?.bootscan_assets?.overall_svg || "").trim()
    || String(currentReportData?.sections?.serotype?.bootscan_assets?.pure_svg || "").trim()
  );
  const hasBandavirusSelection = isBandavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.bandavirus_selection?.rows)
    && currentReportData.sections.serotype.bandavirus_selection.rows.length > 0;
  const hasBandavirusAf = isBandavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.af_segment_typing?.rows)
    && currentReportData.sections.serotype.af_segment_typing.rows.length > 0;
  const hasBandavirusCj = isBandavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.cj_segment_typing?.rows)
    && currentReportData.sections.serotype.cj_segment_typing.rows.length > 0;
  const hasBandavirusConsensus = isBandavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.consensus_typing?.rows)
    && currentReportData.sections.serotype.consensus_typing.rows.length > 0;
  const hasOrthohantavirusSelection = isOrthohantavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.orthohantavirus_selection?.rows)
    && currentReportData.sections.serotype.orthohantavirus_selection.rows.length > 0;
  const hasEbolaSelection = isEbolaNextcladeReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.orthoebolavirus_selection?.rows)
    && currentReportData.sections.serotype.orthoebolavirus_selection.rows.length > 0;
  const hasOrthohantavirusBroad = isOrthohantavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.broad_typing?.rows)
    && currentReportData.sections.serotype.broad_typing.rows.length > 0;
  const hasOrthohantavirusSegments = isOrthohantavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.segment_typing?.rows)
    && currentReportData.sections.serotype.segment_typing.rows.length > 0;
  const hasOrthohantavirusConsensus = isOrthohantavirusTypingReport(currentReportData) && Array.isArray(currentReportData?.sections?.serotype?.consensus_typing?.rows)
    && currentReportData.sections.serotype.consensus_typing.rows.length > 0;
  const hasMutations = Array.isArray(currentReportData?.sections?.serotype?.mutation_table?.rows)
    && currentReportData.sections.serotype.mutation_table.rows.length > 0;
  const hasIgv = String(currentReportData?.sections?.serotype?.igv?.status || "") === "ready";
  const hasFunctionalAnnotation = isHpiv && Array.isArray(currentReportData?.sections?.serotype?.functional_annotation?.rows)
    && currentReportData.sections.serotype.functional_annotation.rows.length > 0;
  const hasNmdcAnnotation = !isHmpv && !isDenv && !isZikav && !isChikv && !isEbola && !isHpiv && !isNorovirus && !isEnterovirus && !isHepatovirus && !isBandavirus && !isRhinovirus && !isSeasonalHcov && Array.isArray(currentReportData?.sections?.serotype?.nmdc_annotation?.rows)
    && currentReportData.sections.serotype.nmdc_annotation.rows.length > 0;
  const hasPhylogeny = (isNorovirus || isEnterovirus || isAstroviridae || isRhinovirus || isSeasonalHcov)
    ? (Array.isArray(currentReportData?.sections?.serotype?.gene_phylogeny?.trees)
      && currentReportData.sections.serotype.gene_phylogeny.trees.some((item) => String(item?.status || "") === "ready"))
    : String(currentReportData?.sections?.serotype?.phylogeny_tree?.status || "") === "ready";
  const typingChildren = [{ href: "#rsv-typing-summary", label: "3.1 分型总表" }];
  let nextTypingIndex = 2;
  if (hasPhf) {
    typingChildren.push({ href: "#rsv-typing-phf", label: `3.${nextTypingIndex} PHF 三基因分型` });
    nextTypingIndex += 1;
  }
  if (hasPhfCoverage) {
    typingChildren.push({ href: "#rsv-typing-phf-coverage", label: `3.${nextTypingIndex} PHF 三基因覆盖度图` });
    nextTypingIndex += 1;
  }
  if (hasPhfSnp) {
    typingChildren.push({ href: "#rsv-typing-phf-snp", label: `3.${nextTypingIndex} PHF 三基因差异SNP数` });
    nextTypingIndex += 1;
  }
  if (hasHivBroad) {
    typingChildren.push({ href: "#rsv-typing-hiv-broad", label: `3.${nextTypingIndex} HIV-1/HIV-2 broad 分型` });
    nextTypingIndex += 1;
  }
  if (hasHivReference) {
    typingChildren.push({ href: "#rsv-typing-hiv-reference", label: `3.${nextTypingIndex} 子亚型代表株筛选` });
    nextTypingIndex += 1;
  }
  if (hasHivSelection) {
    typingChildren.push({ href: "#rsv-typing-hiv-selection", label: `3.${nextTypingIndex} 参考选择总表` });
    nextTypingIndex += 1;
  }
  if (hasQuality) {
    typingChildren.push({ href: "#rsv-typing-quality", label: `3.${nextTypingIndex} 质量指标` });
    nextTypingIndex += 1;
  }
  if (hasHivBootscan) {
    typingChildren.push({ href: "#rsv-typing-hiv-subtyping", label: `3.${nextTypingIndex} 子亚型 / 重组摘要` });
    nextTypingIndex += 1;
  }
  if (hasHivResistance) {
    typingChildren.push({ href: "#rsv-typing-hiv-resistance", label: `3.${nextTypingIndex} HIVDB 耐药解释` });
    nextTypingIndex += 1;
  }
  if (hasHivMutationEvidence) {
    typingChildren.push({ href: "#rsv-typing-hiv-mutations", label: `3.${nextTypingIndex} PR/RT/IN 线索` });
    nextTypingIndex += 1;
  }
  if (hasBandavirusSelection) {
    typingChildren.push({ href: "#rsv-typing-bandavirus-selection", label: `3.${nextTypingIndex} 参考筛选结果` });
    nextTypingIndex += 1;
  }
  if (hasBandavirusAf) {
    typingChildren.push({ href: "#rsv-typing-bandavirus-af", label: `3.${nextTypingIndex} A_F 三片段证据` });
    nextTypingIndex += 1;
  }
  if (hasBandavirusCj) {
    typingChildren.push({ href: "#rsv-typing-bandavirus-cj", label: `3.${nextTypingIndex} CJ 三片段证据` });
    nextTypingIndex += 1;
  }
  if (hasBandavirusConsensus) {
    typingChildren.push({ href: "#rsv-typing-bandavirus-consensus", label: `3.${nextTypingIndex} Consensus 复核` });
    nextTypingIndex += 1;
  }
  if (hasOrthohantavirusSelection) {
    typingChildren.push({ href: "#rsv-typing-orthohantavirus-selection", label: `3.${nextTypingIndex} 参考筛选结果` });
    nextTypingIndex += 1;
  }
  if (hasEbolaSelection) {
    typingChildren.push({ href: "#rsv-typing-ebola-selection", label: `3.${nextTypingIndex} 参考筛选结果` });
    nextTypingIndex += 1;
  }
  if (hasOrthohantavirusBroad) {
    typingChildren.push({ href: "#rsv-typing-orthohantavirus-broad", label: `3.${nextTypingIndex} Broad 分型支持` });
    nextTypingIndex += 1;
  }
  if (hasOrthohantavirusSegments) {
    typingChildren.push({ href: "#rsv-typing-orthohantavirus-segments", label: `3.${nextTypingIndex} L/M/S 三片段证据` });
    nextTypingIndex += 1;
  }
  if (hasOrthohantavirusConsensus) {
    typingChildren.push({ href: "#rsv-typing-orthohantavirus-consensus", label: `3.${nextTypingIndex} Consensus 复核` });
    nextTypingIndex += 1;
  }
  if (hasMutations) {
    typingChildren.push({ href: "#rsv-typing-mutations", label: `3.${nextTypingIndex} 突变位点表` });
    nextTypingIndex += 1;
  }
  if (hasIgv) {
    typingChildren.push({ href: "#rsv-typing-igv", label: `3.${nextTypingIndex} IGV 比对结果` });
    nextTypingIndex += 1;
  }
  if (hasFunctionalAnnotation) {
    typingChildren.push({ href: "#rsv-typing-functional-impact", label: `3.${nextTypingIndex} 突变功能影响` });
    nextTypingIndex += 1;
  }
  if (hasNmdcAnnotation) {
    typingChildren.push({ href: "#rsv-typing-nmdc", label: `3.${nextTypingIndex} 数据库注释` });
    nextTypingIndex += 1;
  }
  if (hasPhylogeny) {
    typingChildren.push({ href: "#rsv-phylogeny", label: `3.${nextTypingIndex} 系统发育树` });
  }
  const groups = [
    {
      section: "section-raw-qc",
      title: "1. 质控与物种鉴定",
      id: "nav-group-rsv-qc",
      children: [
        { href: "#section-raw-qc", label: "1.1 原始数据质控" },
        ...(hasFastp ? [{ href: "#section-fastp", label: "1.2 fastp 结果可视化" }] : []),
        ...((hasSpecies || hasSubspecies) ? [{ href: "#section-species-identification", label: hasFastp ? "1.3 物种鉴定" : "1.2 物种鉴定" }] : []),
      ],
    },
    {
      section: "section-assembly",
      title: "2. 组装与覆盖",
      id: "nav-group-rsv-assembly",
      children: [
        ...(hasAssemblySummary ? [{ href: "#section-assembly-summary", label: "2.1 组装后信息统计" }] : []),
        ...(hasCoverage ? [{ href: "#section-assembly", label: hasAssemblySummary ? "2.2 测序深度" : "2.1 测序深度" }] : []),
      ],
    },
    {
      section: "section-serotype",
      title: `3. ${virusLabel} 分型`,
      id: "nav-group-rsv-typing",
      children: typingChildren,
    },
  ].filter((group) => group.children.length);
  nav.innerHTML = `
    ${groups.map((group) => `
      <div class="report-nav-group has-children" data-nav-group>
        <button class="report-nav-link report-nav-toggle" type="button" data-nav-toggle data-nav-section="${group.section}" aria-expanded="false" aria-controls="${group.id}">
          <span>${escapeHtml(group.section === "section-serotype" ? "🦠 " : "")}${escapeHtml(group.title)}</span>
        </button>
        <div id="${group.id}" class="report-subnav" hidden>
          ${group.children.map((child) => `<a class="report-nav-link subnav-link" href="${child.href}">${escapeHtml(child.label)}</a>`).join("")}
        </div>
      </div>
    `).join("")}
  `;
}

function buildMonkeypoxReportNav() {
  const nav = document.querySelector(".report-nav");
  if (!nav) return;
  const hasFastp = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "") === "ready";
  const hasSpecies = Array.isArray(currentReportData?.sections?.species_identification?.species?.rows)
    && currentReportData.sections.species_identification.species.rows.length > 0;
  const hasSubspecies = Array.isArray(currentReportData?.sections?.species_identification?.subspecies?.rows)
    && currentReportData.sections.species_identification.subspecies.rows.length > 0;
  const hasCoverage = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
    && currentReportData.sections.assembly.coverage.points.length > 0;
  const hasAssemblySummary = Array.isArray(currentReportData?.sections?.assembly?.summary?.rows)
    && currentReportData.sections.assembly.summary.rows.length > 0;
  const hasQuality = Array.isArray(currentReportData?.sections?.serotype?.quality_metrics)
    && currentReportData.sections.serotype.quality_metrics.length > 0;
  const hasMutations = Array.isArray(currentReportData?.sections?.serotype?.mutation_table?.rows)
    && currentReportData.sections.serotype.mutation_table.rows.length > 0;
  const hasIgv = String(currentReportData?.sections?.serotype?.igv?.status || "") === "ready";
  const groups = [
    {
      section: "section-raw-qc",
      title: "1. 质控与物种鉴定",
      id: "nav-group-hmpxv-qc",
      children: [
        { href: "#section-raw-qc", label: "1.1 原始数据质控" },
        ...(hasFastp ? [{ href: "#section-fastp", label: "1.2 fastp 结果可视化" }] : []),
        ...((hasSpecies || hasSubspecies) ? [{ href: "#section-species-identification", label: hasFastp ? "1.3 物种鉴定" : "1.2 物种鉴定" }] : []),
      ],
    },
    {
      section: "section-assembly",
      title: "2. 组装与覆盖",
      id: "nav-group-hmpxv-assembly",
      children: [
        ...(hasAssemblySummary ? [{ href: "#section-assembly-summary", label: "2.1 组装后信息统计" }] : []),
        ...(hasCoverage ? [{ href: "#section-assembly", label: hasAssemblySummary ? "2.2 测序深度" : "2.1 测序深度" }] : []),
      ],
    },
    {
      section: "section-serotype",
      title: "3. 猴痘分型",
      id: "nav-group-hmpxv-typing",
      children: [
        { href: "#monkeypox-typing-summary", label: "3.1 分型总表" },
        ...(hasQuality ? [{ href: "#monkeypox-typing-quality", label: "3.2 质量指标" }] : []),
        ...(hasMutations ? [{ href: "#monkeypox-typing-mutations", label: hasQuality ? "3.3 突变位点表" : "3.2 突变位点表" }] : []),
        ...(hasIgv ? [{ href: "#monkeypox-typing-igv", label: hasQuality ? (hasMutations ? "3.4 IGV 比对结果" : "3.3 IGV 比对结果") : (hasMutations ? "3.3 IGV 比对结果" : "3.2 IGV 比对结果") }] : []),
      ],
    },
  ].filter((group) => group.children.length);
  nav.innerHTML = `
    ${groups.map((group) => `
      <div class="report-nav-group has-children" data-nav-group>
        <button class="report-nav-link report-nav-toggle" type="button" data-nav-toggle data-nav-section="${group.section}" aria-expanded="false" aria-controls="${group.id}">
          <span>${escapeHtml(group.title)}</span>
        </button>
        <div id="${group.id}" class="report-subnav" hidden>
          ${group.children.map((child) => `<a class="report-nav-link subnav-link" href="${child.href}">${escapeHtml(child.label)}</a>`).join("")}
        </div>
      </div>
    `).join("")}
  `;
}

function applyInfluenzaReportChrome(data) {
  if (!isInfluenzaTypingReport(data)) return;
  const task = data?.task || {};
  const backNode = document.querySelector(".report-back");
  const subtitleNode = document.querySelector(".report-subtitle");
  const kickerNode = document.querySelector(".report-title-block .report-kicker");
  const titleNode = document.querySelector(".report-title-block h1");
  const statusNode = document.querySelector(".report-status");
  const sampleTitleNode = document.getElementById("report-sample-title");
  const sampleCopyNode = document.getElementById("report-sample-copy");
  const sampleName = task.sample_display_name || task.sample_name || task.name || task.id || "Influenza";
  if (backNode) backNode.textContent = "返回任务";
  if (kickerNode) kickerNode.textContent = "Influenza Report";
  if (titleNode) titleNode.textContent = "流感分型报告";
  if (subtitleNode) {
    subtitleNode.textContent = "围绕流感 A/B 判断、HA/NA 亚型选择、8 个 segment 参考组成和有参覆盖情况组织结果，更接近病毒分型报告而不是单菌流程报告。";
  }
  if (statusNode) {
    const rawStatus = String(task.status || "").trim();
    statusNode.textContent = rawStatus || "分析完成";
  }
  if (sampleTitleNode) sampleTitleNode.textContent = sampleName;
  if (sampleCopyNode) {
    sampleCopyNode.textContent = "该页面已收敛为流感专用单页视图，重点展示质控、物种鉴定、流感分型结论和 8 segment 参考集合。";
  }
  const assemblySection = document.getElementById("section-assembly");
  const assemblySummarySection = document.getElementById("section-assembly-summary");
  const serotypeSection = document.getElementById("section-serotype");
  const assemblyKicker = assemblySection?.querySelector(".section-heading .report-kicker");
  const assemblyTitle = assemblySection?.querySelector(".section-heading h2");
  const assemblyCopy = assemblySection?.querySelector(".section-heading p:last-child");
  const assemblySummaryKicker = assemblySummarySection?.querySelector(".section-heading .report-kicker");
  const assemblySummaryTitle = assemblySummarySection?.querySelector(".section-heading h2");
  const serotypeKicker = serotypeSection?.querySelector(".section-heading .report-kicker");
  const serotypeTitle = serotypeSection?.querySelector(".section-heading h2");
  const serotypeCopy = serotypeSection?.querySelector(".section-heading p:last-child");
  const coverageCardTitle = document.querySelector("#assembly-coverage-card .card-head h3");
  const coverageCardCopy = document.querySelector("#assembly-coverage-chart p");
  if (assemblyKicker) assemblyKicker.textContent = "Section 2";
  if (assemblyTitle) assemblyTitle.textContent = "组装与覆盖";
  if (assemblyCopy) assemblyCopy.textContent = "围绕最终 8 个 segment 的组装统计和测序深度进行整理，用于支持流感分型和变异判读。";
  if (assemblySummaryKicker) assemblySummaryKicker.textContent = "Section 2.1";
  if (assemblySummaryTitle) assemblySummaryTitle.textContent = "组装后信息统计";
  if (serotypeKicker) serotypeKicker.textContent = "Section 3";
  if (serotypeTitle) serotypeTitle.textContent = "流感分型";
  if (serotypeCopy) serotypeCopy.textContent = "基于 wf_flu 风格的 IRMA reference set 初筛与 HA/NA 最优亚型选择，展示甲/乙流判断、亚型组合和最终 8 segment 参考组成。";
  if (coverageCardTitle) coverageCardTitle.textContent = "测序深度";
  if (coverageCardCopy) coverageCardCopy.textContent = "按流感 segment 展示测序深度与覆盖情况。";
  if (typeof document !== "undefined") {
    document.title = `${sampleName} - 流感分型报告`;
  }
  buildInfluenzaReportNav();
}

function buildInfluenzaReportNav() {
  const nav = document.querySelector(".report-nav");
  if (!nav) return;
  const hasFastp = String(currentReportData?.sections?.raw_qc?.fastp?.status || currentReportData?.sections?.raw_qc?.status || "") === "ready";
  const hasSpecies = Array.isArray(currentReportData?.sections?.species_identification?.species?.rows)
    && currentReportData.sections.species_identification.species.rows.length > 0;
  const hasSubspecies = Array.isArray(currentReportData?.sections?.species_identification?.subspecies?.rows)
    && currentReportData.sections.species_identification.subspecies.rows.length > 0;
  const hasCoverage = Array.isArray(currentReportData?.sections?.assembly?.coverage?.points)
    && currentReportData.sections.assembly.coverage.points.length > 0;
  const hasAssemblySummary = Array.isArray(currentReportData?.sections?.assembly?.summary?.rows)
    && currentReportData.sections.assembly.summary.rows.length > 0;
  const hasMutations = Array.isArray(currentReportData?.sections?.serotype?.mutation_table?.rows)
    && currentReportData.sections.serotype.mutation_table.rows.length > 0;
  const hasIgv = String(currentReportData?.sections?.serotype?.igv?.status || "") === "ready";
  const hasResistance = String(currentReportData?.sections?.serotype?.resistance_annotation?.status || "") === "ready";
  const groups = [
    {
      section: "section-raw-qc",
      title: "1. 质控与物种鉴定",
      id: "nav-group-flu-qc",
      children: [
        { href: "#section-raw-qc", label: "1.1 原始数据质控" },
        ...(hasFastp ? [{ href: "#section-fastp", label: "1.2 fastp 结果可视化" }] : []),
        ...((hasSpecies || hasSubspecies) ? [{ href: "#section-species-identification", label: hasFastp ? "1.3 物种鉴定" : "1.2 物种鉴定" }] : []),
      ],
    },
    {
      section: "section-assembly",
      title: "2. 组装与覆盖",
      id: "nav-group-flu-assembly",
      children: [
        ...(hasAssemblySummary ? [{ href: "#section-assembly-summary", label: "2.1 组装后信息统计" }] : []),
        ...(hasCoverage ? [{ href: "#section-assembly", label: hasAssemblySummary ? "2.2 测序深度" : "2.1 测序深度" }] : []),
      ],
    },
    {
      section: "section-serotype",
      title: "3. 流感分型",
      id: "nav-group-flu-typing",
      children: [
        { href: "#influenza-typing-summary", label: "3.1 分型总表" },
        { href: "#influenza-typing-manifest", label: "3.2 8 Segment 参考组成" },
        ...(hasMutations ? [{ href: "#influenza-typing-mutations", label: "3.3 变异注释表" }] : []),
        ...(hasIgv ? [{ href: "#influenza-typing-igv", label: "3.4 IGV 比对结果" }] : []),
        ...(hasResistance ? [{ href: "#influenza-typing-resistance", label: "3.5 耐药突变注释结果" }] : []),
      ],
    },
  ].filter((group) => group.children.length);
  nav.innerHTML = `
    ${groups.map((group) => `
      <div class="report-nav-group has-children" data-nav-group>
        <button class="report-nav-link report-nav-toggle" type="button" data-nav-toggle data-nav-section="${group.section}" aria-expanded="false" aria-controls="${group.id}">
          <span>${escapeHtml(group.title)}</span>
        </button>
        <div id="${group.id}" class="report-subnav" hidden>
          ${group.children.map((child) => `<a class="report-nav-link subnav-link" href="${child.href}">${escapeHtml(child.label)}</a>`).join("")}
        </div>
      </div>
    `).join("")}
  `;
}

function renderBinningSection(task, section) {
  const isMeta = getTaskMethod(task) === "meta";
  const wrapper = document.getElementById("nav-group-binning-wrapper");
  const mainSection = document.getElementById("section-binning");
  const qualitySection = document.getElementById("section-binning-quality");
  const taxonomySection = document.getElementById("section-binning-taxonomy");
  const viralSection = document.getElementById("section-binning-viral");
  [wrapper, mainSection, qualitySection, taxonomySection, viralSection].forEach((node) => {
    if (!node) return;
    node.classList.toggle("hidden", !isMeta);
  });
  if (!isMeta) return;
  const qualitySummary = section?.quality?.summary || {};
  const taxonomySummary = section?.taxonomy?.summary || {};
  const viralSummary = section?.viral_assembly?.summary || {};
  renderBinningMetricCards("binning-quality-summary", [
    { title: "bin 总数", value: qualitySummary.total_bins ?? "--", label: "已纳入完整性评估的 bin 数量", tag: "总览" },
    { title: "平均完整性", value: qualitySummary.avg_completeness != null ? `${qualitySummary.avg_completeness}%` : "--", label: "全部 bin 平均完整性", tag: "质量" },
    { title: "平均污染率", value: qualitySummary.avg_contamination != null ? `${qualitySummary.avg_contamination}%` : "--", label: "全部 bin 平均污染率", tag: "质量" },
    { title: "高/中/低质量", value: `${qualitySummary.hq_bins ?? 0}/${qualitySummary.mq_bins ?? 0}/${qualitySummary.lq_bins ?? 0}`, label: "高质量 / 中质量 / 低质量 bin", tag: "分层" },
  ]);
  renderBinningChartCards("binning-quality-charts", [
    section?.quality?.charts?.completeness || {},
    section?.quality?.charts?.contamination || {},
    section?.quality?.charts?.quality_tier || {},
  ]);
  renderBinningMetricCards("binning-taxonomy-summary", [
    { title: "已分类 bin", value: taxonomySummary.classified_bins ?? "--", label: "获得 GTDB-Tk 分类结果的 bin 数量", tag: "分类" },
    { title: "未分类 bin", value: taxonomySummary.unclassified_bins ?? "--", label: "仍为 Unclassified 的 bin 数量", tag: "分类" },
    { title: "优势门", value: taxonomySummary.top_phylum || "--", label: "bin 数量最多的门水平类群", tag: "门水平" },
    { title: "优势属", value: taxonomySummary.top_genus || "--", label: "bin 数量最多的属水平类群", tag: "属水平" },
  ]);
  renderBinningChartCards("binning-taxonomy-charts", [
    section?.taxonomy?.charts?.phylum || {},
    section?.taxonomy?.charts?.genus || {},
    section?.taxonomy?.charts?.method || {},
  ]);
  buildTableCard("binning-quality-table", "bin 完整性分析", section?.quality?.table?.columns || [], section?.quality?.table?.rows || []);
  buildTableCard("binning-taxonomy-table", "bin 物种鉴定结果", section?.taxonomy?.table?.columns || [], section?.taxonomy?.table?.rows || []);
  renderBinningMetricCards("binning-viral-summary", [
    { title: "候选 contig", value: viralSummary.candidate_contigs ?? "--", label: "进入病毒筛选链的组装 contig 数量", tag: "候选" },
    { title: "最终保留", value: viralSummary.retained_contigs ?? "--", label: "通过 VirSorter2 / CheckV 过滤后的 contig 数量", tag: "保留" },
    { title: "保留总长度", value: viralSummary.retained_length ? formatBases(viralSummary.retained_length) : "--", label: "最终病毒 contig 总长度", tag: "长度" },
    { title: "最佳质量", value: viralSummary.best_quality || "--", label: "CheckV 最高质量等级", tag: "质量" },
  ]);
  buildTableCard("binning-viral-table", "病毒组装 contig 筛选明细", section?.viral_assembly?.table?.columns || [], section?.viral_assembly?.table?.rows || []);
}

function renderSampleSwitcher(task) {
  const container = document.getElementById("report-sample-switcher");
  if (!container) return;
  const samples = Array.isArray(task?.samples) ? task.samples.filter(Boolean) : [];
  if (samples.length <= 1) {
    container.classList.add("hidden");
    container.innerHTML = "";
    return;
  }
  const explicitSample = new URLSearchParams(window.location.search).get("sample") || "";
  const isBatchLanding = String(task?.report_mode || "").trim() === "multi" && !explicitSample;
  const currentSample = explicitSample || task.sample_display_name || task.sample_name || "";
  const sampleUrl = (sample) => `${window.location.pathname}?sample=${encodeURIComponent(sample)}`;
  const batchUrl = () => window.location.pathname;
  const buildSampleOptions = () => samples.map((sample) => `
    <option value="${escapeHtml(sample)}"${sample === currentSample ? " selected" : ""}>${escapeHtml(sample)}</option>
  `).join("");
  const getVisibleSamples = (query = "") => {
    const keyword = String(query || "").trim().toLowerCase();
    if (keyword) {
      return samples.filter((sample) => sample.toLowerCase().includes(keyword));
    }
    if (samples.length <= 18) return samples;
    const visible = [];
    if (currentSample) visible.push(currentSample);
    samples.forEach((sample) => {
      if (visible.length >= 18) return;
      if (!visible.includes(sample)) visible.push(sample);
    });
    return visible;
  };
  const buildSampleLinks = (query = "") => {
    const visibleSamples = getVisibleSamples(query);
    const shownSamples = visibleSamples.slice(0, 40);
    if (shownSamples.length === 0) {
      return `<span class="report-sample-empty">未找到匹配样本</span>`;
    }
    const hiddenCount = visibleSamples.length - shownSamples.length;
    return `
      <a class="report-sample-chip report-sample-chip-overview${isBatchLanding ? " is-active" : ""}" href="${batchUrl()}">批次概览</a>
      ${shownSamples.map((sample) => {
        const active = !isBatchLanding && sample === currentSample;
        return `<a class="report-sample-chip${active ? " is-active" : ""}" href="${sampleUrl(sample)}">${escapeHtml(sample)}</a>`;
      }).join("")}
      ${hiddenCount > 0 ? `<span class="report-sample-overflow-note">还有 ${hiddenCount} 个</span>` : ""}
    `;
  };
  container.classList.remove("hidden");
  container.innerHTML = `
    <div class="report-sample-switcher-head">
      <span class="report-sample-switcher-label">样本列表</span>
      <span class="report-sample-switcher-meta">当前共 ${samples.length} 个样本</span>
    </div>
    <div class="report-sample-picker-row">
      <label class="report-sample-picker-field">
        <span>当前样本</span>
        <select class="report-sample-select" data-report-sample-select>
          <option value=""${isBatchLanding ? " selected" : ""}>批次概览</option>
          ${buildSampleOptions()}
        </select>
      </label>
      <label class="report-sample-picker-field">
        <span>筛选样本</span>
        <input class="report-sample-search" type="search" data-report-sample-search placeholder="输入样本名" autocomplete="off">
      </label>
    </div>
    <div class="report-sample-chip-row" data-report-sample-results>
      ${buildSampleLinks()}
    </div>
  `;
  const sampleSelect = container.querySelector("[data-report-sample-select]");
  const sampleSearch = container.querySelector("[data-report-sample-search]");
  const sampleResults = container.querySelector("[data-report-sample-results]");
  if (sampleSelect) {
    sampleSelect.addEventListener("change", () => {
      const nextSample = sampleSelect.value;
      if (!nextSample) {
        window.location.href = batchUrl();
        return;
      }
      if (nextSample && nextSample !== currentSample) {
        window.location.href = sampleUrl(nextSample);
      }
    });
  }
  if (sampleSearch && sampleResults) {
    sampleSearch.addEventListener("input", () => {
      sampleResults.innerHTML = buildSampleLinks(sampleSearch.value);
    });
  }
}

function isMultiSampleLanding(data) {
  const task = data?.task || {};
  const samples = Array.isArray(task.samples) ? task.samples : [];
  const explicitSample = new URLSearchParams(window.location.search).get("sample") || "";
  return String(task.report_mode || "").trim() === "multi" && samples.length > 1 && !explicitSample;
}

function renderMultiSampleOverview(data) {
  const container = document.getElementById("multi-sample-overview");
  const metrics = document.getElementById("overview-metrics");
  if (!container) return false;
  const summary = data?.task?.multi_sample_summary || data?.sections?.overview?.multi_sample_summary || {};
  const table = summary?.table || {};
  const columns = Array.isArray(table.columns) ? table.columns : [];
  const rows = Array.isArray(table.rows) ? table.rows : [];
  if (!rows.length || !columns.length) {
    container.classList.add("hidden");
    container.innerHTML = "";
    if (metrics) metrics.classList.remove("hidden");
    return false;
  }
  const analysisLabel = summary.analysis_target === "virus" ? "病毒任务" : "细菌任务";
  const readyCount = Number(summary.ready_count || 0);
  const sampleCount = Number(summary.sample_count || rows.length || 0);
  const speciesRank = Array.isArray(summary.species_rank) ? summary.species_rank : [];
  const topSpecies = speciesRank[0]?.name ? `${speciesRank[0].name} (${speciesRank[0].count || 0})` : "尚未形成";
  const exportColumns = columns.map((column) => column.label || column.key || "");
  const exportRows = rows.map((row) => columns.map((column) => {
    const key = String(column.key || "");
    const value = row?.[key];
    return value == null || value === "" ? "-" : String(value);
  }));
  const exportTitle = `${data?.task?.name || data?.task?.id || "多样本任务"}_${analysisLabel}多样本概览`;
  const rowClass = (row) => {
    const completeness = Number(row?.completeness);
    const contamination = Number(row?.contamination);
    const qcStatus = String(row?.qc_status || "").trim().toLowerCase();
    if (!row?.ready) return "is-missing";
    if ((Number.isFinite(completeness) && completeness < 90) || (Number.isFinite(contamination) && contamination >= 5) || ["bad", "mediocre", "failed", "fail"].includes(qcStatus)) {
      return "needs-attention";
    }
    return "";
  };
  const renderCell = (row, column) => {
    const key = String(column.key || "");
    const value = row?.[key];
    if (key === "sample") {
      return `<a class="multi-sample-link" href="${window.location.pathname}?sample=${encodeURIComponent(String(row.sample || ""))}">${escapeHtml(row.sample || "-")}</a>`;
    }
    if (key === "note") {
      return `<span class="multi-sample-note" title="${escapeHtml(value || "")}">${escapeHtml(value || "-")}</span>`;
    }
    if (key === "ready_label") {
      return `<span class="multi-sample-status ${row?.ready ? "is-ready" : "is-missing"}">${escapeHtml(value || "-")}</span>`;
    }
    return escapeHtml(value == null || value === "" ? "-" : value);
  };
  const sampleUrl = (row) => `${window.location.pathname}?sample=${encodeURIComponent(String(row?.sample || ""))}`;
  container.classList.remove("hidden");
  if (metrics) metrics.classList.add("hidden");
  container.innerHTML = `
    <div class="multi-sample-command">
      <div>
        <p class="report-kicker">Batch Overview</p>
        <h3>${escapeHtml(analysisLabel)}多样本概览</h3>
      </div>
      <div class="multi-sample-command-stats">
        <span><strong>${escapeHtml(String(readyCount))}</strong> / ${escapeHtml(String(sampleCount))} 已生成</span>
        <span>主导物种：<strong>${escapeHtml(topSpecies)}</strong></span>
      </div>
      ${renderTableExportToolbar()}
    </div>
    <div class="multi-sample-table-frame">
      <table class="multi-sample-table">
        <thead>
          <tr>${columns.map((column) => `<th>${escapeHtml(column.label || column.key || "")}</th>`).join("")}</tr>
        </thead>
        <tbody>
          ${rows.map((row) => `
            <tr class="${rowClass(row)}" data-sample-url="${escapeHtml(sampleUrl(row))}" title="双击进入 ${escapeHtml(row.sample || "样本")} 的完整报告">
              ${columns.map((column) => `<td data-column="${escapeHtml(column.key || "")}">${renderCell(row, column)}</td>`).join("")}
            </tr>
          `).join("")}
        </tbody>
      </table>
    </div>
    <div class="multi-sample-footnote">
      <span>默认先显示批次判读表；点击样本名或双击整行进入完整单样本报告。</span>
      <span>关注行会按低完整性、高污染率或病毒分型 QC 异常自动标记。</span>
    </div>
  `;
  bindTableExportButtons(container, exportTitle, exportColumns, exportRows);
  container.querySelectorAll(".multi-sample-table tbody tr[data-sample-url]").forEach((row) => {
    row.addEventListener("dblclick", (event) => {
      const target = event.target;
      if (target instanceof HTMLElement && target.closest("button, a, input, select, textarea")) return;
      const url = row.getAttribute("data-sample-url");
      if (url) window.location.href = url;
    });
  });
  return true;
}

function applyMultiSampleLandingLayout(data) {
  const shell = document.querySelector(".report-shell");
  const landing = isMultiSampleLanding(data);
  if (shell) shell.dataset.reportMode = landing ? "multi-overview" : String(data?.task?.report_mode || "single");
  document.querySelector(".report-summary-strip")?.classList.toggle("hidden", landing);
  document.querySelectorAll(".report-section").forEach((section) => {
    if (!(section instanceof HTMLElement)) return;
    section.classList.toggle("hidden", landing && section.id !== "section-overview");
  });
  document.querySelectorAll(".report-nav-group").forEach((group) => {
    const overviewLink = group.querySelector?.('a[href="#section-overview"]');
    group.classList.toggle("hidden", landing && !overviewLink);
  });
  if (landing) {
    const heading = document.querySelector("#section-overview .section-heading");
    const title = heading?.querySelector("h2");
    const copy = heading?.querySelector("p:last-child");
    if (title) title.textContent = "多样本结果概览";
    if (copy) copy.textContent = "一行一个样本，先比较 QC、物种鉴定、组装和分型/血清型等关键判读信息，再下钻单样本详情。";
  }
  return landing;
}

function renderChartLegend(seriesList) {
  if (!Array.isArray(seriesList) || seriesList.length === 0) return "";
  return `
    <div class="chart-legend" aria-label="图例开关">
      ${seriesList.map((series, index) => `
        <button class="chart-legend-item" type="button" data-series-toggle="${index}" aria-pressed="true">
          <i class="chart-legend-swatch" style="--legend-color:${escapeHtml(series.color)}"></i>
          <span>${escapeHtml(series.label)}</span>
        </button>
      `).join("")}
    </div>
  `;
}

function createSeriesSvg(seriesList, options = {}) {
  const width = options.width || 1100;
  const height = options.height || 320;
  const padX = options.padX ?? 72;
  const padTop = options.padTop ?? 20;
  const padBottom = options.padBottom ?? 58;
  const innerWidth = width - padX * 2;
  const innerHeight = height - padTop - padBottom;
  const maxLength = Math.max(...seriesList.map((series) => (series.values || []).length), 0);
  const allValues = seriesList.flatMap((series) => (series.values || []).map((value) => Number(value) || 0));
  const computedMax = Math.max(...allValues, 1);
  const minValue = options.minValue ?? 0;
  const roundedMax = options.maxValue ?? (computedMax <= 1 ? 1 : Math.ceil(computedMax * 1.1));
  const yTicks = 4;
  const grid = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = minValue + (roundedMax - minValue) * ratio;
    return `
      <line class="chart-grid-line" x1="${padX}" y1="${y}" x2="${width - padX}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padX - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value, options.yFormatter))}</text>
    `;
  }).join("");
  const xValues = Array.isArray(options.xValues) ? options.xValues : [];
  const xTickValues = options.xTicks || [1, Math.max(1, Math.round((maxLength + 1) / 2)), Math.max(1, maxLength)];
  const xTicks = [0, 0.5, 1].map((ratio, index) => {
    const x = padX + innerWidth * ratio;
    return `
      <line class="chart-axis-tick" x1="${x}" y1="${height - padBottom}" x2="${x}" y2="${height - padBottom + 6}"></line>
      <text class="chart-axis-label x-axis" x="${x}" y="${height - 26}">${escapeHtml(String(xTickValues[index] ?? ""))}</text>
    `;
  }).join("");
  const endLabels = Array.isArray(options.endLabels) ? options.endLabels : [];
  const pathChunks = [];
  const endLabelChunks = [];
  seriesList.forEach((series, seriesIndex) => {
    const values = series.values || [];
    const coords = values.map((value, pointIndex) => {
      const x = padX + (maxLength <= 1 ? 0 : (pointIndex / (maxLength - 1)) * innerWidth);
      const rawNumber = Number(value) || 0;
      const scaledValue = options.yFormatter === "percent" ? rawNumber * 100 : rawNumber;
      const normalized = (scaledValue - minValue) / Math.max(roundedMax - minValue, 1);
      const y = padTop + innerHeight - normalized * innerHeight;
      return { x, y };
    });
    const points = coords.map(({ x, y }) => `${x},${y}`).join(" ");
    const areaPoints = [
      `${padX},${height - padBottom}`,
      ...coords.map(({ x, y }) => `${x},${y}`),
      `${padX + innerWidth},${height - padBottom}`,
    ].join(" ");
    pathChunks.push(`
      <polygon class="chart-area-fill" data-series-index="${seriesIndex}" fill="${series.color}" fill-opacity="${seriesList.length > 2 ? "0.04" : "0.08"}" points="${areaPoints}"></polygon>
      <polyline class="chart-series-line" data-series-index="${seriesIndex}" fill="none" stroke="${series.color}" stroke-width="${seriesList.length > 2 ? "2.1" : "2.4"}" points="${points}"></polyline>
    `);
    const endLabel = endLabels[seriesIndex];
    const last = coords[coords.length - 1];
    if (endLabel && last) {
      endLabelChunks.push(`
        <g class="chart-end-label" data-series-index="${seriesIndex}">
          <circle cx="${last.x}" cy="${last.y}" r="4.2" fill="${series.color}"></circle>
          <rect x="${Math.max(last.x - 18, padX + 8)}" y="${Math.max(last.y - 34, padTop + 4)}" width="54" height="22" rx="11" fill="rgba(255,252,247,0.95)" stroke="${series.color}" stroke-opacity="0.28"></rect>
          <text x="${Math.max(last.x + 9, padX + 35)}" y="${Math.max(last.y - 19, padTop + 18)}" class="chart-end-label-text">${escapeHtml(String(endLabel))}</text>
        </g>
      `);
    }
  });
  const paths = pathChunks.join("");
  const endLabelsSvg = endLabelChunks.join("");
  const focusDots = seriesList.map((series, seriesIndex) => (
    `<circle class="chart-focus-dot" data-series-index="${seriesIndex}" fill="${series.color}" r="4" cx="${padX}" cy="${height - padBottom}"></circle>`
  )).join("");
  return `
    <div class="chart-shell">
      ${renderChartLegend(seriesList)}
      <div
        class="interactive-chart"
        data-chart-padx="${padX}"
        data-chart-padtop="${padTop}"
        data-chart-padbottom="${padBottom}"
        data-chart-innerwidth="${innerWidth}"
        data-chart-innerheight="${innerHeight}"
        data-chart-maxlength="${maxLength}"
        data-chart-minvalue="${minValue}"
        data-chart-maxvalue="${roundedMax}"
        data-chart-xlabel="${escapeHtml(options.xLabel || "位置") }"
        data-chart-ylabel="${escapeHtml(options.yLabel || "数值") }"
        data-chart-yformatter="${escapeHtml(options.yFormatter || "number") }"
        data-chart-xvalues="${escapeHtml(JSON.stringify(xValues))}"
        data-chart-series="${escapeHtml(JSON.stringify(seriesList.map((series) => ({ label: series.label, color: series.color, values: series.values || [] }))))}"
      >
        <div class="chart-canvas" style="--chart-height:${height}px">
          <svg class="sparkline-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(options.label || "结果曲线")}">
            ${grid}
            <line class="chart-axis-line" x1="${padX}" y1="${padTop}" x2="${padX}" y2="${height - padBottom}"></line>
            <line class="chart-axis-line" x1="${padX}" y1="${height - padBottom}" x2="${width - padX}" y2="${height - padBottom}"></line>
            ${xTicks}
            ${paths}
            ${endLabelsSvg}
            <line class="chart-focus-line" x1="${padX}" y1="${padTop}" x2="${padX}" y2="${height - padBottom}"></line>
            ${focusDots}
            <text class="chart-axis-title" x="${width / 2}" y="${height - 6}">${escapeHtml(options.xLabel || "位置")}</text>
            <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(options.yLabel || "数值")}</text>
          </svg>
        </div>
        <div class="chart-tooltip" hidden></div>
      </div>
    </div>
  `;
}

function renderStatsWithCharts(containerId, section, readLabel) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.classList.remove("empty-box");
  const qualitySvg = createSeriesSvg([
    { label: `${readLabel} mean`, color: "#4e6177", values: section.quality_curves?.mean || [] },
  ], { label: `${readLabel} 数据质量`, width: 820, height: 280, xLabel: "碱基位置", yLabel: "质量分数", maxValue: 45 });
  const cards = [
    { label: "总 reads", value: String(section.before_summary?.total_reads ?? "-") },
    { label: "总数据量", value: formatBases(section.before_summary?.total_bases) },
    { label: "平均长度", value: String(section.before_summary?.mean_length ?? "-") },
    { label: "Q20", value: formatRate(section.before_summary?.q20_rate) },
    { label: "Q30", value: formatRate(section.before_summary?.q30_rate) },
    { label: "GC", value: formatRate(section.before_summary?.gc_content) },
  ];
  container.innerHTML = `
    <div class="mini-stat-grid">
      ${cards.map((item) => `
        <div class="mini-stat-card">
          <span>${escapeHtml(item.label)}</span>
          <strong>${escapeHtml(item.value)}</strong>
        </div>
      `).join("")}
    </div>
    <div class="mini-chart-stack">
      <div class="mini-chart-card">
        <span class="mini-chart-title">数据质量</span>
        ${buildChartInsight(summarizeRawQc(section, readLabel))}
        ${qualitySvg}
      </div>
    </div>
  `;
}

function renderAdapterTab(panel, adapter, isPaired) {
  const hasAdapterData = Number(adapter?.adapter_trimmed_reads || 0) > 0
    || Number(adapter?.adapter_trimmed_bases || 0) > 0
    || String(adapter?.read1_adapter_sequence || "").trim()
    || String(adapter?.read2_adapter_sequence || "").trim();
  panel.innerHTML = hasAdapterData ? `
    <div class="adapter-grid">
      <div class="adapter-card"><span>接头去除 reads</span><strong>${escapeHtml(String(adapter.adapter_trimmed_reads ?? "-"))}</strong></div>
      <div class="adapter-card"><span>接头去除碱基</span><strong>${escapeHtml(formatBases(adapter.adapter_trimmed_bases))}</strong></div>
      <div class="adapter-card"><span>Read1 接头序列</span><code>${escapeHtml(adapter.read1_adapter_sequence || "-")}</code></div>
      ${isPaired ? `<div class="adapter-card"><span>Read2 接头序列</span><code>${escapeHtml(adapter.read2_adapter_sequence || "-")}</code></div>` : ""}
    </div>
  ` : `
    <div class="empty-box">
      <p>当前 fastp 结果未生成可用的接头统计。</p>
    </div>
  `;
}

function renderFastpTabs(fastp, isPaired) {
  const insertPanel = document.getElementById("fastp-tab-insert-size");
  const basePanel = document.getElementById("fastp-tab-base-content");
  const adapterPanel = document.getElementById("fastp-tab-adapter");
  if (insertPanel) {
    const hasInsertData = Array.isArray(fastp.insert_size?.histogram) && fastp.insert_size.histogram.length;
    insertPanel.innerHTML = isPaired && hasInsertData ? `
      <div class="mini-chart-card">
        <span class="mini-chart-title">插入片段长度</span>
        ${buildChartInsight(summarizeInsertSize(fastp))}
        ${createSeriesSvg([{ label: "insert size", color: "#4e6177", values: fastp.insert_size?.histogram || [] }], { label: "插入片段长度", width: 820, height: 240, xLabel: "插入片段长度区间", yLabel: "read 数" })}
        <p class="empty-copy">峰值：${escapeHtml(String(fastp.insert_size?.peak ?? "-"))}；未知配对：${escapeHtml(String(fastp.insert_size?.unknown ?? "-"))}</p>
      </div>
    ` : `
      <div class="empty-box">
        <p>${isPaired ? "当前 fastp 结果未生成插入片段长度分布。" : "单端数据不生成插入片段长度分布。"}</p>
      </div>
    `;
  }
  if (basePanel) {
    basePanel.innerHTML = isPaired ? `
      <div class="mini-chart-card">
        <div class="subreport-tabs" role="tablist" aria-label="碱基分布切换">
          <button class="subreport-tab-button active" type="button" data-base-tab="read1">R1</button>
          <button class="subreport-tab-button" type="button" data-base-tab="read2">R2</button>
        </div>
        <div class="subreport-tab-panel active" data-base-panel="read1">
          <span class="mini-chart-title">R1 碱基分布</span>
          ${buildChartInsight(summarizeBaseDistribution(fastp.base_distribution?.read1 || {}, "R1"))}
          ${createSeriesSvg([
            { label: "A", color: "#8a6654", values: fastp.base_distribution?.read1?.A || [] },
            { label: "T", color: "#7a7158", values: fastp.base_distribution?.read1?.T || [] },
            { label: "C", color: "#5d7c83", values: fastp.base_distribution?.read1?.C || [] },
            { label: "G", color: "#6d6481", values: fastp.base_distribution?.read1?.G || [] },
            { label: "GC", color: "#3e546f", values: fastp.base_distribution?.read1?.GC || [] },
          ], { label: "R1 碱基分布", width: 820, height: 240, xLabel: "碱基位置", yLabel: "比例", yFormatter: "percent", maxValue: 100 })}
        </div>
        <div class="subreport-tab-panel" data-base-panel="read2">
          <span class="mini-chart-title">R2 碱基分布</span>
          ${buildChartInsight(summarizeBaseDistribution(fastp.base_distribution?.read2 || {}, "R2"))}
          ${createSeriesSvg([
            { label: "A", color: "#8a6654", values: fastp.base_distribution?.read2?.A || [] },
            { label: "T", color: "#7a7158", values: fastp.base_distribution?.read2?.T || [] },
            { label: "C", color: "#5d7c83", values: fastp.base_distribution?.read2?.C || [] },
            { label: "G", color: "#6d6481", values: fastp.base_distribution?.read2?.G || [] },
            { label: "GC", color: "#3e546f", values: fastp.base_distribution?.read2?.GC || [] },
          ], { label: "R2 碱基分布", width: 820, height: 240, xLabel: "碱基位置", yLabel: "比例", yFormatter: "percent", maxValue: 100 })}
        </div>
      </div>
    ` : `
      <div class="mini-chart-card">
        <span class="mini-chart-title">单端碱基分布</span>
        ${buildChartInsight(summarizeBaseDistribution(fastp.base_distribution?.read1 || {}, "单端 reads"))}
        ${createSeriesSvg([
          { label: "A", color: "#8a6654", values: fastp.base_distribution?.read1?.A || [] },
          { label: "T", color: "#7a7158", values: fastp.base_distribution?.read1?.T || [] },
          { label: "C", color: "#5d7c83", values: fastp.base_distribution?.read1?.C || [] },
          { label: "G", color: "#6d6481", values: fastp.base_distribution?.read1?.G || [] },
          { label: "GC", color: "#3e546f", values: fastp.base_distribution?.read1?.GC || [] },
        ], { label: "单端碱基分布", width: 820, height: 240, xLabel: "碱基位置", yLabel: "比例", yFormatter: "percent", maxValue: 100 })}
      </div>
    `;
  }
  if (adapterPanel) {
    renderAdapterTab(adapterPanel, fastp.adapter_cutting || {}, isPaired);
  }
}

function initializeBaseTabs(root = document) {
  root.querySelectorAll(".subreport-tabs").forEach((tabGroup) => {
    if (tabGroup.dataset.initialized === "true") return;
    tabGroup.dataset.initialized = "true";
    const container = tabGroup.closest(".mini-chart-card");
    const buttons = Array.from(tabGroup.querySelectorAll(".subreport-tab-button"));
    const panels = Array.from(container?.querySelectorAll(".subreport-tab-panel") || []);
    buttons.forEach((button) => {
      button.addEventListener("click", () => {
        const target = button.dataset.baseTab;
        buttons.forEach((item) => item.classList.toggle("active", item === button));
        panels.forEach((panel) => {
          panel.classList.toggle("active", panel.dataset.basePanel === target);
        });
        initializeInteractiveCharts(container || document);
      });
    });
  });
}

function initializeInteractiveCharts(root = document) {
  root.querySelectorAll(".interactive-chart").forEach((chart) => {
    if (chart.dataset.initialized === "true") return;
    chart.dataset.initialized = "true";
    const svg = chart.querySelector(".sparkline-svg");
    const tooltip = chart.querySelector(".chart-tooltip");
    const focusLine = chart.querySelector(".chart-focus-line");
    const dots = Array.from(chart.querySelectorAll(".chart-focus-dot"));
    const seriesList = JSON.parse(chart.dataset.chartSeries || "[]");
    const maxLength = Number(chart.dataset.chartMaxlength || 0);
    const padX = Number(chart.dataset.chartPadx || 0);
    const padTop = Number(chart.dataset.chartPadtop || 0);
    const padBottom = Number(chart.dataset.chartPadbottom || 0);
    const innerWidth = Number(chart.dataset.chartInnerwidth || 0);
    const innerHeight = Number(chart.dataset.chartInnerheight || 0);
    const minValue = Number(chart.dataset.chartMinvalue || 0);
    const maxValue = Number(chart.dataset.chartMaxvalue || 1);
    const yFormatter = chart.dataset.chartYformatter || "number";
    const xLabel = chart.dataset.chartXlabel || "位置";
    const xValues = JSON.parse(chart.dataset.chartXvalues || "[]");
    const visibleSeries = new Set(seriesList.map((_, index) => index));
    if (!svg || !tooltip || maxLength < 1) return;

    const applySeriesVisibility = () => {
      chart.querySelectorAll("[data-series-index]").forEach((node) => {
        const index = Number(node.dataset.seriesIndex);
        const visible = visibleSeries.has(index);
        node.classList.toggle("is-series-hidden", !visible);
      });
      chart.querySelectorAll("[data-series-toggle]").forEach((button) => {
        const index = Number(button.dataset.seriesToggle);
        const active = visibleSeries.has(index);
        button.classList.toggle("is-inactive", !active);
        button.setAttribute("aria-pressed", active ? "true" : "false");
      });
    };

    const update = (clientX) => {
      const rect = svg.getBoundingClientRect();
      const relativeX = Math.min(Math.max(clientX - rect.left, 0), rect.width);
      const ratio = rect.width <= 0 ? 0 : relativeX / rect.width;
      const index = Math.min(Math.max(Math.round(ratio * Math.max(maxLength - 1, 0)), 0), Math.max(maxLength - 1, 0));
      const chartX = padX + (maxLength <= 1 ? 0 : (index / (maxLength - 1)) * innerWidth);
      focusLine.setAttribute("x1", chartX);
      focusLine.setAttribute("x2", chartX);
      const rows = [];
      dots.forEach((dot, seriesIndex) => {
        if (!visibleSeries.has(seriesIndex)) {
          dot.classList.add("is-series-hidden");
          return;
        }
        const series = seriesList[seriesIndex] || {};
        const rawValue = Number(series.values?.[index] ?? 0);
        const scaledValue = yFormatter === "percent" ? rawValue * 100 : rawValue;
        const normalized = (scaledValue - minValue) / Math.max(maxValue - minValue, 1);
        const y = padTop + innerHeight - normalized * innerHeight;
        dot.setAttribute("cx", chartX);
        dot.setAttribute("cy", y);
        dot.classList.remove("is-series-hidden");
        rows.push(`<span><i style="background:${series.color}"></i>${escapeHtml(series.label)}: ${escapeHtml(formatChartValue(rawValue, yFormatter))}</span>`);
      });
      const xValue = xValues[index] ?? (index + 1);
      if (!rows.length) {
        tooltip.hidden = true;
        return;
      }
      tooltip.innerHTML = `<strong>${escapeHtml(`${xLabel}: ${xValue}`)}</strong>${rows.join("")}`;
      tooltip.hidden = false;
    };

    chart.querySelectorAll("[data-series-toggle]").forEach((button) => {
      button.addEventListener("click", () => {
        const index = Number(button.dataset.seriesToggle);
        if (visibleSeries.has(index) && visibleSeries.size > 1) visibleSeries.delete(index);
        else visibleSeries.add(index);
        applySeriesVisibility();
        tooltip.hidden = true;
      });
    });

    svg.addEventListener("mousemove", (event) => update(event.clientX));
    svg.addEventListener("mouseenter", (event) => update(event.clientX));
    svg.addEventListener("mouseleave", () => {
      tooltip.hidden = true;
    });
    applySeriesVisibility();
  });
}

function renderRiskBars(items, type) {
  return `
    <div class="risk-bar-list risk-bar-list-${type}">
      ${items.map((item) => {
        const width = Math.max(Math.min(Number(item.ratio) || 0, 100), 2);
        return `
          <div class="risk-bar-item">
            <div class="risk-bar-meta">
              <span>${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(String(item.reads))} reads</strong>
            </div>
            <div class="risk-bar-track" aria-hidden="true">
              <div class="risk-bar-fill" style="width:${width}%"></div>
            </div>
            <div class="risk-bar-foot">
              <span>占比 ${escapeHtml(formatRate((Number(item.ratio) || 0) / 100))}</span>
              <span>分类项 ${escapeHtml(String(item.records || 0))}</span>
            </div>
          </div>
        `;
      }).join("")}
    </div>
  `;
}

function renderKingdomCompareBars(items) {
  const totalReads = items.reduce((sum, item) => sum + (Number(item.reads) || 0), 0);
  return `
    <div class="kingdom-compare-list">
      ${items.map((item) => {
        const reads = Number(item.reads) || 0;
        const ratio = totalReads ? (reads / totalReads) * 100 : Number(item.ratio) || 0;
        const width = Math.max(Math.min(ratio, 100), reads > 0 ? 8 : 2);
        return `
          <div class="kingdom-compare-item">
            <div class="kingdom-compare-meta">
              <span class="kingdom-compare-label">${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(String(reads))} reads</strong>
            </div>
            <div class="kingdom-compare-track" aria-hidden="true">
              <div class="kingdom-compare-fill" style="width:${width}%"></div>
            </div>
            <div class="kingdom-compare-foot">
              <span>占比 ${escapeHtml(formatRate(ratio / 100))}</span>
              <span>${escapeHtml(String(item.records || 0))} 项分类</span>
            </div>
          </div>
        `;
      }).join("")}
    </div>
  `;
}

function renderRvCategoryBars(items, type) {
  return `
    <div class="risk-bar-list risk-bar-list-${type}">
      ${items.map((item) => {
        const width = Math.max(Math.min((Number(item.count) || 0) * 10, 100), 4);
        return `
          <div class="risk-bar-item">
            <div class="risk-bar-meta">
              <span>${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(String(item.count))} 项</strong>
            </div>
            <div class="risk-bar-track" aria-hidden="true">
              <div class="risk-bar-fill" style="width:${width}%"></div>
            </div>
          </div>
        `;
      }).join("")}
    </div>
  `;
}

function renderResistanceVirulenceOverview(section) {
  const resistanceContainer = document.getElementById("resistance-summary-panel");
  const virulenceContainer = document.getElementById("virulence-summary-panel");
  if (!resistanceContainer || !virulenceContainer) return;
  const resistance = section?.resistance || {};
  const virulence = section?.virulence || {};

  resistanceContainer.classList.remove("empty-box");
  virulenceContainer.classList.remove("empty-box");

  resistanceContainer.innerHTML = `
    ${buildChartInsight(resistance.note || "当前未检出可汇总的耐药元件结果。")}
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">耐药药物类别分布</span>
        <span class="risk-summary-count">${escapeHtml(String(resistance.gene_count || 0))} 个基因</span>
      </div>
      <div class="mini-stat-grid">
        <div class="mini-stat-card"><span>元件记录</span><strong>${escapeHtml(String(resistance.hit_count || 0))}</strong></div>
        <div class="mini-stat-card"><span>汇总基因数</span><strong>${escapeHtml(String(resistance.summary_count || 0))}</strong></div>
      </div>
      ${(resistance.top_categories || []).length ? renderRvCategoryBars(resistance.top_categories || [], "hazard") : `<p class="empty-copy">未检出可展示的耐药药物类别分布。</p>`}
    </div>
  `;

  virulenceContainer.innerHTML = `
    ${buildChartInsight(virulence.note || "当前未检出可汇总的毒力元件结果。")}
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">VF 分类分布</span>
        <span class="risk-summary-count">${escapeHtml(String(virulence.gene_count || 0))} 个基因</span>
      </div>
      <div class="mini-stat-grid">
        <div class="mini-stat-card"><span>元件记录</span><strong>${escapeHtml(String(virulence.hit_count || 0))}</strong></div>
        <div class="mini-stat-card"><span>汇总基因数</span><strong>${escapeHtml(String(virulence.summary_count || 0))}</strong></div>
      </div>
      ${(virulence.top_categories || []).length ? renderRvCategoryBars(virulence.top_categories || [], "pathogenicity") : `<p class="empty-copy">未检出可展示的 VF 分类分布。</p>`}
    </div>
  `;
}

function renderMgeSummaryBlock(containerId, block, typeLabel) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const defaultRiskLevels = [
    { label: "Level A 高迁移风险", count: 0 },
    { label: "Level B 中迁移风险", count: 0 },
    { label: "Level C 低到中等风险", count: 0 },
    { label: "Level D 弱证据", count: 0 },
  ];
  const normalizedRiskLevels = defaultRiskLevels.map((item) => {
    const matched = (block?.risk_levels || []).find((risk) => risk?.label === item.label);
    return { label: item.label, count: Number(matched?.count) || 0 };
  });
  if (!block || !Array.isArray(block.rows) || !block.rows.length) {
    container.innerHTML = `
      ${buildChartInsight(`${typeLabel}未检出可展示的移动元件监测结果。`)}
      <div class="risk-summary-panel">
        <div class="risk-summary-head">
          <span class="mini-chart-title">${escapeHtml(typeLabel)}监测摘要</span>
          <span class="risk-summary-count">0 条记录</span>
        </div>
        <div class="mini-stat-grid">
          <div class="mini-stat-card"><span>监测记录</span><strong>0</strong></div>
          <div class="mini-stat-card"><span>Level A/B</span><strong>0</strong></div>
        </div>
        ${renderRvCategoryBars(normalizedRiskLevels, "hazard")}
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  container.innerHTML = `
    ${buildChartInsight(block.note || "")}
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">${escapeHtml(typeLabel)}监测摘要</span>
        <span class="risk-summary-count">${escapeHtml(String(block.hit_count || 0))} 条记录</span>
      </div>
      <div class="mini-stat-grid">
        <div class="mini-stat-card"><span>监测记录</span><strong>${escapeHtml(String(block.hit_count || 0))}</strong></div>
        <div class="mini-stat-card"><span>Level A/B</span><strong>${escapeHtml(String(block.high_risk_count || 0))}</strong></div>
      </div>
      ${renderRvCategoryBars(normalizedRiskLevels, "hazard")}
    </div>
  `;
}

function renderMgeBarChart(containerId, title, items, xLabel, yLabel) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (!Array.isArray(items) || !items.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>${escapeHtml(title)}暂无可展示数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  container.innerHTML = renderBarSvg(
    items.map((item) => Number(item.count) || 0),
    {
      label: title,
      xLabel,
      yLabel,
      xValues: items.map((item) => item.label),
      width: 1120,
      height: 360,
      padBottom: 118,
    },
  );
}

function renderMgePieChart(containerId, title, items) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (!Array.isArray(items) || !items.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>${escapeHtml(title)}暂无可展示数据。</p>
      </div>
    `;
    return;
  }
  const total = items.reduce((sum, item) => sum + (Number(item.count) || 0), 0);
  if (!total) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>${escapeHtml(title)}暂无可展示数据。</p>
      </div>
    `;
    return;
  }
  const colors = ["#355c94", "#7b8f48", "#a56a43", "#6d6481", "#4d7d76", "#9a4f5b", "#c68c3a", "#58708f"];
  let current = 0;
  const segments = items.map((item, index) => {
    const value = Number(item.count) || 0;
    const ratio = value / total;
    const dash = (ratio * 251.2).toFixed(2);
    const offset = (-current * 251.2).toFixed(2);
    current += ratio;
    return {
      label: item.label,
      count: value,
      ratio,
      color: colors[index % colors.length],
      dash,
      offset,
    };
  });
  container.classList.remove("empty-box");
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card">
      <span class="mini-chart-title">${escapeHtml(title)}</span>
      <div class="mge-pie-layout">
        <div class="donut-chart mge-donut-chart" aria-label="${escapeHtml(title)}">
          <svg viewBox="0 0 120 120" role="img">
            <circle class="donut-track" cx="60" cy="60" r="40"></circle>
            ${segments.map((segment) => `
              <circle
                class="donut-segment"
                cx="60"
                cy="60"
                r="40"
                stroke="${segment.color}"
                stroke-dasharray="${segment.dash} 251.2"
                stroke-dashoffset="${segment.offset}">
                <title>${escapeHtml(`${segment.label}: ${segment.count} (${(segment.ratio * 100).toFixed(1)}%)`)}</title>
              </circle>
            `).join("")}
          </svg>
          <div class="donut-total">${escapeHtml(String(total))}</div>
        </div>
        <div class="mge-pie-legend">
          ${segments.map((segment) => `
            <div class="mge-pie-legend-item">
              <span class="mge-pie-swatch" style="--mge-color:${escapeHtml(segment.color)}"></span>
              <span>${escapeHtml(segment.label)}</span>
              <strong>${escapeHtml(String(segment.count))}</strong>
            </div>
          `).join("")}
        </div>
      </div>
    </div>
  `;
}

function renderMgeOverview(section) {
  const summaryContainer = document.getElementById("mge-overview-summary");
  if (!summaryContainer) return;
  const overview = section?.overview || {};
  if (!section || section.status === "empty") {
    summaryContainer.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>当前任务未检出可展示的移动元件监测结果。</p>
      </div>
    `;
    renderMgePieChart("mge-type-chart", "移动元件类型统计", []);
    renderMgeBarChart("mge-risk-chart", "基因转移风险等级分布", [], "风险等级", "记录数");
    renderMgeSummaryBlock("mge-resistance-summary", null, "耐药相关位点");
    renderMgeSummaryBlock("mge-virulence-summary", null, "毒力相关位点");
    renderCategoryGeneRelationship("mge-resistance-relationship-chart", {});
    renderCategoryGeneRelationship("mge-virulence-relationship-chart", {});
    return;
  }

  summaryContainer.classList.remove("empty-box");
  summaryContainer.innerHTML = `
    ${buildChartInsight(overview.note || "")}
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">移动元件监测概览</span>
        <span class="risk-summary-count">${escapeHtml(String(overview.total_hits || 0))} 条记录</span>
      </div>
      <div class="mini-stat-grid">
        <div class="mini-stat-card"><span>总记录数</span><strong>${escapeHtml(String(overview.total_hits || 0))}</strong></div>
        <div class="mini-stat-card"><span>耐药相关</span><strong>${escapeHtml(String(overview.resistance_hits || 0))}</strong></div>
        <div class="mini-stat-card"><span>毒力相关</span><strong>${escapeHtml(String(overview.virulence_hits || 0))}</strong></div>
      </div>
    </div>
  `;

  renderMgePieChart("mge-type-chart", "移动元件类型统计", overview.mge_types || []);
  renderMgeBarChart("mge-risk-chart", "基因转移风险等级分布", overview.risk_levels || [], "风险等级", "记录数");

  renderMgeSummaryBlock("mge-resistance-summary", section.resistance || {}, "耐药相关位点");
  renderMgeSummaryBlock("mge-virulence-summary", section.virulence || {}, "毒力相关位点");
  renderCategoryGeneRelationship("mge-resistance-relationship-chart", section?.resistance?.gene_mge_relationship || {});
  renderCategoryGeneRelationship("mge-virulence-relationship-chart", section?.virulence?.gene_mge_relationship || {});
}

function renderTaxonomyRiskSummary(section) {
  const summaryContainer = document.getElementById("taxonomy-risk-summary");
  const hazardContainer = document.getElementById("taxonomy-hazard-chart");
  if (!summaryContainer || !hazardContainer) return;
  const kingdomSummary = Array.isArray(section?.kingdom_summary) ? section.kingdom_summary : [];
  const pathogenicity = Array.isArray(section?.pathogenicity) ? section.pathogenicity : [];
  const hazard = Array.isArray(section?.hazard) ? section.hazard : [];
  const narrative = section?.narrative || "";
  const dominantVirus = section?.dominant_virus || null;
  const dominantFungus = section?.dominant_fungus || null;

  if (!pathogenicity.length && !hazard.length) {
    summaryContainer.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>当前序列物种鉴定结果未提供可汇总的致病性与危害等级信息。</p>
      </div>
    `;
    hazardContainer.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的危害程度等级分布。</p>
      </div>
    `;
    return;
  }

  summaryContainer.classList.remove("empty-box");
  hazardContainer.classList.remove("empty-box");

  summaryContainer.innerHTML = `
    <div class="taxonomy-summary-stack">
      ${buildChartInsight(narrative)}
      <div class="taxonomy-summary-feature">
        <div class="risk-summary-head taxonomy-summary-feature-head">
          <div>
            <span class="mini-chart-title">细菌 / 病毒 / 真菌对比</span>
            <p class="taxonomy-summary-caption">基于 3.1 序列物种鉴定结果，对不同大类微生物的序列数量与分类记录进行汇总。</p>
          </div>
          <span class="risk-summary-count">${escapeHtml(String(section?.total_records || 0))} 条分类记录</span>
        </div>
        ${renderKingdomCompareBars(kingdomSummary)}
      </div>
      ${dominantVirus ? `
        <div class="taxonomy-summary-subpanel">
          <div class="risk-summary-head">
            <span class="mini-chart-title">主导病毒标准身份</span>
            <span class="risk-summary-count">${escapeHtml(dominantVirus.taxid || "--")}</span>
          </div>
          <div class="metric-card paired-metric-card taxonomy-virus-identity-card">
            <span class="metric-label">${escapeHtml(dominantVirus.species || dominantVirus.scientific_name || "病毒")}</span>
            <div class="paired-metric-grid">
              <div class="paired-metric-item metric-state-neutral">
                <span>NCBI学名</span>
                <strong>${escapeHtml(dominantVirus.scientific_name || "--")}</strong>
              </div>
              <div class="paired-metric-item metric-state-neutral">
                <span>科 / 属</span>
                <strong>${escapeHtml(`${dominantVirus.family || "-"} / ${dominantVirus.genus || "-"}`)}</strong>
              </div>
              <div class="paired-metric-item metric-state-neutral">
                <span>种</span>
                <strong>${escapeHtml(dominantVirus.species_rank || dominantVirus.species || "--")}</strong>
              </div>
              <div class="paired-metric-item metric-state-neutral">
                <span>序列占比 / 读取</span>
                <strong>${escapeHtml(`${Number(dominantVirus.ratio || 0).toFixed(2)}% / ${dominantVirus.reads || 0}`)}</strong>
              </div>
            </div>
          </div>
        </div>
      ` : ""}
      ${dominantFungus ? `
        <div class="taxonomy-summary-subpanel">
          <div class="risk-summary-head">
            <span class="mini-chart-title">主导真菌标准身份</span>
            <span class="risk-summary-count">${escapeHtml(dominantFungus.taxid || "--")}</span>
          </div>
          <div class="metric-card paired-metric-card taxonomy-virus-identity-card">
            <span class="metric-label">${escapeHtml(dominantFungus.species || dominantFungus.scientific_name || "真菌")}</span>
            <div class="paired-metric-grid">
              <div class="paired-metric-item metric-state-neutral">
                <span>NCBI学名</span>
                <strong>${escapeHtml(dominantFungus.scientific_name || "--")}</strong>
              </div>
              <div class="paired-metric-item metric-state-neutral">
                <span>科 / 属</span>
                <strong>${escapeHtml(`${dominantFungus.family || "-"} / ${dominantFungus.genus || "-"}`)}</strong>
              </div>
              <div class="paired-metric-item metric-state-neutral">
                <span>种</span>
                <strong>${escapeHtml(dominantFungus.species_rank || dominantFungus.species || "--")}</strong>
              </div>
              <div class="paired-metric-item metric-state-neutral">
                <span>序列占比 / 读取</span>
                <strong>${escapeHtml(`${Number(dominantFungus.ratio || 0).toFixed(2)}% / ${dominantFungus.reads || 0}`)}</strong>
              </div>
            </div>
          </div>
        </div>
      ` : ""}
      <div class="taxonomy-summary-subpanel">
        <div class="risk-summary-head">
          <span class="mini-chart-title">致病性分布</span>
          <span class="risk-summary-count">${escapeHtml(String(pathogenicity.length))} 类</span>
        </div>
        ${pathogenicity.length ? renderRiskBars(pathogenicity, "pathogenicity") : `<p class="empty-copy">未检出可展示的致病性结果。</p>`}
      </div>
    </div>
  `;

  hazardContainer.innerHTML = `
    <div class="risk-summary-panel">
      <div class="risk-summary-head">
        <span class="mini-chart-title">危害程度等级</span>
        <span class="risk-summary-count">${escapeHtml(String(hazard.length))} 类</span>
      </div>
      ${hazard.length ? renderRiskBars(hazard, "hazard") : `<p class="empty-copy">未检出可展示的危害程度等级结果。</p>`}
    </div>
  `;
}

function renderTaxonomyInterpretation(section) {
  const confusionContainer = document.getElementById("taxonomy-confusion-hint");
  const mixtureContainer = document.getElementById("taxonomy-mixture-hint");
  if (!confusionContainer || !mixtureContainer) return;

  const renderHintCard = (container, block, emptyCopy) => {
    const metrics = Array.isArray(block?.metrics) ? block.metrics : [];
    const evidence = Array.isArray(block?.evidence) ? block.evidence : [];
    const tone = block?.tone || "neutral";
    if (!block || block.status === "empty") {
      container.innerHTML = `
        <div class="empty-box coverage-empty">
          <p>${escapeHtml(emptyCopy)}</p>
        </div>
      `;
      return;
    }
    container.classList.remove("empty-box");
    container.innerHTML = `
      <div class="taxonomy-hint-card taxonomy-hint-${escapeHtml(tone)}">
        <div class="taxonomy-hint-head">
          <span class="taxonomy-hint-badge taxonomy-hint-badge-${escapeHtml(tone)}">${escapeHtml(block?.badge || "提示")}</span>
          <strong>${escapeHtml(block?.headline || "--")}</strong>
        </div>
        ${buildChartInsight(block?.summary || "")}
        ${metrics.length ? `
          <div class="mini-stat-grid">
            ${metrics.map((item) => `
              <div class="mini-stat-card">
                <span>${escapeHtml(item.label || "--")}</span>
                <strong>${escapeHtml(item.value || "--")}</strong>
              </div>
            `).join("")}
          </div>
        ` : ""}
        ${evidence.length ? `
          <ul class="taxonomy-hint-list">
            ${evidence.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
          </ul>
        ` : ""}
      </div>
    `;
  };

  renderHintCard(
    confusionContainer,
    section?.confusion_hint || {},
    "当前物种分类结果不足，暂无法判断是否存在相近物种混淆。",
  );
  renderHintCard(
    mixtureContainer,
    section?.mixture_hint || {},
    "当前物种分类结果不足，暂无法判断样本更偏向单菌还是混菌。",
  );
}

function renderTaxonomyRarefaction(section) {
  const container = document.getElementById("taxonomy-rarefaction-chart");
  if (!container) return;
  const speciesPoints = Array.isArray(section?.species_points) ? section.species_points : [];
  const subspeciesPoints = Array.isArray(section?.subspecies_points) ? section.subspecies_points : [];
  if (!speciesPoints.length && !subspeciesPoints.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的分类稀释曲线数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const maxLength = Math.max(speciesPoints.length, subspeciesPoints.length, 1);
  const xValues = Array.from({ length: maxLength }, (_, index) => {
    const speciesX = Number(speciesPoints[index]?.x || 0);
    const subspeciesX = Number(subspeciesPoints[index]?.x || 0);
    return Math.max(speciesX, subspeciesX, 0);
  });
  const series = [];
  if (speciesPoints.length) {
    series.push({
      label: '种',
      color: '#4e6177',
      values: speciesPoints.map((point) => Number(point.y) || 0),
    });
  }
  if (subspeciesPoints.length) {
    series.push({
      label: '亚种',
      color: '#8a6654',
      values: subspeciesPoints.map((point) => Number(point.y) || 0),
    });
  }
  const note = speciesPoints.length && subspeciesPoints.length
    ? (section?.note || `种与亚种分类均已形成累积曲线，可用于观察随着序列数量增加，分类检出数的变化趋势。`)
    : speciesPoints.length
      ? (section?.note || `当前仅形成种水平的分类稀释曲线。`)
      : (section?.note || `当前仅形成亚种水平的分类稀释曲线。`);
  const endLabels = [];
  if (speciesPoints.length) endLabels.push(section?.species_final_expected ?? speciesPoints[speciesPoints.length - 1]?.y ?? "");
  if (subspeciesPoints.length) endLabels.push(section?.subspecies_final_expected ?? subspeciesPoints[subspeciesPoints.length - 1]?.y ?? "");
  container.innerHTML = `
    <div class="mini-chart-card coverage-chart-card">
      <span class="mini-chart-title">分类稀释曲线</span>
      ${buildChartInsight(note)}
      ${createSeriesSvg(series, {
        label: '分类稀释曲线',
        width: 920,
        height: 280,
        xLabel: section?.x_label || '累计序列数量',
        yLabel: section?.y_label || '累计检出分类数',
        xValues,
        xTicks: [
          formatChartValue(xValues[0] || 0),
          formatChartValue(xValues[Math.max(0, Math.floor((xValues.length - 1) / 2))] || 0),
          formatChartValue(xValues[xValues.length - 1] || 0),
        ],
        maxValue: Math.max(...series.flatMap((item) => item.values), 1),
        endLabels,
      })}
    </div>
  `;
}

function renderSpeciesIdentification(section) {
  const sourceContainer = document.getElementById("taxonomy-source-table");
  const tagsContainer = document.getElementById("taxonomy-rank-tags");
  const tableContainer = document.getElementById("taxonomy-summary-table");
  if (!sourceContainer || !tagsContainer || !tableContainer) return;
  applyTableTone(sourceContainer, "taxonomy-source-table");
  applyTableTone(tableContainer, "taxonomy-summary-table");

  const datasets = {
    species: section?.species || { rows: [], rank_options: [], terminal_column: "种" },
    subspecies: section?.subspecies || { rows: [], rank_options: [], terminal_column: "亚种" },
  };
  const state = {
    tab: datasets.species.rows?.length ? "species" : "subspecies",
    rank: null,
  };

  const aggregateRows = (rows, rank) => {
    const groups = new Map();
    rows.forEach((row) => {
      const key = String(row?.[rank] ?? "").trim() || "未注释";
      const current = groups.get(key) || { ratio: 0, reads: 0 };
      current.ratio += Number(row?.["比例数值"] || 0);
      current.reads += Number(row?.["序列数量数值"] || 0);
      groups.set(key, current);
    });
    return Array.from(groups.entries())
      .map(([name, values]) => [rank, name, `${values.ratio.toFixed(2)}%`, String(values.reads)])
      .sort((left, right) => Number(right[3]) - Number(left[3]));
  };

  const orderColumns = (columns, preferred) => {
    const normalizedColumns = columns.map((column) => String(column || "").toLowerCase());
    const used = new Set();
    const ordered = [];
    preferred.forEach((alias) => {
      const target = String(alias || "").toLowerCase();
      const index = normalizedColumns.findIndex((column, position) => !used.has(position) && column === target);
      if (index >= 0) {
        used.add(index);
        ordered.push(columns[index]);
      }
    });
    columns.forEach((column, index) => {
      if (!used.has(index)) ordered.push(column);
    });
    return ordered;
  };

  const renderSourceTable = (dataset) => {
    const hiddenColumns = new Set(["界", "门", "纲", "目", "科", "比例数值", "序列数量数值"]);
    const rows = dataset.rows || [];
    if (!rows.length) {
      sourceContainer.innerHTML = `
        <div class="empty-table-state">
          <strong>未检出分类明细结果</strong>
          <p class="empty-copy">未检出对应的种或亚种分类结果文件。</p>
        </div>
      `;
      return;
    }
    const sample = rows[0] || {};
    const columns = orderColumns(
      Object.keys(sample).filter((column) => !hiddenColumns.has(column)),
      [
        dataset.terminal_column || "种",
        "属",
        "比例",
        "序列数量",
        "NCBI TaxID",
        "NCBI学名",
        "NCBI分类等级",
        "NCBI目",
        "NCBI科",
        "NCBI属",
        "NCBI种",
        "致病性",
        "危害程度等级",
      ],
    );
    sourceContainer.dataset.exportTitle = state.tab === "species" ? "序列物种鉴定_种" : "序列物种鉴定_亚种";
    renderInteractiveContigTable(
      sourceContainer,
      columns,
      rows.map((row) => columns.map((column) => row?.[column] ?? "")),
    );
  };

  const render = () => {
    const dataset = datasets[state.tab] || { rows: [], rank_options: [], terminal_column: "" };
    const rankOptions = dataset.rank_options || [];
    if (!state.rank || !rankOptions.includes(state.rank)) {
      state.rank = rankOptions[rankOptions.length - 1] || null;
    }

    document.querySelectorAll('[data-taxonomy-tab]').forEach((button) => {
      button.classList.toggle('active', button.dataset.taxonomyTab === state.tab);
    });

    renderSourceTable(dataset);

    tagsContainer.innerHTML = rankOptions.map((rank) => `
      <button class="taxonomy-rank-tag ${state.rank === rank ? 'active' : ''}" type="button" data-taxonomy-rank="${escapeHtml(rank)}">${escapeHtml(rank)}</button>
    `).join('');

    if (!dataset.rows?.length || !state.rank) {
      tableContainer.innerHTML = `
        <div class="empty-table-state">
          <strong>未检出序列物种鉴定结果</strong>
          <p class="empty-copy">未检出对应的种或亚种分类结果文件。</p>
        </div>
      `;
    } else {
      tableContainer.dataset.exportTitle = `${state.tab === "species" ? "种" : "亚种"}_${state.rank}_分类聚合`;
      renderInteractiveContigTable(
        tableContainer,
        ["分类水平", "分类名称", "比例", "序列数量"],
        aggregateRows(dataset.rows, state.rank),
      );
    }

    tagsContainer.querySelectorAll('[data-taxonomy-rank]').forEach((button) => {
      button.addEventListener('click', () => {
        state.rank = button.dataset.taxonomyRank;
        render();
      });
    });
  };

  document.querySelectorAll('[data-taxonomy-tab]').forEach((button) => {
    button.onclick = () => {
      state.tab = button.dataset.taxonomyTab;
      state.rank = null;
      render();
    };
  });

  render();
}

function renderTaxonomyAbundance(section) {
  const container = document.getElementById("taxonomy-abundance-chart");
  const topnSelect = document.getElementById("taxonomy-abundance-topn");
  const tag = document.getElementById("taxonomy-abundance-tag");
  if (!container) return;
  const ranks = Array.isArray(section?.ranks) ? section.ranks : [];
  if (!ranks.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的分类丰度数据。</p>
      </div>
    `;
    return;
  }
  const state = {
    topn: Number(topnSelect?.value || 10),
  };

  const buildSegments = (segments) => {
    const topn = Math.max(1, state.topn || 10);
    const topItems = segments.slice(0, topn);
    const otherItems = segments.slice(topn);
    const merged = topItems.map((segment, index) => ({
      ...segment,
      color: segment.color || ["#526a86", "#76834f", "#8a6654", "#6d6481", "#4e7b75", "#9b7a3f", "#8a4d47", "#5d7c83", "#7b6d5a", "#697789"][index % 10] || '#526a86',
    }));
    if (otherItems.length) {
      merged.push({
        name: 'Other',
        ratio: otherItems.reduce((sum, item) => sum + (Number(item.ratio) || 0), 0),
        reads: otherItems.reduce((sum, item) => sum + (Number(item.reads) || 0), 0),
        color: '#8d8d8d',
      });
    }
    return merged;
  };

  const render = () => {
    container.classList.remove("empty-box");
    if (tag) tag.textContent = `Top ${state.topn} + Other`;
    container.innerHTML = `
      <div class="abundance-stack-list">
        ${ranks.map((rankBlock) => {
          const segments = buildSegments(Array.isArray(rankBlock.segments) ? rankBlock.segments : []);
          return `
            <article class="abundance-rank-card">
              <div class="abundance-rank-head">
                <div>
                  <span class="abundance-rank-label">${escapeHtml(rankBlock.rank)}</span>
                  <strong>分类丰度堆积</strong>
                </div>
                <span class="abundance-rank-total">总占比 ${escapeHtml(formatRate(rankBlock.total_ratio))}</span>
              </div>
              <div class="abundance-bar" role="img" aria-label="${escapeHtml(rankBlock.rank)}分类丰度堆积图">
                ${segments.map((segment) => `
                  <span
                    class="abundance-segment"
                    style="width:${Math.max(Number(segment.ratio) || 0, 0)}%; --segment-color:${escapeHtml(segment.color || '#526a86')}"
                    title="${escapeHtml(`${segment.name} | ${formatRate(segment.ratio)} | ${segment.reads} reads`)}"
                  ></span>
                `).join("")}
              </div>
              <div class="abundance-legend">
                ${segments.map((segment) => `
                  <div class="abundance-legend-item">
                    <i style="--legend-color:${escapeHtml(segment.color || '#526a86')}"></i>
                    <span class="abundance-legend-name">${escapeHtml(segment.name)}</span>
                    <strong>${escapeHtml(formatRate(segment.ratio))}</strong>
                  </div>
                `).join("")}
              </div>
            </article>
          `;
        }).join("")}
      </div>
    `;
  };

  if (topnSelect && topnSelect.dataset.bound !== 'true') {
    topnSelect.dataset.bound = 'true';
    topnSelect.addEventListener('change', () => {
      state.topn = Number(topnSelect.value || 10);
      render();
    });
  }

  render();
}

function renderMlstSection(section) {
  const tableContainer = document.getElementById("mlst-table");
  const detailContainer = document.getElementById("mlst-gene-show");
  const summaryContainer = document.getElementById("mlst-knowledge-summary-panel");
  const sectionNode = document.getElementById("section-mlst");
  const navLinks = Array.from(document.querySelectorAll('.report-nav-link[href="#section-mlst"]'));
  if (!tableContainer || !detailContainer || !summaryContainer) return;

  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  const geneShowMap = section?.gene_show_map || {};
  const cardTitle = String(section?.title || "MLST 分析结果").trim();
  const tagLabel = String(section?.tag_label || "总结提示").trim();
  const emptyMessage = String(section?.empty_message || "未检出可用于知识库判读的 MLST 结果。").trim();
  const detailEmptyMessage = String(section?.detail_empty_message || "未检出 MLST 结果文件或 host gene 比对文件。").trim();
  const genericDetailNote = String(section?.generic_detail_note || "当前没有额外的基因比对详情。").trim();
  const hostGeneIdIndex = columns.indexOf("Host Gene ID");
  const hostGeneDisplayIndex = columns.indexOf("Host Gene 展示");
  const defaultGene = section?.default_gene || (hostGeneIdIndex >= 0 ? (rows[0]?.[hostGeneIdIndex] ?? "") : "");
  const setMlstVisibility = (visible) => {
    if (sectionNode) sectionNode.classList.toggle("hidden", !visible);
    navLinks.forEach((node) => {
      if (!(node instanceof HTMLElement)) return;
      node.classList.toggle("hidden", !visible);
    });
  };

  if (!rows.length || !columns.length) {
    setMlstVisibility(false);
    buildTableCard("mlst-table", cardTitle, [], []);
    summaryContainer.classList.add("empty-box");
    summaryContainer.innerHTML = `<p>${escapeHtml(emptyMessage)}</p>`;
    detailContainer.innerHTML = `
      <div class="empty-box">
        <p>${escapeHtml(detailEmptyMessage)}</p>
      </div>
    `;
    return;
  }

  setMlstVisibility(true);

  const knowledgeSummary = section?.knowledge_summary || {};
  const summaryHeadline = String(knowledgeSummary?.headline || "").trim();
  const summaryItems = Array.isArray(knowledgeSummary?.items) ? knowledgeSummary.items : [];
  if (summaryHeadline || summaryItems.length) {
    summaryContainer.classList.remove("empty-box");
    summaryContainer.innerHTML = `
      <div class="card-head">
        <h3>${escapeHtml(cardTitle)}</h3>
        <span class="card-tag">${escapeHtml(tagLabel)}</span>
      </div>
      ${summaryHeadline ? `<p class="mlst-knowledge-summary-lead">${escapeHtml(summaryHeadline)}</p>` : ""}
      ${summaryItems.length ? `
        <ul class="report-bullet-list">
          ${summaryItems.map((item) => {
            const fragments = [];
            if (item?.lineage_text && item.lineage_text !== "-") fragments.push(`克隆复合群/Lineage：${item.lineage_text}`);
            if (Array.isArray(item?.virulence) && item.virulence.length) fragments.push(`毒力提示：${item.virulence.join("；")}`);
            if (Array.isArray(item?.resistance) && item.resistance.length) fragments.push(`耐药提示：${item.resistance.join("；")}`);
            if (Array.isArray(item?.regional) && item.regional.length) fragments.push(`地域分布：${item.regional.join("；")}`);
            if (item?.interpretation) fragments.push(`判读提示：${item.interpretation}`);
            return `<li>${escapeHtml(fragments.join("。"))}</li>`;
          }).join("")}
        </ul>
      ` : ""}
    `;
  } else {
    summaryContainer.classList.add("empty-box");
    summaryContainer.innerHTML = `<p>${escapeHtml(emptyMessage)}</p>`;
  }

  const displayColumns = columns.filter((column) => column !== "Host Gene ID" && column !== "Host Gene 展示");
  tableContainer.dataset.exportTitle = cardTitle;
  tableContainer.innerHTML = `
    ${renderTableExportToolbar()}
    <div class="table-frame">
      <table class="report-table mlst-report-table">
        <thead><tr>${displayColumns.map((column) => `<th>${escapeHtml(humanizeReportColumnLabel(column))}</th>`).join("")}</tr></thead>
        <tbody>
          ${rows.map((row) => `
            <tr class="mlst-clickable-row" data-host-gene-id="${escapeHtml(hostGeneIdIndex >= 0 ? (row[hostGeneIdIndex] ?? "") : "")}">
              ${displayColumns.map((column) => {
                const index = columns.indexOf(column);
                const value = row[index] ?? "";
                if (column === "Host Gene") {
                  const displayValue = hostGeneDisplayIndex >= 0 ? (row[hostGeneDisplayIndex] ?? value) : value;
                  return `<td><button class="mlst-hostgene-button" type="button" data-host-gene-id="${escapeHtml(hostGeneIdIndex >= 0 ? (row[hostGeneIdIndex] ?? "") : "")}">${escapeHtml(String(displayValue))}</button></td>`;
                }
                return `<td>${renderTableCellContent(value, column)}</td>`;
              }).join("")}
            </tr>
          `).join("")}
        </tbody>
      </table>
    </div>
  `;
  bindTableExportButtons(
    tableContainer,
    cardTitle,
    displayColumns,
    rows.map((row) => displayColumns.map((column) => row[columns.indexOf(column)] ?? "")),
  );

  if (hostGeneIdIndex < 0 || !Object.keys(geneShowMap || {}).length) {
    detailContainer.classList.remove("empty-box");
    detailContainer.innerHTML = `
      <div class="empty-box">
        <p>${escapeHtml(genericDetailNote)}</p>
      </div>
    `;
    return;
  }

  const renderDetail = (hostGeneId) => {
    const detail = geneShowMap[hostGeneId] || "未检出对应的 host gene 比对详情。";
    detailContainer.classList.remove("empty-box");
    detailContainer.innerHTML = `
      <div class="gene-show-meta">当前 Host Gene：${escapeHtml(hostGeneId || "-")}</div>
      <pre class="gene-show-pre">${escapeHtml(detail)}</pre>
    `;
    tableContainer.querySelectorAll('tbody tr').forEach((rowNode) => {
      rowNode.classList.toggle('mlst-row-active', (rowNode.dataset.hostGeneId || '') === hostGeneId);
    });
  };

  tableContainer.querySelectorAll('[data-host-gene-id]').forEach((node) => {
    node.addEventListener('click', () => renderDetail(node.dataset.hostGeneId || ''));
  });
  renderDetail(defaultGene);
}

function renderNeisseriaAmrSection(amrSection) {
  const summaryContainer = document.getElementById("resistance-neisseria-amr-summary-panel");
  const tableContainer = document.getElementById("resistance-neisseria-amr-table");
  const sectionNode = document.getElementById("section-resistance-mutation");
  const sectionCopy = sectionNode?.querySelector(".section-heading p:last-child");
  if (!summaryContainer || !tableContainer) return;
  const amrHeadline = String(amrSection?.headline || "").trim();
  const amrHighlights = Array.isArray(amrSection?.highlights) ? amrSection.highlights : [];
  const amrInterpretationItems = Array.isArray(amrSection?.interpretation_items) ? amrSection.interpretation_items : [];
  const amrPositiveCount = Number.isFinite(Number(amrSection?.positive_count)) ? Number(amrSection?.positive_count) : 0;
  const amrReviewCount = Number.isFinite(Number(amrSection?.review_count)) ? Number(amrSection?.review_count) : 0;
  const amrSourceLabel = String(amrSection?.source_label || "").trim();
  const cardTitle = String(amrSection?.title || "脑膜炎奈瑟菌耐药突变识别").trim();
  const tagLabel = String(amrSection?.tag_label || "知识库判读").trim();
  const isTbCard = cardTitle.includes("结核分枝杆菌");
  const focusCalls = Array.isArray(amrSection?.focus_calls) ? amrSection.focus_calls : [];
  const otherCalls = Array.isArray(amrSection?.other_calls) ? amrSection.other_calls : [];
  const matchedVariantCount = Number.isFinite(Number(amrSection?.matched_variant_count)) ? Number(amrSection?.matched_variant_count) : 0;
  const tbResistanceGrade = amrSection?.resistance_grade && typeof amrSection.resistance_grade === "object"
    ? amrSection.resistance_grade
    : null;
  const tbLineOrder = ["一线药物", "二线药物", "未分层"];
  const tbEvidenceTone = (value) => {
    const text = String(value || "").trim();
    if (text.includes("Interim")) return "watch";
    if (text.includes("Assoc w R")) return "strong";
    return "neutral";
  };
  const tbMutationPills = (mutations) => {
    const items = Array.isArray(mutations) ? mutations.filter(Boolean).slice(0, 4) : [];
    if (!items.length) return `<span class="tb-amr-mutation-empty">未提供关键突变</span>`;
    return items.map((mutation) => `<span class="tb-amr-mutation-pill">${escapeHtml(String(mutation))}</span>`).join("");
  };
  const tbFocusCards = (() => {
    if (!isTbCard || !focusCalls.length) return "";
    const grouped = new Map();
    focusCalls.forEach((item) => {
      const line = String(item?.drug_line || "未分层").trim() || "未分层";
      if (!grouped.has(line)) grouped.set(line, []);
      grouped.get(line).push(item);
    });
    const groups = tbLineOrder
      .filter((line) => grouped.has(line))
      .map((line) => ({
        line,
        items: (grouped.get(line) || []).slice().sort((left, right) => String(left?.drug || "").localeCompare(String(right?.drug || ""), "zh-CN")),
      }));
    return `
      <div class="tb-amr-focus-panel">
        ${groups.map((group) => `
          <section class="tb-amr-line-group">
            <div class="tb-amr-line-head">
              <h4>${escapeHtml(group.line)}</h4>
              <span>${escapeHtml(String(group.items.length))} 个重点药物</span>
            </div>
            <div class="tb-amr-card-grid">
              ${group.items.map((item) => {
                const evidenceText = String(item?.evidence_grade || item?.verdict || "-").trim() || "-";
                const verdictText = String(item?.verdict || evidenceText || "-").trim() || "-";
                const tone = tbEvidenceTone(evidenceText);
                return `
                  <article class="tb-amr-drug-card">
                    <div class="tb-amr-drug-card-head">
                      <div>
                        <h5>${escapeHtml(String(item?.drug || "-"))}</h5>
                        <p>${escapeHtml(verdictText)}</p>
                      </div>
                      <span class="tb-amr-evidence-badge is-${tone}">${escapeHtml(evidenceText)}</span>
                    </div>
                    <div class="tb-amr-mutation-row">
                      ${tbMutationPills(item?.mutations)}
                    </div>
                  </article>
                `;
              }).join("")}
            </div>
          </section>
        `).join("")}
      </div>
    `;
  })();
  const sectionDescription = cardTitle.includes("结核分枝杆菌")
    ? "当物种鉴定提示为结核分枝杆菌时，这里将展示基于 H37Rv 有参 SNP 与 WHO mutation catalogue 的耐药位点判读。"
    : "当 MLST 或物种鉴定提示为脑膜炎奈瑟菌时，这里将展示基于组装注释结果的耐药突变识别与知识库判读。";
  if (sectionCopy) {
    sectionCopy.textContent = sectionDescription;
  }
  if (amrSection?.status === "ready" || amrHeadline) {
    summaryContainer.hidden = false;
    tableContainer.hidden = false;
    summaryContainer.classList.remove("empty-box");
    summaryContainer.innerHTML = `
      <div class="card-head">
        <h3>${escapeHtml(cardTitle)}</h3>
        <span class="card-tag">${escapeHtml(tagLabel)}</span>
      </div>
      ${isTbCard && tbResistanceGrade?.label ? `
        <div class="taxonomy-summary-grid taxonomy-interpretation-grid" style="margin-bottom:18px;">
          <article class="summary-metric-card">
            <span class="summary-metric-label">结核耐药分级</span>
            <strong>${escapeHtml(String(tbResistanceGrade.label || "-"))}</strong>
            ${tbResistanceGrade?.reason ? `<p>${escapeHtml(String(tbResistanceGrade.reason || ""))}</p>` : ""}
          </article>
        </div>
      ` : ""}
      ${amrHeadline ? `<p class="mlst-knowledge-summary-lead">${escapeHtml(amrHeadline)}</p>` : ""}
      <div class="summary-metric-strip">
        <span class="mini-pill">${escapeHtml(isTbCard ? "重点药物" : "命中位点")} ${escapeHtml(String(amrPositiveCount))}</span>
        <span class="mini-pill">${escapeHtml(isTbCard ? "目录命中总药物" : "待复核")} ${escapeHtml(String(amrReviewCount))}</span>
        ${isTbCard ? `<span class="mini-pill">目录命中总条目 ${escapeHtml(String(matchedVariantCount))}</span>` : ""}
      </div>
      ${amrInterpretationItems.length ? `
        <ul class="report-bullet-list${isTbCard ? " tb-amr-summary-list" : ""}">
          ${amrInterpretationItems.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
        </ul>
      ` : ""}
      ${tbFocusCards}
      ${!isTbCard && amrHighlights.length ? `
        <ul class="report-bullet-list">
          ${amrHighlights.map((item) => `<li>${escapeHtml(item)}</li>`).join("")}
        </ul>
      ` : ""}
      ${amrSourceLabel ? `<p class="section-note">数据库来源：${escapeHtml(amrSourceLabel)}</p>` : ""}
    `;
    buildTableCard(
      "resistance-neisseria-amr-table",
      cardTitle,
      Array.isArray(amrSection?.columns) ? amrSection.columns : [],
      Array.isArray(amrSection?.rows) ? amrSection.rows : [],
    );
  } else {
    summaryContainer.hidden = true;
    tableContainer.hidden = true;
    summaryContainer.innerHTML = "";
    tableContainer.innerHTML = "";
  }
}

function formatGenomeCoordinate(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return String(value ?? "--");
  return numeric.toLocaleString("zh-CN");
}

function buildCoveragePairs(section) {
  const points = Array.isArray(section?.points) ? section.points : [];
  const xValues = Array.isArray(section?.x_values) ? section.x_values : [];
  return points.map((depth, index) => {
    const x = Number(xValues[index] ?? index + 1);
    const y = Number(depth) || 0;
    return { x, y };
  }).filter((item) => Number.isFinite(item.x) && Number.isFinite(item.y));
}

function buildCoverageWindowStats(pairs) {
  if (!Array.isArray(pairs) || !pairs.length) {
    return { mean: 0, max: 0, min: 0 };
  }
  const values = pairs.map((item) => Number(item.y) || 0);
  return {
    mean: values.reduce((sum, value) => sum + value, 0) / values.length,
    max: Math.max(...values),
    min: Math.min(...values),
  };
}

function buildCoverageThresholdFraction(pairs, threshold) {
  if (!Array.isArray(pairs) || !pairs.length || !Number.isFinite(threshold) || threshold <= 0) return null;
  const hitCount = pairs.filter((item) => (Number(item.y) || 0) >= threshold).length;
  return pairs.length ? hitCount / pairs.length : null;
}

function layoutNcovGenomeFeatures(features, totalBases) {
  const lanes = [];
  const laneGap = Math.max(110, totalBases * 0.01);
  return (Array.isArray(features) ? features : []).map((feature) => {
    const start = Number(feature?.start) || 1;
    const end = Number(feature?.end) || start;
    let laneIndex = 0;
    while (laneIndex < lanes.length && start <= lanes[laneIndex] + laneGap) {
      laneIndex += 1;
    }
    lanes[laneIndex] = end;
    return {
      ...feature,
      laneIndex,
    };
  });
}

function buildNcovCoveragePlotSvg(pairs, options = {}) {
  if (!Array.isArray(pairs) || !pairs.length) return "";
  const width = options.width || 1140;
  const height = options.height || 360;
  const padLeft = 80;
  const padRight = 30;
  const padTop = 22;
  const padBottom = 52;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const domainStart = Number(options.domainStart) || pairs[0].x;
  const domainEnd = Number(options.domainEnd) || pairs[pairs.length - 1].x;
  const maxDepth = Math.max(...pairs.map((item) => item.y), 1);
  const roundedMax = Math.max(10, Number(options.maxDepth) || Math.ceil(maxDepth * 1.08));
  const scaleX = (value) => padLeft + ((value - domainStart) / Math.max(domainEnd - domainStart, 1)) * innerWidth;
  const scaleY = (value) => padTop + innerHeight - (value / Math.max(roundedMax, 1)) * innerHeight;
  const linePoints = pairs.map((item) => `${scaleX(item.x)},${scaleY(item.y)}`).join(" ");
  const areaPoints = [
    `${scaleX(pairs[0].x)},${height - padBottom}`,
    ...pairs.map((item) => `${scaleX(item.x)},${scaleY(item.y)}`),
    `${scaleX(pairs[pairs.length - 1].x)},${height - padBottom}`,
  ].join(" ");
  const yGrid = Array.from({ length: 5 }, (_, index) => {
    const ratio = index / 4;
    const value = roundedMax * ratio;
    const y = scaleY(value);
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 12}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const tickCount = 5;
  const xTicks = Array.from({ length: tickCount }, (_, index) => {
    const ratio = tickCount === 1 ? 0 : index / (tickCount - 1);
    const value = Math.round(domainStart + (domainEnd - domainStart) * ratio);
    const x = scaleX(value);
    return `
      <line class="chart-axis-tick" x1="${x}" y1="${height - padBottom}" x2="${x}" y2="${height - padBottom + 6}"></line>
      <text class="chart-axis-label x-axis" x="${x}" y="${height - 24}">${escapeHtml(formatGenomeCoordinate(value))}</text>
    `;
  }).join("");
  return `
    <div class="chart-canvas ncov-depth-canvas" style="--chart-height:${height}px">
      <svg class="sparkline-svg ncov-depth-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="新冠基因组测序深度图">
        ${yGrid}
        <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
        <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
        <polygon class="ncov-depth-area" points="${areaPoints}"></polygon>
        <polyline class="chart-series-line" fill="none" stroke="#2f6fd6" stroke-width="2.6" points="${linePoints}"></polyline>
        ${xTicks}
        <text class="chart-axis-title" x="${width / 2}" y="${height - 4}">基因组坐标（nt）</text>
        <text class="chart-axis-title chart-axis-title-y" x="26" y="${height / 2}">测序深度</text>
      </svg>
    </div>
  `;
}

function renderAnnotatedNcovCoverage(container, section) {
  const pairs = buildCoveragePairs(section);
  const features = Array.isArray(section?.genome_features) ? section.genome_features : [];
  const totalBases = Number(section?.total_bases) || Math.max(...pairs.map((item) => item.x), 1);
  const overviewFeatures = features.filter((item) => ["gene", "CDS", "five_prime_UTR", "three_prime_UTR"].includes(String(item?.feature_type || "").trim()));
  const fullRange = { start: 1, end: totalBases, label: "全基因组" };
  let currentRange = fullRange;
  let depthMode = "raw";

  const getPairsInRange = (range) => {
    const filtered = pairs.filter((item) => item.x >= range.start && item.x <= range.end);
    return filtered.length ? filtered : pairs;
  };

  const getRangeWithPadding = (feature) => {
    const start = Number(feature?.start) || 1;
    const end = Number(feature?.end) || start;
    const padding = Math.max(80, Math.round((end - start + 1) * 0.08));
    return {
      start: Math.max(1, start - padding),
      end: Math.min(totalBases, end + padding),
      label: String(feature?.label || "局部区域"),
      feature,
    };
  };

  const getDepthModeConfig = () => {
    if (depthMode === "10x") return { threshold: 10, maxDepth: 10, label: "10x 深度" };
    if (depthMode === "100x") return { threshold: 100, maxDepth: 100, label: "100x 深度" };
    return { threshold: null, maxDepth: null, label: "原始深度" };
  };

  const render = () => {
    const rawVisiblePairs = getPairsInRange(currentRange);
    const depthModeConfig = getDepthModeConfig();
    const visiblePairs = rawVisiblePairs.map((item) => ({
      ...item,
      y: depthModeConfig.maxDepth != null ? Math.min(item.y, depthModeConfig.maxDepth) : item.y,
    }));
    const stats = buildCoverageWindowStats(rawVisiblePairs);
    const activeFeature = currentRange.feature || null;
    const thresholdFraction = buildCoverageThresholdFraction(rawVisiblePairs, depthModeConfig.threshold);
    const windowMetrics = [
      { label: "当前范围", value: `${formatGenomeCoordinate(currentRange.start)} - ${formatGenomeCoordinate(currentRange.end)}` },
      { label: "局部平均深度", value: formatChartValue(stats.mean) },
      { label: "局部最大深度", value: formatChartValue(stats.max) },
      {
        label: depthModeConfig.threshold ? `达到 ${depthModeConfig.threshold}x 比例` : "参考注释",
        value: depthModeConfig.threshold ? `${((thresholdFraction || 0) * 100).toFixed(1)}%` : (section?.reference_name || "NC_045512.2"),
      },
    ];
    const metricMarkup = windowMetrics.map((item) => `
      <div class="coverage-focus-metric">
        <span>${escapeHtml(item.label)}</span>
        <strong>${escapeHtml(String(item.value))}</strong>
      </div>
    `).join("");
    const laidOutFeatures = layoutNcovGenomeFeatures(overviewFeatures, totalBases);
    const laneCount = Math.max(...laidOutFeatures.map((item) => Number(item.laneIndex) || 0), 0) + 1;
    const laneHeight = 42;
    const trackHeight = laneCount * laneHeight + 26;
    const geneTrackMarkup = laidOutFeatures.map((feature) => {
      const start = Number(feature?.start) || 1;
      const end = Number(feature?.end) || start;
      const left = ((start - 1) / Math.max(totalBases, 1)) * 100;
      const width = ((end - start + 1) / Math.max(totalBases, 1)) * 100;
      const label = String(feature?.label || feature?.feature_type || "feature");
      const active = activeFeature && activeFeature.label === feature.label && activeFeature.start === feature.start && activeFeature.end === feature.end;
      const featureKind = String(feature?.category || feature?.feature_type || "gene").toLowerCase();
      const title = `${label} ${formatGenomeCoordinate(start)}-${formatGenomeCoordinate(end)}`;
      const fillColor = getSarsCov2GeneColor(label, feature?.feature_type);
      const textColor = getGenomeFeatureTextColor(label);
      const strandClass = String(feature?.strand || "+").trim() === "-" ? " is-reverse" : " is-forward";
      const top = 18 + (Number(feature?.laneIndex) || 0) * laneHeight;
      const compactClass = width < 4.2 ? " is-compact" : "";
      return `
        <button
          type="button"
          class="ncov-gene-block is-${escapeHtml(featureKind)}${strandClass}${compactClass}${active ? " is-active" : ""}"
          data-gene-start="${start}"
          data-gene-end="${end}"
          data-gene-label="${escapeHtml(label)}"
          style="left:${left}%;top:${top}px;width:${Math.max(width, 0.55)}%;--gene-fill:${escapeHtml(fillColor)};--gene-text:${escapeHtml(textColor)};"
          title="${escapeHtml(title)}"
        >
          <span>${escapeHtml(label)}</span>
        </button>
      `;
    }).join("");
    const overviewWindowLeft = ((currentRange.start - 1) / Math.max(totalBases, 1)) * 100;
    const overviewWindowWidth = ((currentRange.end - currentRange.start + 1) / Math.max(totalBases, 1)) * 100;
    container.classList.remove("empty-box");
    container.innerHTML = `
      <div class="coverage-summary-bar">
        <div class="coverage-summary-item">
          <span>总位点</span>
          <strong>${escapeHtml(String(section.total_bases ?? "-"))}</strong>
        </div>
        <div class="coverage-summary-item">
          <span>平均深度</span>
          <strong>${escapeHtml(String(section.mean_depth ?? "-"))}</strong>
        </div>
        <div class="coverage-summary-item">
          <span>最大深度</span>
          <strong>${escapeHtml(String(section.max_depth ?? "-"))}</strong>
        </div>
        <div class="coverage-summary-item">
          <span>Contig 数</span>
          <strong>${escapeHtml(String(section.contig_count ?? "-"))}</strong>
        </div>
      </div>
      <div class="mini-chart-card coverage-chart-card ncov-coverage-card">
        <div class="ncov-coverage-head">
          <div>
            <span class="mini-chart-title">${escapeHtml(section.label || "基因组覆盖度")}</span>
            <h3>${escapeHtml(currentRange.label || "全基因组")} 深度视图</h3>
            <p>下方注释轨使用 ${escapeHtml(String(section?.annotation_label || section?.annotation_source || "参考基因组 GFF"))}。点击不同基因区域可快速放大到对应区段。</p>
          </div>
        </div>
        <div class="ncov-coverage-focus-bar">
          ${metricMarkup}
        </div>
        <div class="ncov-coverage-overview">
          <div class="ncov-coverage-overview-track"></div>
          <div class="ncov-coverage-overview-window" style="left:${overviewWindowLeft}%;width:${Math.max(overviewWindowWidth, 0.8)}%;"></div>
        </div>
        <div class="ncov-plot-frame">
          <div class="ncov-plot-toolbar">
            <div class="subreport-tabs ncov-depth-tabs" role="tablist" aria-label="覆盖度模式切换">
              <button class="subreport-tab-button${depthMode === "raw" ? " active" : ""}" type="button" data-ncov-depth-mode="raw">原始深度</button>
              <button class="subreport-tab-button${depthMode === "10x" ? " active" : ""}" type="button" data-ncov-depth-mode="10x">10x</button>
              <button class="subreport-tab-button${depthMode === "100x" ? " active" : ""}" type="button" data-ncov-depth-mode="100x">100x</button>
            </div>
            <button type="button" class="table-export-button ncov-coverage-reset"${currentRange.start === 1 && currentRange.end === totalBases ? " disabled" : ""}>返回全长</button>
          </div>
          ${buildNcovCoveragePlotSvg(visiblePairs, { domainStart: currentRange.start, domainEnd: currentRange.end, maxDepth: depthModeConfig.maxDepth })}
          <div class="ncov-inline-annotation">
            <div class="ncov-annotation-scale">
              <span>1</span>
              <span>${escapeHtml(formatGenomeCoordinate(Math.round(totalBases / 2)))}</span>
              <span>${escapeHtml(formatGenomeCoordinate(totalBases))}</span>
            </div>
            <div class="ncov-gene-track" style="height:${trackHeight}px">
              <div class="ncov-gene-track-axis"></div>
              ${geneTrackMarkup}
            </div>
          </div>
        </div>
      </div>
    `;
    const resetButton = container.querySelector(".ncov-coverage-reset");
    if (resetButton) {
      resetButton.addEventListener("click", () => {
        currentRange = fullRange;
        render();
      });
    }
    container.querySelectorAll("[data-ncov-depth-mode]").forEach((button) => {
      button.addEventListener("click", () => {
        depthMode = String(button.getAttribute("data-ncov-depth-mode") || "raw");
        render();
      });
    });
    container.querySelectorAll(".ncov-gene-block").forEach((button) => {
      button.addEventListener("click", () => {
        const start = Number(button.getAttribute("data-gene-start") || 1);
        const end = Number(button.getAttribute("data-gene-end") || start);
        const label = button.getAttribute("data-gene-label") || "区域";
        currentRange = getRangeWithPadding({ start, end, label });
        render();
      });
    });
  };

  render();
}

function renderAssemblyCoverage(section) {
  const container = document.getElementById("assembly-coverage-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出覆盖度结果文件。</p>
      </div>
    `;
    return;
  }
  if (String(section?.view_mode || "").trim() === "ncov_annotated" && Array.isArray(section?.genome_features) && section.genome_features.length) {
    renderAnnotatedNcovCoverage(container, section);
    return;
  }
  const isInfluenza = isInfluenzaTypingReport(currentReportData);
  const isBandavirus = isBandavirusTypingReport(currentReportData);
  if ((isInfluenza || isBandavirus) && Array.isArray(section?.segments) && section.segments.length > 1) {
    container.classList.remove("empty-box");
    container.innerHTML = `
      <div class="coverage-summary-bar">
        <div class="coverage-summary-item">
          <span>${escapeHtml(isBandavirus ? "片段数" : "Segment 数")}</span>
          <strong>${escapeHtml(String(section.segments.length))}</strong>
        </div>
        <div class="coverage-summary-item">
          <span>总位点</span>
          <strong>${escapeHtml(String(section.total_bases ?? "-"))}</strong>
        </div>
        <div class="coverage-summary-item">
          <span>平均深度</span>
          <strong>${escapeHtml(String(section.mean_depth ?? "-"))}</strong>
        </div>
        <div class="coverage-summary-item">
          <span>整体覆盖度</span>
          <strong>${escapeHtml(formatRate(Number(section.coverage_fraction || 0)))}</strong>
        </div>
      </div>
      <div class="influenza-segment-coverage-grid">
        ${section.segments.map((segment) => {
          const segmentPoints = Array.isArray(segment?.points) ? segment.points : [];
          const bandavirusMeta = isBandavirus
            ? [
              `A_F ${String(segment?.af_group || "--").trim() || "--"}`,
              `CJ ${String(segment?.cj_group || "--").trim() || "--"}`,
              String(segment?.accession || "").trim() || "--",
            ].join(" · ")
            : "";
          const chartSvg = createSeriesSvg(
            [{ label: "测序深度", color: "#4e6177", values: segmentPoints }],
            {
              label: String(segment?.label || segment?.name || "segment"),
              width: 560,
              height: 240,
              padX: 40,
              padTop: 12,
              padBottom: 38,
              xLabel: "segment 位置",
              yLabel: "深度",
              xTicks: Array.isArray(segment?.x_ticks) ? segment.x_ticks : [1, Math.max(1, Math.round(segmentPoints.length / 2)), segmentPoints.length],
              xValues: Array.isArray(segment?.x_values) ? segment.x_values : [],
            },
          );
          return `
            <article class="mini-chart-card influenza-segment-coverage-card">
              <div class="influenza-segment-coverage-head">
                <div>
                  <span class="mini-chart-title">${escapeHtml(String(segment?.label || segment?.name || "segment"))}</span>
                  <p class="empty-copy">覆盖度 ${escapeHtml(formatRate(Number(segment?.coverage_fraction || 0)))} · 10x ${escapeHtml(formatRate(Number(segment?.coverage_10x_fraction || 0)))} · 100x ${escapeHtml(formatRate(Number(segment?.coverage_100x_fraction || 0)))}</p>
                  ${bandavirusMeta ? `<p class="empty-copy">${escapeHtml(bandavirusMeta)}</p>` : ""}
                </div>
                <span class="card-tag">${escapeHtml(`${segment?.total_bases ?? "-"} bp`)}</span>
              </div>
              ${chartSvg}
            </article>
          `;
        }).join("")}
      </div>
    `;
    return;
  }
  const chartSvg = createSeriesSvg(
    [{ label: "测序深度", color: "#4e6177", values: points }],
    {
      label: section.label || "基因组覆盖度",
      width: 1120,
      height: 420,
      padX: 44,
      padTop: 10,
      padBottom: 42,
      xLabel: section.x_label || "基因组位置",
      yLabel: section.y_label || "测序深度",
      xTicks: section.x_ticks || [1, Math.max(1, Math.round(points.length / 2)), points.length],
      xValues: section.x_values || [],
    },
  );
  container.classList.remove("empty-box");
  container.innerHTML = `
    <div class="coverage-summary-bar">
      <div class="coverage-summary-item">
        <span>总位点</span>
        <strong>${escapeHtml(String(section.total_bases ?? "-"))}</strong>
      </div>
      <div class="coverage-summary-item">
        <span>平均深度</span>
        <strong>${escapeHtml(String(section.mean_depth ?? "-"))}</strong>
      </div>
      <div class="coverage-summary-item">
        <span>最大深度</span>
        <strong>${escapeHtml(String(section.max_depth ?? "-"))}</strong>
      </div>
      <div class="coverage-summary-item">
        <span>Contig 数</span>
        <strong>${escapeHtml(String(section.contig_count ?? "-"))}</strong>
      </div>
    </div>
    <div class="mini-chart-card coverage-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || "基因组覆盖度")}</span>
      ${buildChartInsight(summarizeCoverage(section))}
      ${chartSvg}
    </div>
  `;
}


function renderBarSvg(values, options = {}) {
  const width = options.width || 1100;
  const height = options.height || 360;
  const padX = options.padX ?? 60;
  const padTop = options.padTop ?? 32;
  const padBottom = options.padBottom ?? 88;
  const innerWidth = width - padX * 2;
  const innerHeight = height - padTop - padBottom;
  const maxValue = Math.max(...values.map((value) => Number(value) || 0), 1);
  const yTicks = 4;
  const valueLabelStep = values.length <= 12 ? 1 : Math.ceil(values.length / 8);
  const grid = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = maxValue * ratio;
    return `
      <line class="chart-grid-line" x1="${padX}" y1="${y}" x2="${width - padX}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padX - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const barGap = 8;
  const barWidth = Math.max((innerWidth - Math.max(values.length - 1, 0) * barGap) / Math.max(values.length, 1), 10);
  const bars = values.map((value, index) => {
    const numeric = Number(value) || 0;
    const barHeight = (numeric / Math.max(maxValue, 1)) * innerHeight;
    const x = padX + index * (barWidth + barGap);
    const y = padTop + innerHeight - barHeight;
    const showValueLabel = index % valueLabelStep === 0 || index === values.length - 1;
    return `
      <rect class="chart-bar" x="${x}" y="${y}" width="${barWidth}" height="${barHeight}" rx="6" ry="6">
        <title>${escapeHtml(`${options.xValues?.[index] ?? index + 1}: ${numeric}`)}</title>
      </rect>
      ${showValueLabel ? `<text class="chart-bar-value" x="${x + barWidth / 2}" y="${Math.max(y - 10, padTop - 8)}">${escapeHtml(String(numeric))}</text>` : ""}
      <text class="chart-axis-label x-axis chart-rotated-label" x="${x + barWidth / 2}" y="${height - 30}">${escapeHtml(String(options.xValues?.[index] ?? index + 1))}</text>
    `;
  }).join("");
  return `
    <div class="mini-chart-card coverage-chart-card">
      <span class="mini-chart-title">${escapeHtml(options.label || "柱状分布图")}</span>
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg bar-chart-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(options.label || "柱状分布图")}">
          ${grid}
          <line class="chart-axis-line" x1="${padX}" y1="${padTop}" x2="${padX}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padX}" y1="${height - padBottom}" x2="${width - padX}" y2="${height - padBottom}"></line>
          ${bars}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 8}">${escapeHtml(options.xLabel || "分类")}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(options.yLabel || "数量")}</text>
        </svg>
      </div>
    </div>
  `;
}


function renderGeneLengthDistribution(section) {
  const container = document.getElementById("gene-length-distribution-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出基因长度统计文件。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  container.innerHTML = `
    ${buildChartInsight(summarizeGeneLengthDistribution(section))}
    ${renderBarSvg(points, {
    label: section.label || "基因长度与数量分布",
    width: 1120,
    height: 420,
    padX: 48,
    padTop: 24,
    padBottom: 96,
    xLabel: section.x_label || "基因长度范围",
    yLabel: section.y_label || "Gene数量",
    xValues: section.x_values || [],
  })}
  `;
}


function renderContigDepthRelationship(section) {
  const container = document.getElementById("contig-depth-relationship-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的序列类型与平均深度关系数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 980;
  const height = 320;
  const padLeft = 84;
  const padRight = 38;
  const padTop = 26;
  const padBottom = 58;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const groups = ["基因组", "质粒"];
  const maxDepth = Math.max(...points.map((point) => Number(point.depth) || 0), 1);
  const roundedMax = Math.ceil(maxDepth * 1.08);
  const xPositions = { "基因组": padLeft + innerWidth * 0.28, "质粒": padLeft + innerWidth * 0.72 };
  const yTicks = 4;
  const grid = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = roundedMax * ratio;
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const dots = points.map((point, index) => {
    const baseX = xPositions[point.type] || xPositions["基因组"];
    const jitter = ((index % 11) - 5) * 7;
    const depth = Number(point.depth) || 0;
    const y = padTop + innerHeight - (depth / Math.max(roundedMax, 1)) * innerHeight;
    const fill = point.type === "质粒" ? "#a88262" : "#4e6177";
    return `<circle cx="${baseX + jitter}" cy="${y}" r="5.5" fill="${fill}" fill-opacity="0.82"><title>${escapeHtml(`${point.name} | ${point.type} | 平均深度 ${depth}`)}</title></circle>`;
  }).join("");
  const labels = groups.map((group) => `<text class="chart-axis-label x-axis" x="${xPositions[group]}" y="${height - 28}">${group}</text>`).join("");
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || "基因组/质粒与平均深度关系图")}</span>
      ${buildChartInsight(summarizeContigDepth(section))}
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(section.label || "基因组/质粒与平均深度关系图")}">
          ${grid}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
          ${labels}
          ${dots}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 6}">${escapeHtml(section.x_label || "序列类型")}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(section.y_label || "平均深度")}</text>
        </svg>
      </div>
    </div>
  `;
}

function summarizeContigLengthDepthScatter(section) {
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    return "未检出可用于长度与平均深度联合展示的 Contig 数据。";
  }
  const genome = points.filter((point) => point.type === "基因组");
  const plasmid = points.filter((point) => point.type === "质粒");
  const maxLengthPoint = points.reduce((best, current) => ((Number(current.length) || 0) > (Number(best.length) || 0) ? current : best), points[0]);
  const maxDepthPoint = points.reduce((best, current) => ((Number(current.depth) || 0) > (Number(best.depth) || 0) ? current : best), points[0]);
  const plasmidMeanDepth = plasmid.length
    ? plasmid.reduce((sum, item) => sum + (Number(item.depth) || 0), 0) / plasmid.length
    : 0;
  const genomeMeanDepth = genome.length
    ? genome.reduce((sum, item) => sum + (Number(item.depth) || 0), 0) / genome.length
    : 0;
  if (plasmid.length && genome.length && plasmidMeanDepth > genomeMeanDepth * 1.35) {
    return `质粒序列整体平均深度高于基因组背景，提示部分质粒可能具有更高拷贝特征；最长Contig为 ${maxLengthPoint?.name || "-"}。`;
  }
  return `Contig 长度与平均深度分布整体可见，最长Contig为 ${maxLengthPoint?.name || "-"}，最高深度序列为 ${maxDepthPoint?.name || "-"}。`;
}

function renderContigLengthDepthScatter(section) {
  const container = document.getElementById("contig-length-depth-scatter-chart");
  if (!container) return;
  const points = Array.isArray(section?.points) ? section.points : [];
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的 Contig 长度与平均测序深度数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 920;
  const height = 320;
  const padLeft = 86;
  const padRight = 38;
  const padTop = 24;
  const padBottom = 66;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const maxLength = Math.max(...points.map((point) => Number(point.length) || 0), 1);
  const maxDepth = Math.max(...points.map((point) => Number(point.depth) || 0), 1);
  const roundedLength = Math.ceil(maxLength * 1.04);
  const roundedDepth = Math.ceil(maxDepth * 1.08);
  const xTicks = 4;
  const yTicks = 4;
  const gridX = Array.from({ length: xTicks + 1 }, (_, index) => {
    const ratio = index / xTicks;
    const x = padLeft + innerWidth * ratio;
    const value = roundedLength * ratio;
    return `
      <line class="chart-grid-line" x1="${x}" y1="${padTop}" x2="${x}" y2="${height - padBottom}"></line>
      <text class="chart-axis-label x-axis" x="${x}" y="${height - 28}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const gridY = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const y = padTop + innerHeight - innerHeight * ratio;
    const value = roundedDepth * ratio;
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const dots = points.map((point) => {
    const length = Number(point.length) || 0;
    const depth = Number(point.depth) || 0;
    const x = padLeft + (length / Math.max(roundedLength, 1)) * innerWidth;
    const y = padTop + innerHeight - (depth / Math.max(roundedDepth, 1)) * innerHeight;
    const fill = point.type === "质粒" ? "#a88262" : "#4e6177";
    return `<circle cx="${x}" cy="${y}" r="5" fill="${fill}" fill-opacity="0.82"><title>${escapeHtml(`${point.name} | ${point.type} | 长度 ${formatChartValue(length)} bp | 平均深度 ${depth}`)}</title></circle>`;
  }).join("");
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || "Contig长度与平均测序深度散点图")}</span>
      ${buildChartInsight(summarizeContigLengthDepthScatter(section))}
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(section.label || "Contig长度与平均测序深度散点图")}">
          ${gridX}
          ${gridY}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
          ${dots}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 6}">${escapeHtml(section.x_label || "Contig长度(bp)")}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(section.y_label || "平均测序深度")}</text>
        </svg>
      </div>
    </div>
  `;
}

function renderCategoryGeneRelationship(containerId, section) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const leftNodes = Array.isArray(section?.nodes_left) ? section.nodes_left : [];
  const rightNodes = Array.isArray(section?.nodes_right) ? section.nodes_right : [];
  const links = Array.isArray(section?.links) ? section.links : [];
  if (!leftNodes.length || !rightNodes.length || !links.length) {
    container.innerHTML = `
      <div class="empty-box coverage-empty">
        <p>未检出可展示的关系数据。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 1040;
  const height = Math.max(300, Math.max(leftNodes.length, rightNodes.length) * 34 + 70);
  const leftX = 176;
  const rightX = width - 176;
  const leftYStep = leftNodes.length > 1 ? (height - 80) / (leftNodes.length - 1) : 0;
  const rightYStep = rightNodes.length > 1 ? (height - 80) / (rightNodes.length - 1) : 0;
  const leftPos = Object.fromEntries(leftNodes.map((node, index) => [node.name, 40 + index * leftYStep]));
  const rightPos = Object.fromEntries(rightNodes.map((node, index) => [node.name, 40 + index * rightYStep]));
  const maxLink = Math.max(...links.map((link) => Number(link.value) || 0), 1);
  const linkSvg = links.map((link) => {
    const y1 = leftPos[link.source];
    const y2 = rightPos[link.target];
    if (y1 === undefined || y2 === undefined) return '';
    const strokeWidth = 1.6 + ((Number(link.value) || 0) / maxLink) * 5.4;
    return `<path d="M ${leftX} ${y1} C ${leftX + 160} ${y1}, ${rightX - 160} ${y2}, ${rightX} ${y2}" class="relation-link" data-source="${escapeHtml(link.source)}" data-target="${escapeHtml(link.target)}" style="stroke-width:${strokeWidth}px"><title>${escapeHtml(`${link.source} -> ${link.target}: ${link.value}`)}</title></path>`;
  }).join('');
  const leftSvg = leftNodes.map((node) => `
    <g class="relation-node-group" data-relation-side="left" data-relation-name="${escapeHtml(node.name)}" transform="translate(0 ${leftPos[node.name]})">
      <text class="relation-label relation-label-left" x="154" y="4">${escapeHtml(node.name)}</text>
      <rect class="relation-node relation-node-left" x="160" y="-10" width="16" height="20" rx="6"></rect>
      <text class="relation-value relation-value-left" x="150" y="-14">${escapeHtml(String(node.value))}</text>
    </g>
  `).join('');
  const rightSvg = rightNodes.map((node) => `
    <g class="relation-node-group" data-relation-side="right" data-relation-name="${escapeHtml(node.name)}" transform="translate(0 ${rightPos[node.name]})">
      <rect class="relation-node relation-node-right" x="${rightX}" y="-10" width="16" height="20" rx="6"></rect>
      <text class="relation-label relation-label-right" x="${rightX + 22}" y="4">${escapeHtml(node.name)}</text>
      <text class="relation-value relation-value-right" x="${rightX + 22}" y="-14">${escapeHtml(String(node.value))}</text>
    </g>
  `).join('');
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card">
      <span class="mini-chart-title">${escapeHtml(section.label || '关系图')}</span>
      ${buildChartInsight(summarizeRelationship(section))}
      <div class="relationship-chart-head">
        <span>${escapeHtml(section.left_label || '左侧分类')}</span>
        <span>${escapeHtml(section.right_label || '右侧基因')}</span>
      </div>
      <div class="chart-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg relation-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(section.label || '关系图')}">
          ${linkSvg}
          ${leftSvg}
          ${rightSvg}
        </svg>
      </div>
    </div>
  `;
  const relationLinks = Array.from(container.querySelectorAll(".relation-link"));
  const relationNodes = Array.from(container.querySelectorAll(".relation-node-group"));
  const clearHighlight = () => {
    relationLinks.forEach((link) => link.classList.remove("is-dimmed", "is-highlighted"));
    relationNodes.forEach((node) => node.classList.remove("is-dimmed", "is-highlighted"));
  };
  const applyHighlight = (side, name) => {
    relationLinks.forEach((link) => {
      const matched = side === "left"
        ? link.dataset.source === name
        : link.dataset.target === name;
      link.classList.toggle("is-highlighted", matched);
      link.classList.toggle("is-dimmed", !matched);
    });
    relationNodes.forEach((node) => {
      const matched = node.dataset.relationName === name
        || relationLinks.some((link) => link.classList.contains("is-highlighted") && (link.dataset.source === node.dataset.relationName || link.dataset.target === node.dataset.relationName));
      node.classList.toggle("is-highlighted", matched);
      node.classList.toggle("is-dimmed", !matched);
    });
  };
  relationNodes.forEach((node) => {
    node.addEventListener("mouseenter", () => applyHighlight(node.dataset.relationSide || "", node.dataset.relationName || ""));
    node.addEventListener("mouseleave", clearHighlight);
  });
}

function renderRawQc(sections) {
  const left = sections?.raw_qc?.paired_end?.left || {};
  const right = sections?.raw_qc?.paired_end?.right || {};
  const fastp = sections?.raw_qc?.fastp || {};
  const rawQc = sections?.raw_qc || {};
  const paired = isPairedEndFastp(fastp, rawQc);
  const rawQcGrid = document.querySelector("#section-raw-qc .two-column");
  const fastpGrid = document.querySelector("#section-fastp .two-column");
  rawQcGrid?.classList.toggle("single-end-layout", !paired);
  fastpGrid?.classList.toggle("single-end-layout", !paired);
  renderStatsWithCharts("raw-qc-left", left, paired ? "R1" : "SE");
  const leftCard = document.getElementById("raw-qc-left")?.closest(".result-card");
  const rightCard = document.getElementById("raw-qc-right")?.closest(".result-card");
  if (leftCard) {
    const title = leftCard.querySelector(".card-head h3");
    const badge = leftCard.querySelector(".card-head .card-tag");
    if (title) title.textContent = paired ? "左端测序数据" : "单端测序数据";
    if (badge) badge.textContent = paired ? "R1" : "SE";
  }
  if (rightCard) {
    rightCard.style.display = paired ? "" : "none";
  }
  if (paired) {
    renderStatsWithCharts("raw-qc-right", right, "R2");
  }
  const summaryItems = [
    { label: "测序模式", value: String(fastp.sequencing || "-") },
    { label: "过滤后 Reads", value: String(fastp.filtering_result?.passed_filter_reads ?? "-") },
    { label: "过短 Reads", value: String(fastp.filtering_result?.too_short_reads ?? "-") },
    { label: "N 过多 Reads", value: String(fastp.filtering_result?.too_many_N_reads ?? "-") },
    { label: "重复率", value: formatRate(fastp.duplication_rate) },
  ];
  const summaryBox = document.getElementById("fastp-summary");
  if (summaryBox) {
    summaryBox.classList.remove("empty-box");
    const [primaryItem, ...secondaryItems] = summaryItems;
    const passed = Number(fastp.filtering_result?.passed_filter_reads || 0);
    const tooShort = Number(fastp.filtering_result?.too_short_reads || 0);
    const tooManyN = Number(fastp.filtering_result?.too_many_N_reads || 0);
    const duplication = Number(fastp.duplication_rate || 0);
    const totalObserved = passed + tooShort + tooManyN;
    const lossRate = totalObserved ? (tooShort + tooManyN) / totalObserved : 0;
    const riskState = duplication >= 0.5 || lossRate >= 0.2 ? 'danger' : duplication >= 0.3 || lossRate >= 0.1 ? 'warning' : 'success';
    const riskLabel = riskState === 'danger' ? '需重点关注' : riskState === 'warning' ? '轻度风险' : '质控平稳';
    summaryBox.innerHTML = `
      <div class="qc-summary-layout">
        <article class="qc-summary-primary qc-summary-primary-${riskState}">
          <div class="qc-summary-primary-head">
            <span class="qc-summary-kicker">质控总览</span>
            <span class="qc-summary-badge qc-summary-badge-${riskState}">${escapeHtml(riskLabel)}</span>
          </div>
          <strong>${escapeHtml(primaryItem.value)}</strong>
          <p>${escapeHtml(primaryItem.label)}</p>
        </article>
        <div class="qc-summary-grid">
          ${secondaryItems.map((item) => `
            <article class="qc-summary-metric">
              <span>${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(item.value)}</strong>
            </article>
          `).join("")}
        </div>
      </div>
    `;
  }
  const fastpTabGroup = document.querySelector(".report-tabs");
  const insertTabButton = fastpTabGroup?.querySelector('[data-report-tab="insert-size"]');
  const baseTabButton = fastpTabGroup?.querySelector('[data-report-tab="base-content"]');
  const adapterTabButton = fastpTabGroup?.querySelector('[data-report-tab="adapter"]');
  if (insertTabButton) {
    insertTabButton.style.display = paired ? "" : "none";
    insertTabButton.classList.toggle("active", paired);
  }
  if (baseTabButton) {
    baseTabButton.classList.toggle("active", !paired);
  }
  if (adapterTabButton && !paired) {
    adapterTabButton.classList.remove("active");
  }
  document.querySelectorAll("#section-fastp .report-tab-panel").forEach((panel) => {
    panel.classList.toggle("active", paired ? panel.dataset.reportPanel === "insert-size" : panel.dataset.reportPanel === "base-content");
  });
  renderFastpTabs(fastp, paired);
  if (paired) {
    initializeBaseTabs();
  }
  bindChartExportButtons();
}

function isPathoSourceReport(data) {
  return String(data?.task?.report_kind || "").trim() === "pathosource_phylogeny"
    || String(data?.task?.workstation_key || "").trim() === "pathosource";
}

function isCommunityReport(data) {
  return String(data?.task?.report_kind || "").trim() === "community_meta_ecology"
    || String(data?.task?.workstation_key || "").trim() === "community";
}

function isSarsCov2NextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "sars_cov_2_nextclade";
}

function isMonkeypoxNextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "monkeypox_nextclade";
}

function isRsvNextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "rsv_nextclade";
}

function isHmpvNextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "hmpv_nextclade";
}

function isDenvNextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "denv_nextclade";
}

function isZikavNextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "zikav_nextclade";
}

function isChikvNextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "chikv_nextclade";
}

function isEbolaNextcladeReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "ebola_nextclade";
}

function isHpivTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "hpiv_typing";
}

function isHadvTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "hadv_typing";
}

function isNorovirusTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "norovirus_typing";
}

function isEnterovirusTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "enterovirus_typing";
}

function isHivTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "hiv_resistance";
}

function resolveSerotypeMode(dataOrSerotype) {
  const serotype = dataOrSerotype?.sections?.serotype && typeof dataOrSerotype.sections.serotype === "object"
    ? dataOrSerotype.sections.serotype
    : (dataOrSerotype && typeof dataOrSerotype === "object" ? dataOrSerotype : {});
  const explicitMode = String(serotype?.mode || "").trim();
  if (explicitMode) return explicitMode;
  const predictedGroup = String(serotype?.predicted_group || "").trim().toUpperCase();
  const predictedClade = String(serotype?.predicted_clade || "").trim().toUpperCase();
  const referenceName = String(serotype?.reference_name || serotype?.selected_reference || "").trim().toLowerCase();
  const summaryCards = Array.isArray(serotype?.summary_cards) ? serotype.summary_cards : [];
  const summaryLabels = summaryCards.map((item) => String(item?.label || "").trim().toLowerCase());
  const summaryValues = summaryCards.map((item) => String(item?.value || "").trim().toLowerCase());
  if (
    summaryLabels.includes("hav 子亚型")
    || summaryLabels.includes("hav子亚型")
    || summaryLabels.includes("子亚型")
    || summaryLabels.includes("大亚型")
    || ["HAV", "HBV", "HCV", "HDV", "HEV"].includes(predictedGroup)
    || ["IA", "IB", "IIA", "IIB", "IIIA", "IIIB"].includes(predictedClade)
    || referenceName.includes("hepatitis a virus")
    || referenceName.includes("hepatitis b virus")
    || referenceName.includes("hepatitis c virus")
    || referenceName.includes("hepatitis d virus")
    || referenceName.includes("hepatitis e virus")
    || referenceName.includes("hepatovirus")
    || summaryValues.some((value) => ["hav", "hbv", "hcv", "hdv", "hev"].includes(value) || /^(hav|hbv|hcv|hdv|hev)\s/.test(value))
  ) {
    return "hepatovirus_typing";
  }
  return "";
}

function isHepatovirusTypingReport(data) {
  return resolveSerotypeMode(data) === "hepatovirus_typing";
}

function isBandavirusTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "bandavirus_typing";
}

function isOrthohantavirusTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "orthohantavirus_typing";
}

function isOrthoebolavirusTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "orthoebolavirus_typing";
}

function isAstroviridaeTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "astroviridae_typing";
}

function isRhinovirusTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "rhinovirus_typing";
}

function isSeasonalHcovTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "seasonal_hcov_typing";
}

function isRotavirusTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "rotavirus_typing";
}

function isInfluenzaTypingReport(data) {
  return String(data?.sections?.serotype?.mode || "").trim() === "influenza_typing";
}

function buildCommunityNav() {
  const groups = [
    {
      section: "section-community-qc",
      title: "1. 数据质控与过滤",
      id: "nav-group-community-qc",
      children: [
        { href: "#community-demux-table", label: "1.1 demux 读数预览" },
        { href: "#community-denoise-table", label: "1.2 去噪结果预览" },
        { href: "#community-qc-rarefaction-chart", label: "1.3 稀释曲线" },
      ],
    },
    {
      section: "section-community-composition",
      title: "2. 物种组成",
      id: "nav-group-community-composition",
      children: [
        { href: "#community-taxonomy-chart", label: "2.1 丰度图" },
        { href: "#community-taxonomy-table", label: "2.2 丰度表" },
      ],
    },
    {
      section: "section-community-alpha",
      title: "3. Alpha 多样性分析",
      id: "nav-group-community-alpha",
      children: [
        { href: "#community-alpha-chart", label: "3.1 Alpha 图形" },
        { href: "#community-alpha-pairwise-table", label: "3.2 两两比较" },
        { href: "#community-alpha-table", label: "3.3 样本预览" },
      ],
    },
    {
      section: "section-community-beta",
      title: "4. Beta 多样性分析",
      id: "nav-group-community-beta",
      children: [
        { href: "#community-beta-pcoa-figure", label: "4.1 PCoA 排序" },
        { href: "#community-beta-nmds-figure", label: "4.2 NMDS 排序" },
        { href: "#community-beta-distance-figure", label: "4.3 距离热图" },
        { href: "#community-beta-stats-table", label: "4.4 显著性统计" },
      ],
    },
    {
      section: "section-community-biomarker",
      title: "5. Biomarker",
      id: "nav-group-community-biomarker",
      children: [
        { href: "#community-biomarker-summary-grid", label: "5.1 结果摘要" },
        { href: "#community-differential-table", label: "5.2 LEfSe 结果" },
        { href: "#community-rf-table", label: "5.3 RF 结果" },
      ],
    },
    {
      section: "section-community-network",
      title: "6. Network",
      id: "nav-group-community-network",
      children: [
        { href: "#community-network-graph", label: "6.1 共现网络" },
        { href: "#community-network-hubs", label: "6.2 Hub taxa" },
        { href: "#community-network-module-table", label: "6.3 模块统计" },
      ],
    },
    {
      section: "section-community-metadata",
      title: "7. 元数据概况",
      id: "nav-group-community-metadata",
      children: [
        { href: "#community-metadata-table", label: "7.1 元数据表" },
      ],
    },
    {
      section: "section-community-notes",
      title: "8. 分析说明",
      id: "nav-group-community-notes",
      children: [
        { href: "#community-modules-table", label: "8.1 模块规划" },
        { href: "#community-outputs-table", label: "8.2 输出产物" },
        { href: "#community-notes-card", label: "8.3 判读说明" },
      ],
    },
  ];
  return groups.map((group) => `
    <div class="report-nav-group has-children" data-nav-group>
      <button class="report-nav-link report-nav-toggle" type="button" data-nav-toggle data-nav-section="${group.section}" aria-expanded="false" aria-controls="${group.id}">
        <span>${group.title}</span>
      </button>
      <div class="report-subnav" id="${group.id}" hidden>
        ${group.children.map((item) => `<a class="report-nav-link subnav-link" href="${item.href}">${item.label}</a>`).join("")}
      </div>
    </div>
  `).join("");
}

function buildCommunityLayout() {
  return `
    <section class="report-intro-card">
      <div class="report-intro-main">
        <p class="report-kicker">Community Ecology Report</p>
        <h2 id="report-sample-title">多样本群落后分析结果</h2>
        <p id="report-sample-copy">围绕 QIIME2 amplicon 实际流程整理测序质控、特征表、taxonomy 注释与差异分析结果，保持和其他样本结果页一致的浏览体验。</p>
        <div id="report-sample-switcher" class="report-sample-switcher hidden" aria-label="样本列表"></div>
      </div>
      <dl class="report-meta-grid">
        <div><dt>任务归属</dt><dd id="meta-owner">-</dd></div>
        <div><dt>用户组</dt><dd id="meta-group">-</dd></div>
        <div><dt>分析模式</dt><dd id="meta-asm-type">-</dd></div>
        <div><dt>核心流程</dt><dd id="meta-method">-</dd></div>
        <div><dt>输入路径</dt><dd id="meta-input">-</dd></div>
        <div><dt>输出目录</dt><dd id="meta-output">-</dd></div>
      </dl>
    </section>

    <section id="section-overview" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Overview</p>
        <h2>数据概览</h2>
        <p>概述样本规模、推荐参数和当前已生成结果，用于快速了解本次群落分析的整体情况。</p>
      </div>
      <div id="overview-metrics" class="metric-grid"></div>
      <div id="community-summary-grid" class="mini-stat-grid"></div>
    </section>

    <section id="section-community-qc" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 1</p>
        <h2>数据质控与过滤</h2>
        <p>展示 demux 统计、DADA2 去噪建议参数和质控相关输出，帮助判断测序数据是否适合继续分析。</p>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>质控摘要</h3>
            <span class="card-tag">Demux / DADA2</span>
          </div>
          <div class="chart-insight" role="note" aria-label="质控摘要说明">
            <span class="chart-insight-label">结果说明</span>
            <p>这里汇总原始测序读长统计和推荐过滤参数，用于快速判断样本测序量是否足够、截断长度是否合理。</p>
          </div>
          <div id="community-demux-grid" class="mini-stat-grid"></div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <h3>质控结果文件</h3>
            <span class="card-tag">输出文件</span>
          </div>
          <div class="chart-insight" role="note" aria-label="质控结果文件说明">
            <span class="chart-insight-label">结果说明</span>
            <p>这里展示可直接打开的原始 QIIME2 质控结果文件，适合继续查看官方可视化细节与导出内容。</p>
          </div>
          <div id="community-qc-assets" class="empty-box">
            <p>显示当前已经生成的 demux 和 DADA2 质控输出。</p>
          </div>
        </article>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>去噪过滤摘要</h3>
            <span class="card-tag">DADA2</span>
          </div>
          <div class="chart-insight" role="note" aria-label="去噪过滤摘要说明">
            <span class="chart-insight-label">结果说明</span>
            <p>这里重点看过滤保留率、非嵌合保留率和 merged 数量，用于判断去噪后可用于下游分析的数据保留情况。</p>
          </div>
          <div id="community-denoise-grid" class="mini-stat-grid"></div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <h3>稀释曲线</h3>
            <span class="card-tag">Rarefaction</span>
          </div>
          <div id="community-qc-rarefaction-chart" class="empty-box">
            <p>显示 observed features 随测序深度变化的稀释曲线。</p>
          </div>
        </article>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>过滤说明</h3>
            <span class="card-tag">判读提示</span>
          </div>
          <div id="community-qc-note" class="empty-box">
            <p>读取 demux.qzv 与 denoising-stats-dada2.qzv 中的结构化结果，直接展示过滤前后统计。</p>
          </div>
        </article>
      </div>
      <div class="chart-insight" role="note" aria-label="demux 表格说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表列出每个样本的正反向 reads 数，用来识别极低测序量样本、0 reads 样本或测序量分布异常的样本。</p>
      </div>
      <div id="community-demux-table" class="report-table-card"></div>
      <div class="chart-insight" role="note" aria-label="去噪结果表格说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表展示每个样本从输入到过滤、合并和去除嵌合体后的保留情况，是判断数据质量和样本可用性的核心明细表。</p>
      </div>
      <div id="community-denoise-table" class="report-table-card"></div>
    </section>

    <section id="section-community-composition" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 2</p>
        <h2>物种组成</h2>
        <p>展示 taxonomy 注释摘要、物种组成预览和 barplot 相关输出，方便直接查看样本群落结构。</p>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>注释摘要</h3>
            <span class="card-tag">Taxonomy</span>
          </div>
          <div class="chart-insight" role="note" aria-label="注释摘要说明">
            <span class="chart-insight-label">结果说明</span>
            <p>这里概括当前分类组成的层级汇总和注释情况，用于判断不同分类水平上哪些类群在样本中更占优势。</p>
          </div>
          <div id="community-taxonomy-grid" class="mini-stat-grid"></div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <h3>物种组成结果文件</h3>
            <span class="card-tag">输出文件</span>
          </div>
          <div class="chart-insight" role="note" aria-label="物种组成结果文件说明">
            <span class="chart-insight-label">结果说明</span>
            <p>这里保留 taxonomy 与 taxa barplot 的原始结果入口，适合继续查看组成结构和官方交互可视化。</p>
          </div>
          <div id="community-composition-assets" class="empty-box">
            <p>显示当前已经生成的 taxonomy 与 taxa barplot 结果文件。</p>
          </div>
        </article>
      </div>
      <div class="chart-insight" role="note" aria-label="taxonomy 表格说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表按不同分类水平汇总平均相对丰度和检出样本数，适合用来快速判断门、科、属、种层面最主要的优势类群。</p>
      </div>
      <div class="report-tabs community-rank-tabs" role="tablist" aria-label="分类水平切换">
        <button class="report-tab-button active" type="button" data-community-rank="门">门</button>
        <button class="report-tab-button" type="button" data-community-rank="科">科</button>
        <button class="report-tab-button" type="button" data-community-rank="属">属</button>
      </div>
      <div class="community-chart-toolbar">
        <label class="community-chart-control">
          <span>显示方式</span>
          <select id="community-view-select">
            <option value="sample" selected>按样本</option>
            <option value="group">按分组合并</option>
          </select>
        </label>
        <label class="community-chart-control">
          <span>Top N</span>
          <select id="community-topn-select">
            <option value="5">5</option>
            <option value="10" selected>10</option>
            <option value="15">15</option>
          </select>
        </label>
        <label class="community-chart-control">
          <span>样本排序</span>
          <select id="community-sort-select">
            <option value="group">按分组</option>
            <option value="sample">按样本名</option>
          </select>
        </label>
      </div>
      <div id="community-taxonomy-chart" class="empty-box">
        <p>显示当前分类水平下的平均相对丰度分布。</p>
      </div>
      <div id="community-taxonomy-table" class="report-table-card"></div>
    </section>

    <section id="section-community-alpha" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 3</p>
        <h2>Alpha 多样性分析</h2>
        <p>展示稀释深度建议与 Alpha 指标变化，便于结合分组信息查看群落多样性变化。</p>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>Alpha 摘要</h3>
            <span class="card-tag">Rarefaction</span>
          </div>
          <div class="chart-insight" role="note" aria-label="Alpha 摘要说明">
            <span class="chart-insight-label">结果说明</span>
            <p>这里先给出稀释深度和显著条目等概览信息，帮助快速判断 alpha 多样性和差异分析的观察起点。</p>
          </div>
          <div id="community-alpha-grid" class="mini-stat-grid"></div>
        </article>
      </div>
      <div class="chart-insight" role="note" aria-label="Alpha 图形说明">
        <span class="chart-insight-label">结果说明</span>
        <p>这里按指标切换查看 Alpha 箱线图。先看组内分布，再看整体显著性，最后用任意两组的 P 值表确认具体是哪些分组之间存在差异。</p>
      </div>
      <div class="report-tabs community-alpha-metric-tabs" role="tablist" aria-label="Alpha 指标切换"></div>
      <div class="community-alpha-workspace">
        <div id="community-alpha-chart" class="report-card-slot"></div>
        <aside id="community-alpha-detail" class="result-card community-alpha-detail-card"></aside>
      </div>
      <div class="chart-insight" role="note" aria-label="Alpha 两两比较说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表列出当前指标下任意两组之间的 P 值，方便快速定位具体是哪几组之间存在显著性差异。</p>
      </div>
      <div id="community-alpha-pairwise-table" class="report-table-card"></div>
      <div class="chart-insight" role="note" aria-label="Alpha 样本预览说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表给出样本级 Alpha 指标预览，适合核对异常样本以及判断组内波动是否过大。</p>
      </div>
      <div id="community-alpha-table" class="report-table-card"></div>
    </section>

    <section id="section-community-beta" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 4</p>
        <h2>Beta 多样性分析</h2>
        <p>接入 microeco 的 beta 多样性结果，在报告页直接展示 PCoA 排序、组内距离分布和组间显著性统计。</p>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>Beta 摘要</h3>
            <span class="card-tag">microeco</span>
          </div>
          <div class="chart-insight" role="note" aria-label="Beta 摘要说明">
            <span class="chart-insight-label">结果说明</span>
            <p>先看 Bray-Curtis 距离下的 PERMANOVA、ANOSIM 和 betadisper，判断组间是否显著分离以及组内离散度是否均衡。</p>
          </div>
          <div id="community-beta-grid" class="mini-stat-grid"></div>
        </article>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>PCoA 排序图</h3>
            <span class="card-tag">Bray-Curtis</span>
          </div>
          <div id="community-beta-pcoa-figure" class="empty-box">
            <p>显示 microeco 输出的 PCoA 排序图。</p>
          </div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <h3>NMDS 排序图</h3>
            <span class="card-tag">Bray-Curtis</span>
          </div>
          <div id="community-beta-nmds-figure" class="empty-box">
            <p>显示 microeco 输出的 NMDS 排序图。</p>
          </div>
        </article>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>距离矩阵热图</h3>
            <span class="card-tag">Clustered Heatmap</span>
          </div>
          <div id="community-beta-distance-figure" class="empty-box">
            <p>显示 Bray-Curtis 距离矩阵热图，并对横纵轴样本同时进行层次聚类排序。</p>
          </div>
        </article>
      </div>
      <article class="result-card">
        <div class="card-head">
          <h3>聚类树与样本物种组成</h3>
          <span class="card-tag">Cluster + Composition</span>
        </div>
        <div id="community-beta-composition-figure" class="empty-box">
          <p>按 Bray-Curtis 聚类顺序展示样本组成堆积条形图。</p>
        </div>
      </article>
      <div class="chart-insight" role="note" aria-label="Beta 统计说明">
        <span class="chart-insight-label">统计说明</span>
        <p>PERMANOVA 反映整体组间差异，ANOSIM 反映排序分离程度，betadisper 用来判断组内离散度是否均衡；三者结合更适合解释 beta 多样性结果。</p>
      </div>
      <div id="community-beta-stats-table" class="report-table-card"></div>
      <div class="chart-insight" role="note" aria-label="组内距离统计说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表总结组内距离差异检验结果和各组纳入样本数，可与聚类热图一起判断分组内部稳定性和样本间相似性结构。</p>
      </div>
      <div id="community-beta-distance-table" class="report-table-card"></div>
    </section>

    <section id="section-community-biomarker" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 5</p>
        <h2>Biomarker</h2>
        <p>单独汇总 LEfSe 与 RF Biomarker 的关键结果，帮助快速定位最具分组判别力的微生物特征。</p>
      </div>
      <div id="community-biomarker-summary-grid" class="mini-stat-grid"></div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>LEfSe LDA 条形图</h3>
            <span class="card-tag">Interactive</span>
          </div>
          <div id="community-biomarker-lefse-chart" class="empty-box">
            <p>显示 LEfSe 的 LDA 排序结果。</p>
          </div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <h3>RF 重要性排序</h3>
            <span class="card-tag">Interactive</span>
          </div>
          <div id="community-biomarker-rf-chart" class="empty-box">
            <p>显示 Random Forest 的特征重要性排序。</p>
          </div>
        </article>
      </div>
      <div class="chart-insight" role="note" aria-label="Biomarker 表格说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表汇总 LEfSe 与 RF 识别出的重点特征，适合快速查看差异分类单元、分组方向以及特征重要性。</p>
      </div>
      <div id="community-differential-table" class="report-table-card"></div>
      <div id="community-rf-table" class="report-table-card"></div>
    </section>

    <section id="section-community-network" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 6</p>
        <h2>Network</h2>
        <p>基于 microeco trans_network 构建共现网络，直接查看模块结构、关键 hub taxa 与节点连接特征。</p>
      </div>
      <div id="community-network-summary-grid" class="mini-stat-grid"></div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>共现网络图</h3>
            <span class="card-tag">Interactive</span>
          </div>
          <div id="community-network-graph" class="empty-box">
            <p>显示按模块分组的共现网络结构。</p>
          </div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <h3>Hub taxa 排序</h3>
            <span class="card-tag">Degree</span>
          </div>
          <div id="community-network-hubs" class="empty-box">
            <p>显示度中心性最高的关键节点。</p>
          </div>
        </article>
      </div>
      <div class="two-column taxonomy-summary-grid">
        <article class="result-card">
          <div class="card-head">
            <h3>节点角色散点图</h3>
            <span class="card-tag">Zi-Pi</span>
          </div>
          <div id="community-network-role-scatter" class="empty-box">
            <p>按 within-module 和 among-module connectivity 展示节点角色分布。</p>
          </div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <h3>Genus 分层角色图</h3>
            <span class="card-tag">Genus</span>
          </div>
          <div id="community-network-role-layered" class="empty-box">
            <p>按 Genus 查看 among-module 与 within-module connectivity 的分层分布。</p>
          </div>
        </article>
      </div>
      <div class="chart-insight" role="note" aria-label="Network 模块统计说明">
        <span class="chart-insight-label">表格说明</span>
        <p>模块表优先看节点数和平均 Degree，判断哪个模块更大、更紧密；节点表和边表适合进一步核对关键属和具体关联关系。</p>
      </div>
      <div id="community-network-module-table" class="report-table-card"></div>
      <div id="community-network-node-table" class="report-table-card"></div>
      <div id="community-network-edge-table" class="report-table-card"></div>
    </section>

    <section id="section-community-metadata" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 7</p>
        <h2>元数据概况</h2>
        <p>确认样本 ID、分组列和元数据字段是否完整，这决定了后续 alpha / beta 和差异比较能否稳定开展。</p>
      </div>
      <div class="chart-insight" role="note" aria-label="元数据表格说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表用于核对样本分组、元数据字段和分析参数是否一致，是解释后续群落差异结果的基础上下文。</p>
      </div>
      <div id="community-metadata-table" class="report-table-card"></div>
    </section>

    <section id="section-community-notes" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 8</p>
        <h2>分析说明</h2>
        <p>补充当前流程说明、结果阅读建议和后续可继续展开的分析方向。</p>
      </div>
      <div class="chart-insight" role="note" aria-label="模块规划说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表记录当前群落分析模块和状态，用来说明本次任务已经跑到哪里，以及还可以继续延伸哪些分析。</p>
      </div>
      <div id="community-modules-table" class="report-table-card"></div>
      <div class="chart-insight" role="note" aria-label="输出产物表格说明">
        <span class="chart-insight-label">表格说明</span>
        <p>这张表列出当前结果目录中已经生成的重要产物，便于快速核对可交付文件和后续复核入口。</p>
      </div>
      <div id="community-outputs-table" class="report-table-card"></div>
      <div id="community-notes-card" class="result-card"></div>
    </section>
  `;
}

function setCommunityReportChrome(task) {
  const shell = document.querySelector(".report-shell");
  if (shell) {
    shell.dataset.reportKind = "community";
  }
  document.querySelectorAll(".report-scene-switcher").forEach((node) => node.classList.add("hidden"));
  const clinicalNav = document.getElementById("nav-group-clinical-wrapper");
  if (clinicalNav) clinicalNav.classList.add("hidden");
  const title = document.querySelector(".report-title-block h1");
  if (title) {
    title.textContent = "多样本群落后分析结果";
  }
  const subtitle = document.querySelector(".report-subtitle");
  if (subtitle) {
    subtitle.textContent = "本报告按 QIIME2 amplicon 主流程组织 demux、DADA2、taxonomy、alpha / beta 多样性与 microeco Biomarker（LEfSe / RF）结果，并沿用统一结果页样式展示。";
  }
  if (task) {
    task.asm_type = "群落后分析";
    task.method = "QIIME2 amplicon";
  }
}

function renderCommunityNotes(containerId, notes) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const items = Array.isArray(notes) ? notes.filter(Boolean) : [];
  if (!items.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有额外说明</strong>
        <p class="empty-copy">后续如果接通真实 QIIME2 / microeco 执行结果，这里会补充运行说明与判读提示。</p>
      </div>
    `;
    return;
  }
  container.innerHTML = `
    <div class="chart-insight" role="note" aria-label="群落分析说明">
      <span class="chart-insight-label">分析说明</span>
      <ul class="report-bullet-list">
        ${items.map((note) => `<li>${escapeHtml(String(note))}</li>`).join("")}
      </ul>
    </div>
  `;
}

function normalizeCommunityAssetRows(items) {
  const outputRoot = String(document.getElementById("meta-output")?.textContent || "").trim();
  return Array.isArray(items)
    ? items
      .map((item) => {
        const fullPath = String(item?.path || "").trim();
        let relativePath = fullPath;
        if (outputRoot && fullPath.startsWith(`${outputRoot}/`)) {
          relativePath = fullPath.slice(outputRoot.length + 1);
        } else if (outputRoot && fullPath === outputRoot) {
          relativePath = "";
        }
        return {
          ...item,
          relativePath,
        };
      })
      .filter((item) => item && item.label && item.relativePath)
    : [];
}

function renderCommunityAssetLinks(containerId, taskId, items, emptyTitle, emptyCopy) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const rows = normalizeCommunityAssetRows(items);
  if (!rows.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>${escapeHtml(emptyTitle)}</strong>
        <p class="empty-copy">${escapeHtml(emptyCopy)}</p>
      </div>
    `;
    return;
  }
  container.innerHTML = `
    <div class="community-asset-grid">
      ${rows.map((item) => `
        <a class="community-asset-link" href="/api/tasks/${encodeURIComponent(taskId)}/report-asset/${encodeURIComponent(item.relativePath)}" target="_blank" rel="noopener noreferrer">
          <span class="community-asset-name">${escapeHtml(item.label)}</span>
          <span class="community-asset-copy">点击打开结果文件，进入原始可视化或导出内容。</span>
          <strong class="community-asset-status">${escapeHtml(item.status === "ready" ? "已生成" : item.statusLabel)}</strong>
        </a>
      `).join("")}
    </div>
  `;
}

function renderCommunityImageAsset(containerId, taskId, item, title, copy) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const normalized = normalizeCommunityAssetRows(item ? [item] : []);
  const target = normalized[0];
  if (!target) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>${escapeHtml(title || "当前暂无可展示图形")}</strong>
        <p class="empty-copy">${escapeHtml(copy || "请先生成对应图形结果后再查看。")}</p>
      </div>
    `;
    return;
  }
  const assetUrl = `/api/tasks/${encodeURIComponent(taskId)}/report-asset/${encodeURIComponent(target.relativePath)}`;
  container.classList.remove("empty-box");
  container.innerHTML = `
    <figure class="community-image-figure">
      <img class="community-report-image" src="${assetUrl}" alt="${escapeHtml(target.label || title || "图形结果")}">
      <figcaption class="community-image-caption">
        <strong>${escapeHtml(title || target.label || "图形结果")}</strong>
        <span>${escapeHtml(copy || "点击结果文件入口可继续下载原图或对应 PDF。")}</span>
      </figcaption>
    </figure>
  `;
}

function renderCommunityBetaOrdination(containerId, points, options = {}) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const plotTitle = String(options?.title || "排序图");
  const emptyLabel = String(options?.emptyLabel || plotTitle);
  const xLabel = String(options?.xLabel || "Axis 1");
  const yLabel = String(options?.yLabel || "Axis 2");
  const insight = String(options?.insight || `当前展示 ${points.length} 个样本在前两排序轴上的位置；悬停点位可查看样本名和分组。`);
  const stress = options?.stress;
  if (!points.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有可展示的${escapeHtml(emptyLabel)}点位</strong>
        <p class="empty-copy">请先生成对应坐标表后再查看交互散点图。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 760;
  const height = 760;
  const padLeft = 92;
  const padRight = 40;
  const padTop = 34;
  const padBottom = 88;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const xValues = points.map((point) => Number(point.x) || 0);
  const yValues = points.map((point) => Number(point.y) || 0);
  const groups = Array.from(new Set(points.map((point) => String(point.group || "未分组")))).sort((left, right) => left.localeCompare(right, "zh-CN"));
  const palette = ["#3e546f", "#8a6654", "#7a7158", "#4f7f6b", "#7c5d83", "#9a6b62", "#557a95", "#a88262", "#5b6b4d"];
  const colorMap = new Map(groups.map((group, index) => [group, palette[index % palette.length]]));
  const chiSquare90 = 4.605170186;
  const ellipseModels = groups.map((group) => {
    const groupPoints = points
      .filter((point) => String(point.group || "未分组") === group)
      .map((point) => ({
        x: Number(point.x) || 0,
        y: Number(point.y) || 0,
      }));
    if (groupPoints.length <= 3) return null;
    const meanX = groupPoints.reduce((sum, point) => sum + point.x, 0) / groupPoints.length;
    const meanY = groupPoints.reduce((sum, point) => sum + point.y, 0) / groupPoints.length;
    let covXX = 0;
    let covYY = 0;
    let covXY = 0;
    groupPoints.forEach((point) => {
      const dx = point.x - meanX;
      const dy = point.y - meanY;
      covXX += dx * dx;
      covYY += dy * dy;
      covXY += dx * dy;
    });
    const divisor = Math.max(groupPoints.length - 1, 1);
    covXX /= divisor;
    covYY /= divisor;
    covXY /= divisor;
    const trace = covXX + covYY;
    const determinant = covXX * covYY - covXY * covXY;
    const discriminant = Math.max((trace * trace) / 4 - determinant, 0);
    const root = Math.sqrt(discriminant);
    const lambda1 = Math.max(trace / 2 + root, 0);
    const lambda2 = Math.max(trace / 2 - root, 0);
    if ((lambda1 <= 0 && lambda2 <= 0) || !Number.isFinite(lambda1) || !Number.isFinite(lambda2)) return null;
    let vectorX = covXY;
    let vectorY = lambda1 - covXX;
    if (Math.abs(vectorX) < 1e-9 && Math.abs(vectorY) < 1e-9) {
      vectorX = 1;
      vectorY = 0;
    }
    const vectorNorm = Math.hypot(vectorX, vectorY) || 1;
    const axis1X = vectorX / vectorNorm;
    const axis1Y = vectorY / vectorNorm;
    const axis2X = -axis1Y;
    const axis2Y = axis1X;
    const radius1 = Math.sqrt(lambda1 * chiSquare90);
    const radius2 = Math.sqrt(lambda2 * chiSquare90);
    return {
      group,
      meanX,
      meanY,
      axis1X,
      axis1Y,
      axis2X,
      axis2Y,
      radius1,
      radius2,
      count: groupPoints.length,
    };
  }).filter(Boolean);
  const ellipseBounds = ellipseModels.reduce((bounds, model) => {
    const xExtent = Math.sqrt((model.radius1 * model.axis1X) ** 2 + (model.radius2 * model.axis2X) ** 2);
    const yExtent = Math.sqrt((model.radius1 * model.axis1Y) ** 2 + (model.radius2 * model.axis2Y) ** 2);
    bounds.minX = Math.min(bounds.minX, model.meanX - xExtent);
    bounds.maxX = Math.max(bounds.maxX, model.meanX + xExtent);
    bounds.minY = Math.min(bounds.minY, model.meanY - yExtent);
    bounds.maxY = Math.max(bounds.maxY, model.meanY + yExtent);
    return bounds;
  }, {
    minX: Math.min(...xValues),
    maxX: Math.max(...xValues),
    minY: Math.min(...yValues),
    maxY: Math.max(...yValues),
  });
  const minX = ellipseBounds.minX;
  const maxX = ellipseBounds.maxX;
  const minY = ellipseBounds.minY;
  const maxY = ellipseBounds.maxY;
  const rangeX = Math.max(maxX - minX, 1e-6);
  const rangeY = Math.max(maxY - minY, 1e-6);
  const xPadding = rangeX * 0.1;
  const yPadding = rangeY * 0.1;
  const domainXMin = minX - xPadding;
  const domainXMax = maxX + xPadding;
  const domainYMin = minY - yPadding;
  const domainYMax = maxY + yPadding;
  const xScale = (value) => padLeft + ((value - domainXMin) / Math.max(domainXMax - domainXMin, 1e-6)) * innerWidth;
  const yScale = (value) => padTop + innerHeight - ((value - domainYMin) / Math.max(domainYMax - domainYMin, 1e-6)) * innerHeight;
  const xTicks = 4;
  const yTicks = 4;
  const zeroX = domainXMin <= 0 && domainXMax >= 0 ? xScale(0) : null;
  const zeroY = domainYMin <= 0 && domainYMax >= 0 ? yScale(0) : null;
  const gridX = Array.from({ length: xTicks + 1 }, (_, index) => {
    const ratio = index / xTicks;
    const value = domainXMin + (domainXMax - domainXMin) * ratio;
    const x = padLeft + innerWidth * ratio;
    return `
      <line class="chart-grid-line" x1="${x}" y1="${padTop}" x2="${x}" y2="${height - padBottom}"></line>
      <text class="chart-axis-label x-axis" x="${x}" y="${height - 28}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const gridY = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const value = domainYMin + (domainYMax - domainYMin) * ratio;
    const y = padTop + innerHeight - innerHeight * ratio;
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const ellipses = ellipseModels.map((model) => {
    const samples = 64;
    const path = Array.from({ length: samples }, (_, index) => {
      const theta = (Math.PI * 2 * index) / samples;
      const dx = model.radius1 * Math.cos(theta) * model.axis1X + model.radius2 * Math.sin(theta) * model.axis2X;
      const dy = model.radius1 * Math.cos(theta) * model.axis1Y + model.radius2 * Math.sin(theta) * model.axis2Y;
      const sx = xScale(model.meanX + dx);
      const sy = yScale(model.meanY + dy);
      return `${index === 0 ? "M" : "L"} ${sx} ${sy}`;
    }).join(" ");
    return `
      <path
        class="community-pcoa-ellipse"
        d="${path} Z"
        fill="${escapeHtml(colorMap.get(model.group) || "#3e546f")}"
        stroke="${escapeHtml(colorMap.get(model.group) || "#3e546f")}"
        data-group="${escapeHtml(model.group)}"
        data-count="${model.count}"
      ></path>
    `;
  }).join("");
  const dots = points.map((point) => {
    const x = xScale(Number(point.x) || 0);
    const y = yScale(Number(point.y) || 0);
    const group = String(point.group || "未分组");
    return `
      <circle
        class="community-pcoa-point"
        cx="${x}"
        cy="${y}"
        r="5.5"
        fill="${escapeHtml(colorMap.get(group) || "#3e546f")}"
        fill-opacity="0.86"
        data-sample="${escapeHtml(String(point.sample || "--"))}"
        data-group="${escapeHtml(group)}"
        data-x="${escapeHtml(String(point.x ?? ""))}"
        data-y="${escapeHtml(String(point.y ?? ""))}"
      >
        <title>${escapeHtml(`${point.sample} | ${group} | PCo1=${point.x} | PCo2=${point.y}`)}</title>
      </circle>
    `;
  }).join("");
  container.innerHTML = `
    <div class="mini-chart-card relation-chart-card community-pcoa-card">
      <span class="mini-chart-title">${escapeHtml(plotTitle)}</span>
      ${buildChartInsight(`${escapeHtml(insight)}${stress != null ? ` 当前 NMDS stress = ${escapeHtml(String(stress))}。` : ""}`)}
      <div class="chart-canvas community-pcoa-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg community-pcoa-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(plotTitle)}">
          ${gridX}
          ${gridY}
          ${zeroX != null ? `<line class="chart-axis-line community-pcoa-zero" x1="${zeroX}" y1="${padTop}" x2="${zeroX}" y2="${height - padBottom}"></line>` : ""}
          ${zeroY != null ? `<line class="chart-axis-line community-pcoa-zero" x1="${padLeft}" y1="${zeroY}" x2="${width - padRight}" y2="${zeroY}"></line>` : ""}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
          ${ellipses}
          ${dots}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 6}">${escapeHtml(xLabel)}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(yLabel)}</text>
        </svg>
        <div class="chart-tooltip" hidden></div>
      </div>
      <div class="community-pcoa-legend">
        ${groups.map((group) => `
          <span class="community-pcoa-legend-item">
            <i style="background:${escapeHtml(colorMap.get(group) || "#3e546f")}"></i>
            <span>${escapeHtml(group)}</span>
          </span>
        `).join("")}
      </div>
    </div>
  `;
  const svg = container.querySelector(".community-pcoa-svg");
  const tooltip = container.querySelector(".chart-tooltip");
  const circles = Array.from(container.querySelectorAll(".community-pcoa-point"));
  if (!svg || !tooltip || !circles.length) return;
  const showTooltip = (event, node) => {
    const sample = node.getAttribute("data-sample") || "--";
    const group = node.getAttribute("data-group") || "未分组";
    const x = node.getAttribute("data-x") || "--";
    const y = node.getAttribute("data-y") || "--";
    tooltip.innerHTML = `<strong>${escapeHtml(sample)}</strong><div>分组: ${escapeHtml(group)}</div><div>PCo1: ${escapeHtml(x)}</div><div>PCo2: ${escapeHtml(y)}</div>`;
    tooltip.hidden = false;
    const rect = container.getBoundingClientRect();
    const offsetX = event.clientX - rect.left;
    const offsetY = event.clientY - rect.top;
    tooltip.style.left = `${Math.min(offsetX + 16, rect.width - 180)}px`;
    tooltip.style.top = `${Math.max(offsetY - 18, 12)}px`;
  };
  circles.forEach((node) => {
    node.addEventListener("mouseenter", (event) => {
      circles.forEach((item) => item.classList.toggle("is-dimmed", item !== node));
      node.classList.add("is-active");
      showTooltip(event, node);
    });
    node.addEventListener("mousemove", (event) => showTooltip(event, node));
    node.addEventListener("mouseleave", () => {
      circles.forEach((item) => {
        item.classList.remove("is-dimmed", "is-active");
      });
      tooltip.hidden = true;
    });
  });
}

function parseCommunityBiomarkerRows(rows = [], kind = "lefse") {
  if (!Array.isArray(rows)) return [];
  return rows
    .map((row) => {
      const cells = Array.isArray(row) ? row : [];
      if (kind === "rf") {
        const label = String(cells[1] || "").trim();
        const model = String(cells[2] || "").trim() || "random_forest";
        const value = Number(cells[4]);
        return Number.isFinite(value) && label
          ? { label, group: model, value, pValue: null, method: "RF" }
          : null;
      }
      const label = String(cells[1] || "").trim();
      const group = String(cells[2] || "").trim() || "未分组";
      const pValue = Number(cells[3]);
      const value = Number(cells[4]);
      return Number.isFinite(value) && label
        ? { label, group, value, pValue: Number.isFinite(pValue) ? pValue : null, method: "LEfSe" }
        : null;
    })
    .filter(Boolean)
    .sort((a, b) => Number(b.value || 0) - Number(a.value || 0));
}

function renderCommunityBiomarkerBars(containerId, rows, options = {}) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const chartId = `${containerId}-topn`;
  const kind = String(options.kind || "lefse");
  const label = String(options.label || "Biomarker 排序");
  const xLabel = String(options.xLabel || "Score");
  const parsed = parseCommunityBiomarkerRows(rows, kind);
  if (!parsed.length) {
    container.classList.add("empty-box");
    container.innerHTML = `<p>当前没有可展示的${escapeHtml(label)}数据。</p>`;
    return;
  }
  container.classList.remove("empty-box");
  const topnOptions = [10, 15, 20, 30].filter((value, index, list) => value <= parsed.length || index === 0 || value === list[list.length - 1]);
  const groupPalette = ["#335c67", "#9e2a2b", "#5f0f40", "#588157", "#bc6c25", "#4361ee", "#7f5539", "#6c757d"];
  const groups = Array.from(new Set(parsed.map((item) => String(item.group || "未分组"))));
  const colorMap = new Map(groups.map((group, index) => [group, groupPalette[index % groupPalette.length]]));
  container.innerHTML = `
    <div class="chart-insight" role="note" aria-label="${escapeHtml(label)}说明">
      <span class="chart-insight-label">结果说明</span>
      <p>${escapeHtml(kind === "rf" ? "按重要性从高到低查看最能区分分组的菌群特征。悬浮可查看具体分值，点击条形可高亮关注项。" : "按 LDA 从高到低查看 LEfSe 筛出的差异菌群。颜色表示富集分组，悬浮可查看 P 值与分组。")}</p>
    </div>
    <div class="community-biomarker-controls">
      <label class="community-chart-control">
        <span>显示条目</span>
        <select id="${chartId}">
          ${topnOptions.map((value) => `<option value="${value}" ${value === topnOptions[0] ? "selected" : ""}>Top ${value}</option>`).join("")}
        </select>
      </label>
    </div>
    <div class="mini-chart-card relation-chart-card community-biomarker-card">
      <span class="mini-chart-title">${escapeHtml(label)}</span>
      <div class="chart-canvas community-biomarker-canvas" style="--chart-height:440px">
        <svg class="sparkline-svg community-biomarker-svg" viewBox="0 0 1120 440" preserveAspectRatio="xMidYMid meet" role="img" aria-label="${escapeHtml(label)}"></svg>
        <div class="chart-tooltip" hidden></div>
      </div>
      ${kind === "lefse" ? `
        <div class="community-pcoa-legend">
          ${groups.map((group) => `
            <span class="community-pcoa-legend-item">
              <i style="background:${escapeHtml(colorMap.get(group) || "#3e546f")}"></i>
              <span>${escapeHtml(group)}</span>
            </span>
          `).join("")}
        </div>
      ` : ""}
    </div>
  `;
  const select = document.getElementById(chartId);
  const svg = container.querySelector(".community-biomarker-svg");
  const tooltip = container.querySelector(".chart-tooltip");
  if (!(select instanceof HTMLSelectElement) || !(svg instanceof SVGElement) || !(tooltip instanceof HTMLElement)) return;

  const renderChart = () => {
    const limit = Math.max(1, Number(select.value) || topnOptions[0] || 10);
    const items = parsed.slice(0, limit).reverse();
    const width = 1120;
    const height = 440;
    const padLeft = 300;
    const padRight = 34;
    const padTop = 28;
    const padBottom = 28;
    const innerWidth = width - padLeft - padRight;
    const rowHeight = Math.max(26, Math.floor((height - padTop - padBottom) / Math.max(items.length, 1)));
    const barHeight = Math.min(20, Math.max(12, rowHeight - 8));
    const maxValue = Math.max(...items.map((item) => Number(item.value) || 0), 1);
    const grid = Array.from({ length: 5 }, (_, index) => {
      const ratio = index / 4;
      const value = maxValue * ratio;
      const x = padLeft + innerWidth * ratio;
      return `
        <line class="chart-grid-line" x1="${x}" y1="${padTop}" x2="${x}" y2="${height - padBottom}"></line>
        <text class="chart-axis-label x-axis" x="${x}" y="${height - 6}">${escapeHtml(formatChartValue(value))}</text>
      `;
    }).join("");
    const bars = items.map((item, index) => {
      const value = Number(item.value) || 0;
      const y = padTop + index * rowHeight + (rowHeight - barHeight) / 2;
      const barWidth = (value / maxValue) * innerWidth;
      const color = kind === "rf" ? "#355070" : (colorMap.get(String(item.group || "未分组")) || "#355070");
      const payload = encodeURIComponent(JSON.stringify(item));
      return `
        <g class="community-biomarker-bar-group" data-community-biomarker="${payload}" tabindex="0">
          <text class="chart-axis-label y-axis" x="${padLeft - 12}" y="${y + barHeight / 2 + 4}" text-anchor="end">${escapeHtml(item.label)}</text>
          <rect class="community-biomarker-bar" x="${padLeft}" y="${y}" width="${barWidth}" height="${barHeight}" rx="8" ry="8" fill="${escapeHtml(color)}"></rect>
          <text class="chart-bar-value" x="${padLeft + barWidth + 8}" y="${y + barHeight / 2 + 4}">${escapeHtml(value.toFixed(3))}</text>
        </g>
      `;
    }).join("");
    svg.innerHTML = `
      ${grid}
      <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
      <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
      ${bars}
      <text class="chart-axis-title" x="${(padLeft + width - padRight) / 2}" y="${height - 18}">${escapeHtml(xLabel)}</text>
    `;
    const nodes = Array.from(container.querySelectorAll(".community-biomarker-bar-group"));
    const hideTooltip = () => {
      tooltip.hidden = true;
      nodes.forEach((node) => node.classList.remove("is-active", "is-dimmed"));
    };
    const showTooltip = (event, node) => {
      let data = null;
      try {
        data = JSON.parse(decodeURIComponent(node.getAttribute("data-community-biomarker") || ""));
      } catch (_error) {
        data = null;
      }
      if (!data) return;
      nodes.forEach((item) => item.classList.toggle("is-dimmed", item !== node));
      node.classList.add("is-active");
      tooltip.innerHTML = kind === "rf"
        ? `<strong>${escapeHtml(data.label || "--")}</strong><div>模型: ${escapeHtml(data.group || "RF")}</div><div>重要性: ${escapeHtml(String(data.value ?? "--"))}</div>`
        : `<strong>${escapeHtml(data.label || "--")}</strong><div>分组: ${escapeHtml(data.group || "--")}</div><div>LDA: ${escapeHtml(String(data.value ?? "--"))}</div><div>P值: ${escapeHtml(data.pValue != null ? String(data.pValue) : "--")}</div>`;
      tooltip.hidden = false;
      const rect = container.getBoundingClientRect();
      const offsetX = event.clientX - rect.left;
      const offsetY = event.clientY - rect.top;
      tooltip.style.left = `${Math.min(offsetX + 16, rect.width - 220)}px`;
      tooltip.style.top = `${Math.max(offsetY - 18, 12)}px`;
    };
    nodes.forEach((node) => {
      node.addEventListener("mouseenter", (event) => showTooltip(event, node));
      node.addEventListener("mousemove", (event) => showTooltip(event, node));
      node.addEventListener("mouseleave", hideTooltip);
      node.addEventListener("focus", () => {
        const rect = node.getBoundingClientRect();
        showTooltip({ clientX: rect.left, clientY: rect.top }, node);
      });
      node.addEventListener("blur", hideTooltip);
      node.addEventListener("click", (event) => showTooltip(event, node));
    });
  };
  select.addEventListener("change", renderChart);
  renderChart();
}

function buildCommunityNetworkPalette(labels = []) {
  const palette = ["#355070", "#6d597a", "#b56576", "#e56b6f", "#eaac8b", "#588157", "#4361ee", "#bc6c25", "#7f5539", "#5f0f40"];
  return new Map(labels.map((label, index) => [label, palette[index % palette.length]]));
}

function renderCommunityNetworkGraph(containerId, networkSection = {}) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const nodes = Array.isArray(networkSection?.nodes) ? networkSection.nodes.filter((item) => item?.id) : [];
  const edges = Array.isArray(networkSection?.edges) ? networkSection.edges.filter((item) => item?.source && item?.target) : [];
  if (!nodes.length || !edges.length) {
    container.classList.add("empty-box");
    container.innerHTML = "<p>当前没有可展示的网络节点或边。</p>";
    return;
  }
  container.classList.remove("empty-box");
  const modules = Array.from(new Set(nodes.map((item) => String(item.module || "未分模块"))));
  const colorMap = buildCommunityNetworkPalette(modules);
  const selectId = `${containerId}-limit`;
  const moduleId = `${containerId}-module`;
  container.innerHTML = `
    <div class="chart-insight" role="note" aria-label="共现网络图说明">
      <span class="chart-insight-label">结果说明</span>
      <p>模块颜色表示共现群落分区，节点大小代表 Degree。悬浮可查看菌群、模块、门分类和连接强度，切换 TopN 可以收紧网络视野。</p>
    </div>
    <div class="community-biomarker-controls">
      <label class="community-chart-control">
        <span>显示节点</span>
        <select id="${selectId}">
          <option value="30">Top 30</option>
          <option value="60" selected>Top 60</option>
          <option value="90">Top 90</option>
          <option value="120">Top 120</option>
        </select>
      </label>
      <label class="community-chart-control">
        <span>模块筛选</span>
        <select id="${moduleId}">
          <option value="all">全部模块</option>
          ${modules.map((module) => `<option value="${escapeHtml(module)}">${escapeHtml(module)}</option>`).join("")}
        </select>
      </label>
    </div>
    <div class="mini-chart-card relation-chart-card community-biomarker-card">
      <span class="mini-chart-title">模块化共现网络</span>
      <div class="chart-canvas community-network-canvas" style="--chart-height:560px">
        <svg class="sparkline-svg community-network-svg" viewBox="0 0 1120 560" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Community network"></svg>
        <div class="chart-tooltip" hidden></div>
      </div>
      <div class="community-pcoa-legend">
        ${modules.map((module) => `
          <span class="community-pcoa-legend-item">
            <i style="background:${escapeHtml(colorMap.get(module) || "#355070")}"></i>
            <span>${escapeHtml(module)}</span>
          </span>
        `).join("")}
      </div>
    </div>
  `;
  const select = document.getElementById(selectId);
  const moduleSelect = document.getElementById(moduleId);
  const svg = container.querySelector(".community-network-svg");
  const tooltip = container.querySelector(".chart-tooltip");
  if (!(select instanceof HTMLSelectElement) || !(moduleSelect instanceof HTMLSelectElement) || !(svg instanceof SVGElement) || !(tooltip instanceof HTMLElement)) return;

  const render = () => {
    const limit = Math.max(10, Number(select.value) || 60);
    const moduleFilter = String(moduleSelect.value || "all");
    const rankedNodes = nodes
      .filter((item) => moduleFilter === "all" || String(item.module || "未分模块") === moduleFilter)
      .sort((left, right) => (Number(right.degree) || 0) - (Number(left.degree) || 0))
      .slice(0, limit);
    const nodeIds = new Set(rankedNodes.map((item) => String(item.id)));
    const rankedEdges = edges
      .filter((item) => nodeIds.has(String(item.source)) && nodeIds.has(String(item.target)))
      .sort((left, right) => Math.abs(Number(right.weight) || 0) - Math.abs(Number(left.weight) || 0))
      .slice(0, Math.max(120, limit * 3));
    if (!rankedNodes.length || !rankedEdges.length) {
      svg.innerHTML = `<text x="560" y="280" text-anchor="middle" class="chart-axis-label">当前筛选条件下没有可绘制的网络结构</text>`;
      return;
    }
    const width = 1120;
    const height = 560;
    const centerX = width / 2;
    const centerY = height / 2;
    const grouped = Array.from(
      rankedNodes.reduce((map, item) => {
        const key = String(item.module || "未分模块");
        const list = map.get(key) || [];
        list.push(item);
        map.set(key, list);
        return map;
      }, new Map()).entries()
    ).sort((left, right) => right[1].length - left[1].length);
    const moduleCenters = new Map();
    grouped.forEach(([module], index) => {
      const angle = (Math.PI * 2 * index) / Math.max(grouped.length, 1) - Math.PI / 2;
      const radius = Math.min(180, 90 + grouped.length * 18);
      moduleCenters.set(module, {
        x: centerX + Math.cos(angle) * radius,
        y: centerY + Math.sin(angle) * radius,
      });
    });
    const positionMap = new Map();
    grouped.forEach(([module, list]) => {
      const center = moduleCenters.get(module) || { x: centerX, y: centerY };
      list.forEach((item, index) => {
        const ring = Math.floor(index / 10);
        const offsetAngle = (Math.PI * 2 * (index % 10)) / 10;
        const orbit = 26 + ring * 24;
        positionMap.set(String(item.id), {
          x: center.x + Math.cos(offsetAngle) * orbit,
          y: center.y + Math.sin(offsetAngle) * orbit,
        });
      });
    });
    const maxDegree = Math.max(...rankedNodes.map((item) => Number(item.degree) || 0), 1);
    const edgeMarkup = rankedEdges.map((edge) => {
      const source = positionMap.get(String(edge.source));
      const target = positionMap.get(String(edge.target));
      if (!source || !target) return "";
      const payload = encodeURIComponent(JSON.stringify(edge));
      const stroke = String(edge.label || "+") === "-" ? "#b56576" : "#90be6d";
      const widthValue = 0.6 + Math.min(3.2, Math.abs(Number(edge.weight) || 0) * 2.4);
      return `<line class="community-network-edge" data-network-edge="${payload}" x1="${source.x}" y1="${source.y}" x2="${target.x}" y2="${target.y}" stroke="${escapeHtml(stroke)}" stroke-width="${widthValue}" stroke-opacity="0.32"></line>`;
    }).join("");
    const moduleMarkup = grouped.map(([module, list]) => {
      const center = moduleCenters.get(module) || { x: centerX, y: centerY };
      const radius = 42 + Math.min(120, list.length * 3.5);
      return `
        <g class="community-network-module-shell">
          <circle cx="${center.x}" cy="${center.y}" r="${radius}" fill="${escapeHtml(colorMap.get(module) || "#355070")}" fill-opacity="0.05" stroke="${escapeHtml(colorMap.get(module) || "#355070")}" stroke-opacity="0.28" stroke-width="1.5"></circle>
          <text x="${center.x}" y="${center.y - radius - 8}" class="chart-axis-label" text-anchor="middle">${escapeHtml(module)}</text>
        </g>
      `;
    }).join("");
    const nodeMarkup = rankedNodes.map((item) => {
      const position = positionMap.get(String(item.id));
      if (!position) return "";
      const radius = 5 + ((Number(item.degree) || 0) / maxDegree) * 10;
      const payload = encodeURIComponent(JSON.stringify(item));
      return `
        <g class="community-network-node-group" data-network-node="${payload}" tabindex="0">
          <circle class="community-network-node" cx="${position.x}" cy="${position.y}" r="${radius}" fill="${escapeHtml(colorMap.get(String(item.module || "未分模块")) || "#355070")}" fill-opacity="0.92" stroke="#fff" stroke-width="1.5"></circle>
        </g>
      `;
    }).join("");
    svg.innerHTML = `${moduleMarkup}${edgeMarkup}${nodeMarkup}`;

    const nodeElements = Array.from(container.querySelectorAll(".community-network-node-group"));
    const edgeElements = Array.from(container.querySelectorAll(".community-network-edge"));
    const resetHighlight = () => {
      tooltip.hidden = true;
      nodeElements.forEach((node) => node.classList.remove("is-active", "is-dimmed"));
      edgeElements.forEach((edge) => edge.classList.remove("is-active", "is-dimmed"));
    };
    const showNodeTooltip = (event, nodeElement) => {
      let data = null;
      try {
        data = JSON.parse(decodeURIComponent(nodeElement.getAttribute("data-network-node") || ""));
      } catch (_error) {
        data = null;
      }
      if (!data) return;
      const currentId = String(data.id || "");
      nodeElements.forEach((node) => {
        let nodeData = null;
        try {
          nodeData = JSON.parse(decodeURIComponent(node.getAttribute("data-network-node") || ""));
        } catch (_error) {
          nodeData = null;
        }
        const linked = !!nodeData && rankedEdges.some((edge) => {
          const source = String(edge.source || "");
          const target = String(edge.target || "");
          return (source === currentId && target === String(nodeData.id || "")) || (target === currentId && source === String(nodeData.id || ""));
        });
        node.classList.toggle("is-active", String(nodeData?.id || "") === currentId);
        node.classList.toggle("is-dimmed", !linked && String(nodeData?.id || "") !== currentId);
      });
      edgeElements.forEach((edge) => {
        let edgeData = null;
        try {
          edgeData = JSON.parse(decodeURIComponent(edge.getAttribute("data-network-edge") || ""));
        } catch (_error) {
          edgeData = null;
        }
        const linked = !!edgeData && (String(edgeData.source || "") === currentId || String(edgeData.target || "") === currentId);
        edge.classList.toggle("is-active", linked);
        edge.classList.toggle("is-dimmed", !linked);
      });
      tooltip.innerHTML = `
        <strong>${escapeHtml(data.label || data.id || "--")}</strong>
        <div>模块: ${escapeHtml(data.module || "未分模块")}</div>
        <div>门: ${escapeHtml(data.phylum || "未注释")}</div>
        <div>Degree: ${escapeHtml(String(data.degree ?? "--"))}</div>
        <div>丰度: ${escapeHtml(data.abundance != null ? `${data.abundance}%` : "--")}</div>
      `;
      tooltip.hidden = false;
      const rect = container.getBoundingClientRect();
      tooltip.style.left = `${Math.min(event.clientX - rect.left + 16, rect.width - 220)}px`;
      tooltip.style.top = `${Math.max(event.clientY - rect.top - 18, 12)}px`;
    };
    nodeElements.forEach((nodeElement) => {
      nodeElement.addEventListener("mouseenter", (event) => showNodeTooltip(event, nodeElement));
      nodeElement.addEventListener("mousemove", (event) => showNodeTooltip(event, nodeElement));
      nodeElement.addEventListener("mouseleave", resetHighlight);
      nodeElement.addEventListener("focus", () => {
        const rect = nodeElement.getBoundingClientRect();
        showNodeTooltip({ clientX: rect.left, clientY: rect.top }, nodeElement);
      });
      nodeElement.addEventListener("blur", resetHighlight);
    });
  };
  select.addEventListener("change", render);
  moduleSelect.addEventListener("change", render);
  render();
}

function renderCommunityNetworkHubBars(containerId, nodes = []) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const items = Array.isArray(nodes)
    ? [...nodes].sort((left, right) => (Number(right.degree) || 0) - (Number(left.degree) || 0)).slice(0, 20)
    : [];
  if (!items.length) {
    container.classList.add("empty-box");
    container.innerHTML = "<p>当前没有可展示的 hub taxa。</p>";
    return;
  }
  container.classList.remove("empty-box");
  const palette = buildCommunityNetworkPalette(Array.from(new Set(items.map((item) => String(item.module || "未分模块")))));
  const maxValue = Math.max(...items.map((item) => Number(item.degree) || 0), 1);
  container.innerHTML = `
    <div class="chart-insight" role="note" aria-label="Hub taxa 说明">
      <span class="chart-insight-label">结果说明</span>
      <p>这里按 Degree 排序查看潜在 hub taxa。颜色继承模块，悬浮可以同时看到门分类和相对丰度。</p>
    </div>
    <div class="mini-chart-card relation-chart-card community-biomarker-card">
      <span class="mini-chart-title">Top hub taxa</span>
      <div class="chart-canvas community-biomarker-canvas" style="--chart-height:440px">
        <svg class="sparkline-svg community-network-hub-svg" viewBox="0 0 1120 440" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Top hub taxa"></svg>
        <div class="chart-tooltip" hidden></div>
      </div>
    </div>
  `;
  const svg = container.querySelector(".community-network-hub-svg");
  const tooltip = container.querySelector(".chart-tooltip");
  if (!(svg instanceof SVGElement) || !(tooltip instanceof HTMLElement)) return;
  const width = 1120;
  const height = 440;
  const padLeft = 280;
  const padRight = 36;
  const padTop = 24;
  const padBottom = 28;
  const innerWidth = width - padLeft - padRight;
  const rowHeight = Math.max(24, Math.floor((height - padTop - padBottom) / Math.max(items.length, 1)));
  const barHeight = Math.min(18, Math.max(12, rowHeight - 8));
  svg.innerHTML = `
    ${items.map((item, index) => {
      const y = padTop + index * rowHeight + (rowHeight - barHeight) / 2;
      const value = Number(item.degree) || 0;
      const payload = encodeURIComponent(JSON.stringify(item));
      return `
        <g class="community-biomarker-bar-group" data-community-biomarker="${payload}" tabindex="0">
          <text class="chart-axis-label y-axis" x="${padLeft - 12}" y="${y + barHeight / 2 + 4}" text-anchor="end">${escapeHtml(item.label || item.id || "--")}</text>
          <rect class="community-biomarker-bar" x="${padLeft}" y="${y}" width="${(value / maxValue) * innerWidth}" height="${barHeight}" rx="8" ry="8" fill="${escapeHtml(palette.get(String(item.module || "未分模块")) || "#355070")}"></rect>
          <text class="chart-bar-value" x="${padLeft + ((value / maxValue) * innerWidth) + 8}" y="${y + barHeight / 2 + 4}">${escapeHtml(String(value))}</text>
        </g>
      `;
    }).join("")}
  `;
  const bars = Array.from(container.querySelectorAll(".community-biomarker-bar-group"));
  const hideTooltip = () => {
    tooltip.hidden = true;
    bars.forEach((bar) => bar.classList.remove("is-active", "is-dimmed"));
  };
  const showTooltip = (event, bar) => {
    let data = null;
    try {
      data = JSON.parse(decodeURIComponent(bar.getAttribute("data-community-biomarker") || ""));
    } catch (_error) {
      data = null;
    }
    if (!data) return;
    bars.forEach((item) => item.classList.toggle("is-dimmed", item !== bar));
    bar.classList.add("is-active");
    tooltip.innerHTML = `
      <strong>${escapeHtml(data.label || data.id || "--")}</strong>
      <div>模块: ${escapeHtml(data.module || "未分模块")}</div>
      <div>门: ${escapeHtml(data.phylum || "未注释")}</div>
      <div>Degree: ${escapeHtml(String(data.degree ?? "--"))}</div>
      <div>丰度: ${escapeHtml(data.abundance != null ? `${data.abundance}%` : "--")}</div>
    `;
    tooltip.hidden = false;
    const rect = container.getBoundingClientRect();
    tooltip.style.left = `${Math.min(event.clientX - rect.left + 16, rect.width - 220)}px`;
    tooltip.style.top = `${Math.max(event.clientY - rect.top - 18, 12)}px`;
  };
  bars.forEach((bar) => {
    bar.addEventListener("mouseenter", (event) => showTooltip(event, bar));
    bar.addEventListener("mousemove", (event) => showTooltip(event, bar));
    bar.addEventListener("mouseleave", hideTooltip);
    bar.addEventListener("focus", () => {
      const rect = bar.getBoundingClientRect();
      showTooltip({ clientX: rect.left, clientY: rect.top }, bar);
    });
    bar.addEventListener("blur", hideTooltip);
  });
}

function renderCommunityNetworkRoleScatter(containerId, nodes = []) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const items = Array.isArray(nodes)
    ? nodes.filter((item) => Number.isFinite(Number(item?.p)) && Number.isFinite(Number(item?.z)))
    : [];
  if (!items.length) {
    container.classList.add("empty-box");
    container.innerHTML = "<p>当前没有可用于绘制 Zi-Pi 角色图的节点。</p>";
    return;
  }
  container.classList.remove("empty-box");
  const rolePalette = new Map([
    ["Module hubs", "#6c63b5"],
    ["Connectors", "#d95f02"],
    ["Peripheral nodes", "#1b9e77"],
    ["Network hubs", "#7b2cbf"],
    ["未分类", "#8c8c8c"],
  ]);
  const roles = Array.from(new Set(items.map((item) => String(item.role || "未分类"))));
  container.innerHTML = `
    <div class="chart-insight" role="note" aria-label="节点角色散点图说明">
      <span class="chart-insight-label">结果说明</span>
      <p>这张图直接复现 microeco 常见的 Zi-Pi 角色判定方式。横轴是 among-module connectivity，纵轴是 within-module connectivity，虚线分别对应 P=0.62 和 Z=2.5。</p>
    </div>
    <div class="mini-chart-card relation-chart-card community-biomarker-card">
      <span class="mini-chart-title">Node roles (Zi-Pi plot)</span>
      <div class="chart-canvas community-biomarker-canvas" style="--chart-height:520px">
        <svg class="sparkline-svg community-network-role-scatter-svg" viewBox="0 0 1120 520" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Node roles scatter"></svg>
        <div class="chart-tooltip" hidden></div>
      </div>
      <div class="community-pcoa-legend">
        ${roles.map((role) => `
          <span class="community-pcoa-legend-item">
            <i style="background:${escapeHtml(rolePalette.get(role) || "#8c8c8c")}"></i>
            <span>${escapeHtml(role)}</span>
          </span>
        `).join("")}
      </div>
    </div>
  `;
  const svg = container.querySelector(".community-network-role-scatter-svg");
  const tooltip = container.querySelector(".chart-tooltip");
  if (!(svg instanceof SVGElement) || !(tooltip instanceof HTMLElement)) return;
  const width = 1120;
  const height = 520;
  const padLeft = 92;
  const padRight = 40;
  const padTop = 24;
  const padBottom = 58;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const xMax = Math.max(1, ...items.map((item) => Number(item.p) || 0));
  const yMax = Math.max(5.5, ...items.map((item) => Number(item.z) || 0));
  const xScale = (value) => padLeft + (Math.max(0, Math.min(xMax, value)) / xMax) * innerWidth;
  const yScale = (value) => padTop + innerHeight - (Math.max(-1.5, Math.min(yMax, value)) + 1.5) / (yMax + 1.5) * innerHeight;
  const roleThresholdX = xScale(0.62);
  const roleThresholdY = yScale(2.5);
  const payload = (item) => encodeURIComponent(JSON.stringify(item));
  svg.innerHTML = `
    <rect x="${padLeft}" y="${padTop}" width="${innerWidth}" height="${innerHeight}" fill="rgba(255,255,255,0.88)" stroke="rgba(62,84,111,0.16)"></rect>
    <line class="chart-axis-line" x1="${padLeft}" y1="${padTop + innerHeight}" x2="${width - padRight}" y2="${padTop + innerHeight}"></line>
    <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${padTop + innerHeight}"></line>
    <line x1="${roleThresholdX}" y1="${padTop}" x2="${roleThresholdX}" y2="${padTop + innerHeight}" stroke="#20232a" stroke-width="2.2" stroke-dasharray="12 10"></line>
    <line x1="${padLeft}" y1="${roleThresholdY}" x2="${width - padRight}" y2="${roleThresholdY}" stroke="#20232a" stroke-width="2.2" stroke-dasharray="12 10"></line>
    ${[0, 0.25, 0.5, 0.75, 1].map((tick) => `
      <text class="chart-axis-label x-axis" x="${xScale(tick)}" y="${height - 18}" text-anchor="middle">${escapeHtml(tick.toFixed(2))}</text>
    `).join("")}
    ${[0, 2, 4].map((tick) => `
      <text class="chart-axis-label y-axis" x="${padLeft - 12}" y="${yScale(tick) + 4}" text-anchor="end">${escapeHtml(String(tick))}</text>
    `).join("")}
    <text class="chart-axis-title" x="${padLeft + innerWidth / 2}" y="${height - 6}" text-anchor="middle">Among-module connectivity</text>
    <text class="chart-axis-title" transform="translate(28 ${padTop + innerHeight / 2}) rotate(-90)" text-anchor="middle">Within-module connectivity</text>
    ${items.map((item) => `
      <circle class="community-network-point" data-network-role-point="${payload(item)}" cx="${xScale(Number(item.p) || 0)}" cy="${yScale(Number(item.z) || 0)}" r="6" fill="${escapeHtml(rolePalette.get(String(item.role || "未分类")) || "#8c8c8c")}" fill-opacity="0.92"></circle>
    `).join("")}
  `;
  const points = Array.from(container.querySelectorAll("[data-network-role-point]"));
  const hideTooltip = () => {
    tooltip.hidden = true;
    points.forEach((point) => point.classList.remove("is-active", "is-dimmed"));
  };
  const showTooltip = (event, point) => {
    let data = null;
    try {
      data = JSON.parse(decodeURIComponent(point.getAttribute("data-network-role-point") || ""));
    } catch (_error) {
      data = null;
    }
    if (!data) return;
    points.forEach((item) => item.classList.toggle("is-dimmed", item !== point));
    point.classList.add("is-active");
    tooltip.innerHTML = `
      <strong>${escapeHtml(data.label || data.id || "--")}</strong>
      <div>Role: ${escapeHtml(data.role || "未分类")}</div>
      <div>Genus: ${escapeHtml(data.genus || "未注释")}</div>
      <div>Pi: ${escapeHtml(String(data.p ?? "--"))}</div>
      <div>Zi: ${escapeHtml(String(data.z ?? "--"))}</div>
    `;
    tooltip.hidden = false;
    const rect = container.getBoundingClientRect();
    tooltip.style.left = `${Math.min(event.clientX - rect.left + 16, rect.width - 220)}px`;
    tooltip.style.top = `${Math.max(event.clientY - rect.top - 18, 12)}px`;
  };
  points.forEach((point) => {
    point.addEventListener("mouseenter", (event) => showTooltip(event, point));
    point.addEventListener("mousemove", (event) => showTooltip(event, point));
    point.addEventListener("mouseleave", hideTooltip);
  });
}

function renderCommunityNetworkLayeredRoles(containerId, nodes = []) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const filtered = Array.isArray(nodes)
    ? nodes.filter((item) => Number.isFinite(Number(item?.p)) && Number.isFinite(Number(item?.z)) && String(item?.genus || "").trim())
    : [];
  if (!filtered.length) {
    container.classList.add("empty-box");
    container.innerHTML = "<p>当前没有可用于绘制 Genus 分层角色图的节点。</p>";
    return;
  }
  const genusRanked = Array.from(
    filtered.reduce((map, item) => {
      const genus = String(item.genus || "未注释").trim() || "未注释";
      const current = map.get(genus) || 0;
      map.set(genus, current + (Number(item.abundance) || 0));
      return map;
    }, new Map()).entries()
  ).sort((left, right) => right[1] - left[1]).slice(0, 12).map(([genus]) => genus);
  const items = filtered
    .filter((item) => genusRanked.includes(String(item.genus || "").trim() || "未注释"))
    .sort((left, right) => genusRanked.indexOf(String(left.genus || "").trim()) - genusRanked.indexOf(String(right.genus || "").trim()));
  const formatGenusName = (value) => {
    const raw = String(value || "").trim() || "未注释";
    const stripped = raw.replace(/^[a-z]__+/i, "").trim() || "未注释";
    return stripped;
  };
  const compactGenusLabel = (value) => {
    const clean = formatGenusName(value);
    if (clean.length <= 16) return clean;
    const segments = clean.split(/[-_]/).filter(Boolean);
    if (segments.length >= 2) {
      const head = segments[0];
      const tail = segments[segments.length - 1];
      const composed = `${head.slice(0, 8)}...${tail.slice(0, 6)}`;
      if (composed.length <= 18) return composed;
    }
    return `${clean.slice(0, 12)}...`;
  };
  const genusPalette = buildCommunityNetworkPalette(genusRanked);
  const roleShape = new Map([
    ["Connectors", "circle"],
    ["Module hubs", "triangle"],
    ["Peripheral nodes", "square"],
    ["Network hubs", "diamond"],
    ["未分类", "circle"],
  ]);
  container.classList.remove("empty-box");
  container.innerHTML = `
    <div class="chart-insight" role="note" aria-label="Genus 分层角色图说明">
      <span class="chart-insight-label">结果说明</span>
      <p>这张图按 Genus 展开节点角色分布，上面一层是 among-module connectivity，下面一层是 within-module connectivity。点大小表示丰度，形状表示角色。</p>
    </div>
    <div class="mini-chart-card relation-chart-card community-biomarker-card">
      <span class="mini-chart-title">Layered node roles by Genus</span>
      <div class="chart-canvas community-biomarker-canvas" style="--chart-height:620px">
        <svg class="sparkline-svg community-network-layered-svg" viewBox="0 0 1180 620" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Layered node roles by genus"></svg>
        <div class="chart-tooltip" hidden></div>
      </div>
    </div>
  `;
  const svg = container.querySelector(".community-network-layered-svg");
  const tooltip = container.querySelector(".chart-tooltip");
  if (!(svg instanceof SVGElement) || !(tooltip instanceof HTMLElement)) return;
  const width = 1180;
  const height = 620;
  const padLeft = 92;
  const padRight = 28;
  const padTop = 18;
  const padBottom = 104;
  const panelGap = 22;
  const panelHeight = (height - padTop - padBottom - panelGap) / 2;
  const panelTopA = padTop;
  const panelTopB = padTop + panelHeight + panelGap;
  const xStep = (width - padLeft - padRight) / Math.max(genusRanked.length, 1);
  const pMax = Math.max(1, ...items.map((item) => Number(item.p) || 0));
  const zMax = Math.max(5.5, ...items.map((item) => Number(item.z) || 0));
  const abundanceMax = Math.max(1, ...items.map((item) => Number(item.abundance) || 0));
  const xFor = (genus, offset) => padLeft + genusRanked.indexOf(genus) * xStep + xStep / 2 + offset;
  const yForP = (value) => panelTopA + panelHeight - (Math.max(0, Math.min(pMax, value)) / pMax) * (panelHeight - 22);
  const yForZ = (value) => panelTopB + panelHeight - ((Math.max(-1.5, Math.min(zMax, value)) + 1.5) / (zMax + 1.5)) * (panelHeight - 22);
  const offsetCount = {};
  const shapeMarkup = (shape, x, y, size, color, payload) => {
    if (shape === "triangle") {
      return `<path class="community-network-point" data-network-role-point="${payload}" d="M ${x} ${y - size} L ${x - size} ${y + size} L ${x + size} ${y + size} Z" fill="${color}" fill-opacity="0.84"></path>`;
    }
    if (shape === "square") {
      return `<rect class="community-network-point" data-network-role-point="${payload}" x="${x - size}" y="${y - size}" width="${size * 2}" height="${size * 2}" fill="${color}" fill-opacity="0.16" stroke="${color}" stroke-width="1.4"></rect><path class="community-network-point-overlay" data-network-role-point="${payload}" d="M ${x - size + 2} ${y - size + 2} L ${x + size - 2} ${y + size - 2} M ${x + size - 2} ${y - size + 2} L ${x - size + 2} ${y + size - 2}" stroke="${color}" stroke-width="1.1" fill="none"></path>`;
    }
    if (shape === "diamond") {
      return `<path class="community-network-point" data-network-role-point="${payload}" d="M ${x} ${y - size} L ${x - size} ${y} L ${x} ${y + size} L ${x + size} ${y} Z" fill="${color}" fill-opacity="0.84"></path>`;
    }
    return `<circle class="community-network-point" data-network-role-point="${payload}" cx="${x}" cy="${y}" r="${size}" fill="${color}" fill-opacity="0.84"></circle>`;
  };
  svg.innerHTML = `
    <rect x="${padLeft}" y="${panelTopA}" width="${width - padLeft - padRight}" height="${panelHeight}" fill="rgba(255,255,255,0.92)" stroke="rgba(62,84,111,0.12)"></rect>
    <rect x="${padLeft}" y="${panelTopB}" width="${width - padLeft - padRight}" height="${panelHeight}" fill="rgba(255,255,255,0.92)" stroke="rgba(62,84,111,0.12)"></rect>
    <text class="chart-axis-title" transform="translate(${width - 16} ${panelTopA + panelHeight / 2}) rotate(90)" text-anchor="middle">Among-module connectivity</text>
    <text class="chart-axis-title" transform="translate(${width - 16} ${panelTopB + panelHeight / 2}) rotate(90)" text-anchor="middle">Within-module connectivity</text>
    ${genusRanked.map((genus, index) => `
      <line class="chart-grid-line" x1="${padLeft + index * xStep}" y1="${panelTopA}" x2="${padLeft + index * xStep}" y2="${panelTopB + panelHeight}"></line>
      <text class="chart-axis-label x-axis" transform="translate(${padLeft + index * xStep + xStep / 2} ${height - 20}) rotate(-36)" text-anchor="end">${escapeHtml(compactGenusLabel(genus))}</text>
    `).join("")}
    <text class="chart-axis-title" x="${padLeft + (width - padLeft - padRight) / 2}" y="${height - 4}" text-anchor="middle">Genus</text>
    ${items.map((item, index) => {
      const genus = String(item.genus || "未注释").trim() || "未注释";
      const key = `${genus}-${index % 7}`;
      const offsetIndex = offsetCount[key] || 0;
      offsetCount[key] = offsetIndex + 1;
      const offset = ((offsetIndex % 7) - 3) * 8;
      const size = 3 + ((Number(item.abundance) || 0) / abundanceMax) * 9;
      const color = genusPalette.get(genus) || "#355070";
      const shape = roleShape.get(String(item.role || "未分类")) || "circle";
      const payload = encodeURIComponent(JSON.stringify(item));
      return `
        ${shapeMarkup(shape, xFor(genus, offset), yForP(Number(item.p) || 0), size, color, payload)}
        ${shapeMarkup(shape, xFor(genus, offset), yForZ(Number(item.z) || 0), size, color, payload)}
      `;
    }).join("")}
  `;
  const points = Array.from(container.querySelectorAll("[data-network-role-point]"));
  const hideTooltip = () => {
    tooltip.hidden = true;
    points.forEach((point) => point.classList.remove("is-active", "is-dimmed"));
  };
  const showTooltip = (event, point) => {
    let data = null;
    try {
      data = JSON.parse(decodeURIComponent(point.getAttribute("data-network-role-point") || ""));
    } catch (_error) {
      data = null;
    }
    if (!data) return;
    points.forEach((item) => item.classList.toggle("is-dimmed", item !== point));
    point.classList.add("is-active");
    tooltip.innerHTML = `
      <strong>${escapeHtml(data.label || data.id || "--")}</strong>
      <div>Genus: ${escapeHtml(formatGenusName(data.genus || "未注释"))}</div>
      <div>Role: ${escapeHtml(data.role || "未分类")}</div>
      <div>Abundance: ${escapeHtml(data.abundance != null ? `${data.abundance}%` : "--")}</div>
      <div>Pi / Zi: ${escapeHtml(String(data.p ?? "--"))} / ${escapeHtml(String(data.z ?? "--"))}</div>
    `;
    tooltip.hidden = false;
    const rect = container.getBoundingClientRect();
    tooltip.style.left = `${Math.min(event.clientX - rect.left + 16, rect.width - 240)}px`;
    tooltip.style.top = `${Math.max(event.clientY - rect.top - 18, 12)}px`;
  };
  points.forEach((point) => {
    point.addEventListener("mouseenter", (event) => showTooltip(event, point));
    point.addEventListener("mousemove", (event) => showTooltip(event, point));
    point.addEventListener("mouseleave", hideTooltip);
  });
}

function clusterCommunityDistanceTree(matrix) {
  const size = Array.isArray(matrix) ? matrix.length : 0;
  if (!size) return null;
  let clusters = Array.from({ length: size }, (_, index) => ({
    members: [index],
    leaf: index,
    height: 0,
  }));
  const averageDistance = (leftMembers, rightMembers) => {
    let total = 0;
    let count = 0;
    leftMembers.forEach((left) => {
      rightMembers.forEach((right) => {
        if (left === right) return;
        total += Number(matrix[left]?.[right]) || 0;
        count += 1;
      });
    });
    return count ? total / count : 0;
  };
  while (clusters.length > 1) {
    let bestLeft = 0;
    let bestRight = 1;
    let bestDistance = Number.POSITIVE_INFINITY;
    for (let left = 0; left < clusters.length; left += 1) {
      for (let right = left + 1; right < clusters.length; right += 1) {
        const distance = averageDistance(clusters[left].members, clusters[right].members);
        if (distance < bestDistance) {
          bestDistance = distance;
          bestLeft = left;
          bestRight = right;
        }
      }
    }
    const merged = {
      members: [...clusters[bestLeft].members, ...clusters[bestRight].members],
      left: clusters[bestLeft],
      right: clusters[bestRight],
      height: bestDistance,
    };
    if ((merged.left?.members?.[0] ?? 0) > (merged.right?.members?.[0] ?? 0)) {
      merged.left = clusters[bestRight];
      merged.right = clusters[bestLeft];
    }
    clusters = clusters.filter((_, index) => index !== bestLeft && index !== bestRight);
    clusters.push(merged);
  }
  return clusters[0] || null;
}

function collectCommunityClusterOrder(node) {
  if (!node) return [];
  if (typeof node.leaf === "number") return [node.leaf];
  return [
    ...collectCommunityClusterOrder(node.left),
    ...collectCommunityClusterOrder(node.right),
  ];
}

function layoutCommunityCluster(node, leafCenters, maxHeight, orientation = "top") {
  if (!node) return null;
  if (typeof node.leaf === "number") {
    const center = leafCenters[node.leaf] || 0;
    return {
      ...node,
      x: orientation === "top" ? center : 0,
      y: orientation === "top" ? 0 : center,
    };
  }
  const left = layoutCommunityCluster(node.left, leafCenters, maxHeight, orientation);
  const right = layoutCommunityCluster(node.right, leafCenters, maxHeight, orientation);
  const scaled = maxHeight > 0 ? (node.height / maxHeight) : 0;
  if (orientation === "top") {
    return {
      ...node,
      left,
      right,
      x: ((left?.x || 0) + (right?.x || 0)) / 2,
      y: scaled,
    };
  }
  return {
    ...node,
    left,
    right,
    x: scaled,
    y: ((left?.y || 0) + (right?.y || 0)) / 2,
  };
}

function renderCommunityClusterBranches(node, width, height, orientation = "top") {
  if (!node || typeof node.leaf === "number") return "";
  const left = renderCommunityClusterBranches(node.left, width, height, orientation);
  const right = renderCommunityClusterBranches(node.right, width, height, orientation);
  if (orientation === "top") {
    const x1 = node.left?.x || 0;
    const x2 = node.right?.x || 0;
    const y = height - (node.y || 0) * height;
    const y1 = height - ((node.left?.y || 0) * height);
    const y2 = height - ((node.right?.y || 0) * height);
    return `
      ${left}
      ${right}
      <path class="community-cluster-branch" d="M ${x1} ${y1} L ${x1} ${y} L ${x2} ${y} L ${x2} ${y2}"></path>
    `;
  }
  const y1 = node.left?.y || 0;
  const y2 = node.right?.y || 0;
  const x = width - (node.x || 0) * width;
  const x1 = width - ((node.left?.x || 0) * width);
  const x2 = width - ((node.right?.x || 0) * width);
  return `
    ${left}
    ${right}
    <path class="community-cluster-branch" d="M ${x1} ${y1} L ${x} ${y1} L ${x} ${y2} L ${x2} ${y2}"></path>
  `;
}

function renderCommunityBetaHeatmap(containerId, section) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const matrixSection = section?.distance_matrix || {};
  const samples = Array.isArray(matrixSection?.samples) ? matrixSection.samples : [];
  const rows = Array.isArray(matrixSection?.rows) ? matrixSection.rows : [];
  const groups = matrixSection?.groups && typeof matrixSection.groups === "object" ? matrixSection.groups : {};
  if (!samples.length || !rows.length || rows.length !== samples.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有可展示的距离矩阵热图</strong>
        <p class="empty-copy">请先生成 bray_distance_matrix.tsv 后再查看交互热图。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const min = Number(matrixSection?.min ?? 0);
  const max = Number(matrixSection?.max ?? 1);
  const tree = clusterCommunityDistanceTree(rows);
  const order = collectCommunityClusterOrder(tree);
  const orderedSamples = order.length === samples.length ? order.map((index) => samples[index]) : samples.slice();
  const orderedRows = (order.length === samples.length ? order : samples.map((_, index) => index)).map((rowIndex) => ({
    sample: samples[rowIndex],
    values: (order.length === samples.length ? order : samples.map((_, index) => index)).map((colIndex) => {
      const value = Number(rows[rowIndex]?.[colIndex] ?? 0);
      return {
        source: samples[rowIndex],
        target: samples[colIndex],
        value,
        display: formatChartValue(value),
      };
    }),
  }));
  const defaultCell = orderedRows[0]?.values?.[0] || null;
  const defaultDetail = defaultCell
    ? `${defaultCell.source} (${groups[defaultCell.source] || "未分组"}) vs ${defaultCell.target} (${groups[defaultCell.target] || "未分组"}): ${defaultCell.display}`
    : "悬停任意格点查看样本对距离。";
  const significanceText = section?.summary?.permanova_p != null
    ? `PERMANOVA p = ${section.summary.permanova_p}`
    : "PERMANOVA 未提供";
  const groupNames = Array.from(new Set(orderedSamples.map((sample) => String(groups[sample] || "未分组"))));
  const palette = ["#3e546f", "#8a6654", "#7a7158", "#4f7f6b", "#7c5d83", "#9a6b62", "#557a95", "#a88262", "#5b6b4d"];
  const colorMap = new Map(groupNames.map((group, index) => [group, palette[index % palette.length]]));
  const viewportWidth = Math.max(window.innerWidth || 1280, 960);
  const viewportHeight = Math.max(window.innerHeight || 860, 720);
  const heatmapCount = orderedSamples.length;
  const containerWidth = Math.max(container.clientWidth || 0, Math.min(viewportWidth - 120, 1820));
  const availableWidth = Math.max(560, Math.min(containerWidth - 24, 1820));
  const availableHeight = Math.max(420, Math.min(viewportHeight - 260, 1040));
  const annotationThickness = heatmapCount > 64 ? 8 : 10;
  const plotSize = Math.max(420, Math.min(availableWidth, availableHeight));
  const treeSize = Math.max(46, Math.min(92, Math.round(plotSize * 0.12)));
  const leftTreeWidth = treeSize;
  const topTreeHeight = treeSize;
  const heatmapSide = Math.max(320, plotSize - treeSize - annotationThickness);
  const cellSize = Math.max(6, Math.min(18, Math.floor(heatmapSide / Math.max(heatmapCount, 1))));
  const heatmapBodyWidth = heatmapCount * cellSize;
  const heatmapBodyHeight = heatmapCount * cellSize;
  const rowHeaderWidth = annotationThickness;
  const headerHeight = annotationThickness;
  const compactClass = cellSize <= 10 ? " community-beta-heatmap-composite-compact" : "";
  const orderedIndices = order.length === samples.length ? order : samples.map((_, index) => index);
  const topLeafCenters = Object.fromEntries(orderedIndices.map((sampleIndex, position) => [sampleIndex, position * cellSize + cellSize / 2]));
  const leftLeafCenters = Object.fromEntries(orderedIndices.map((sampleIndex, position) => [sampleIndex, position * cellSize + cellSize / 2]));
  const maxTreeHeight = Number(tree?.height || 0);
  const topTree = layoutCommunityCluster(tree, topLeafCenters, maxTreeHeight, "top");
  const leftTree = layoutCommunityCluster(tree, leftLeafCenters, maxTreeHeight, "left");
  const topTreeSvg = renderCommunityClusterBranches(topTree, heatmapBodyWidth, topTreeHeight, "top");
  const leftTreeSvg = renderCommunityClusterBranches(leftTree, leftTreeWidth, heatmapBodyHeight, "left");
  const compositeWidth = leftTreeWidth + rowHeaderWidth + heatmapBodyWidth;
  const compositeHeight = topTreeHeight + headerHeight + heatmapBodyHeight;
  container.innerHTML = `
    <div class="patho-heatmap-card community-beta-heatmap-card">
      <div class="patho-heatmap-head">
        <div>
          <strong>Bray-Curtis 距离矩阵聚类热图</strong>
          <p>横轴和纵轴都按样本间距离进行层次聚类排序；悬停查看样本对距离，点击后可锁定当前说明。</p>
        </div>
        <div class="patho-heatmap-scale" aria-hidden="true">
          <span>${escapeHtml(String(min))}</span>
          <div class="patho-heatmap-scale-bar patho-heatmap-scale-bar-208"></div>
          <span>${escapeHtml(String(max))}</span>
        </div>
      </div>
      <div class="community-beta-heatmap-legend">
        <span class="community-significance-badge ${(section?.summary?.permanova_p != null && Number(section.summary.permanova_p) < 0.05) ? "is-significant" : "is-neutral"}">${escapeHtml(significanceText)}</span>
        ${groupNames.map((group) => `
          <span class="community-pcoa-legend-item">
            <i style="background:${escapeHtml(colorMap.get(group) || "#3e546f")}"></i>
            <span>${escapeHtml(group)}</span>
          </span>
        `).join("")}
      </div>
      <div class="patho-heatmap-detail" data-heatmap-detail>${escapeHtml(defaultDetail)}</div>
      <div class="community-beta-heatmap-frame">
        <div
          class="community-beta-heatmap-composite${compactClass}"
          style="width:${compositeWidth}px; height:${compositeHeight}px; grid-template-columns:${leftTreeWidth}px ${rowHeaderWidth}px repeat(${orderedSamples.length}, ${cellSize}px); grid-template-rows:${topTreeHeight}px ${headerHeight}px repeat(${orderedSamples.length}, ${cellSize}px);"
        >
          <div class="community-beta-heatmap-top-tree" style="grid-column: 3 / span ${orderedSamples.length}; grid-row: 1; width:${heatmapBodyWidth}px; height:${topTreeHeight}px;">
            <svg viewBox="0 0 ${heatmapBodyWidth} ${topTreeHeight}" preserveAspectRatio="none" role="img" aria-label="顶部层次聚类树">
              ${topTreeSvg}
            </svg>
          </div>
          <div class="community-beta-heatmap-left-tree" style="grid-column: 1; grid-row: 3 / span ${orderedSamples.length}; width:${leftTreeWidth}px; height:${heatmapBodyHeight}px;">
            <svg viewBox="0 0 ${leftTreeWidth} ${heatmapBodyHeight}" preserveAspectRatio="none" role="img" aria-label="左侧层次聚类树">
              ${leftTreeSvg}
            </svg>
          </div>
          <div class="community-beta-heatmap-corner community-beta-heatmap-annotation-corner" style="grid-column: 2; grid-row: 2;" aria-hidden="true"></div>
          ${orderedSamples.map((sample, sampleIndex) => `
            <div
              class="community-beta-heatmap-col-header"
              style="grid-column: ${sampleIndex + 3}; grid-row: 2;"
              data-patho-sample-name="${escapeHtml(normalizePathoSampleName(sample))}"
              title="${escapeHtml(`${sample} | ${groups[sample] || "未分组"}`)}"
              aria-label="${escapeHtml(`${sample}，分组：${groups[sample] || "未分组"}`)}"
            >
              <i class="community-beta-heatmap-group-dot" style="background:${escapeHtml(colorMap.get(String(groups[sample] || "未分组")) || "#3e546f")}"></i>
            </div>
          `).join("")}
          ${orderedRows.map((row, rowIndex) => `
            <div
              class="community-beta-heatmap-row-header"
              style="grid-column: 2; grid-row: ${rowIndex + 3};"
              data-patho-sample-name="${escapeHtml(normalizePathoSampleName(row.sample))}"
              title="${escapeHtml(`${row.sample} | ${groups[row.sample] || "未分组"}`)}"
              aria-label="${escapeHtml(`${row.sample}，分组：${groups[row.sample] || "未分组"}`)}"
            >
              <i class="community-beta-heatmap-group-dot" style="background:${escapeHtml(colorMap.get(String(groups[row.sample] || "未分组")) || "#3e546f")}"></i>
            </div>
            ${row.values.map((cell, colIndex) => {
              const isDiagonal = cell.source === cell.target;
              const bg = isDiagonal ? "rgba(62, 84, 111, 0.08)" : buildHeatValueColor(cell.value, min, max, 208);
              const color = isDiagonal ? "var(--report-ink-soft)" : buildHeatTextColor(cell.value, min, max);
              const detail = `${cell.source} (${groups[cell.source] || "未分组"}) vs ${cell.target} (${groups[cell.target] || "未分组"}): ${cell.display}`;
              return `
                <button
                  class="community-beta-heatmap-cell${isDiagonal ? " is-diagonal" : ""}"
                  type="button"
                  style="grid-column: ${colIndex + 3}; grid-row: ${rowIndex + 3}; --cell-bg:${bg}; --cell-color:${color};"
                  data-heatmap-detail-text="${escapeHtml(detail)}"
                  data-patho-sample-a="${escapeHtml(normalizePathoSampleName(cell.source))}"
                  data-patho-sample-b="${escapeHtml(normalizePathoSampleName(cell.target))}"
                  title="${escapeHtml(detail)}"
                ><span class="sr-only">${escapeHtml(cell.display)}</span></button>
              `;
            }).join("")}
          `).join("")}
        </div>
      </div>
    </div>
  `;
  const detailNode = container.querySelector("[data-heatmap-detail]");
  let lockedButton = null;
  container.querySelectorAll(".community-beta-heatmap-cell").forEach((button) => {
    button.addEventListener("mouseenter", () => {
      if (lockedButton) return;
      if (detailNode) detailNode.textContent = button.dataset.heatmapDetailText || defaultDetail;
    });
    button.addEventListener("focus", () => {
      if (lockedButton) return;
      if (detailNode) detailNode.textContent = button.dataset.heatmapDetailText || defaultDetail;
    });
    button.addEventListener("click", () => {
      if (lockedButton === button) {
        button.classList.remove("is-locked");
        lockedButton = null;
        if (detailNode) detailNode.textContent = defaultDetail;
        return;
      }
      if (lockedButton) lockedButton.classList.remove("is-locked");
      lockedButton = button;
      button.classList.add("is-locked");
      if (detailNode) detailNode.textContent = button.dataset.heatmapDetailText || defaultDetail;
    });
  });
}

function renderCommunityBetaClusterComposition(containerId, betaSection, taxaSection) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const matrixSection = betaSection?.distance_matrix || {};
  const samples = Array.isArray(matrixSection?.samples) ? matrixSection.samples : [];
  const rows = Array.isArray(matrixSection?.rows) ? matrixSection.rows : [];
  const groups = matrixSection?.groups && typeof matrixSection.groups === "object" ? matrixSection.groups : {};
  const levels = taxaSection?.levels && typeof taxaSection.levels === "object" ? taxaSection.levels : {};
  const rankOptions = ["目", "科", "属"].filter((level) => Array.isArray(levels?.[level]?.sample_series) && levels[level].sample_series.length);
  if (!samples.length || !rows.length || rows.length !== samples.length || !rankOptions.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有可展示的聚类组成图</strong>
        <p class="empty-copy">需要同时存在距离矩阵和目 / 科 / 属水平样本组成数据后，才能绘制聚类树与堆积条形图。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const tree = clusterCommunityDistanceTree(rows);
  const order = collectCommunityClusterOrder(tree);
  const orderedIndices = order.length === samples.length ? order : samples.map((_, index) => index);
  const orderedSamples = orderedIndices.map((index) => samples[index]);
  const viewportWidth = Math.max(window.innerWidth || 1280, 960);
  const containerWidth = Math.max(container.clientWidth || 0, Math.min(viewportWidth - 120, 1820));
  const availableWidth = Math.max(620, Math.min(containerWidth - 24, 1820));
  const treeWidth = Math.max(74, Math.min(118, Math.round(availableWidth * 0.13)));
  const stripWidth = 10;
  const maxTreeHeight = Number(tree?.height || 0);
  container.innerHTML = `
    <div class="community-beta-composition-card">
      <div class="community-beta-composition-head">
        <div>
          <strong>层次聚类树与样本组成</strong>
          <p>按 Bray-Curtis 层次聚类顺序排列样本，同时展示不同分类层级下的相对丰度堆积组成。</p>
        </div>
        <div class="community-beta-composition-controls">
          <label class="community-beta-control">
            <span>层级</span>
            <select data-community-beta-rank>
              ${rankOptions.map((level) => `<option value="${escapeHtml(level)}">${escapeHtml(level)}</option>`).join("")}
            </select>
          </label>
          <label class="community-beta-control">
            <span>TopN</span>
            <select data-community-beta-topn>
              ${[5, 8, 10, 12, 15].map((count) => `<option value="${count}"${count === 8 ? " selected" : ""}>Top ${count}</option>`).join("")}
            </select>
          </label>
          <span class="card-tag" data-community-beta-composition-tag>Bray + ${escapeHtml(rankOptions[rankOptions.length - 1] || rankOptions[0])}</span>
        </div>
      </div>
      <div class="community-beta-composition-frame" data-community-beta-composition-stage></div>
      <div class="community-beta-composition-legend" data-community-beta-composition-legend></div>
    </div>
  `;
  const rankSelect = container.querySelector("[data-community-beta-rank]");
  const topnSelect = container.querySelector("[data-community-beta-topn]");
  const stageNode = container.querySelector("[data-community-beta-composition-stage]");
  const legendNode = container.querySelector("[data-community-beta-composition-legend]");
  const tagNode = container.querySelector("[data-community-beta-composition-tag]");
  const groupPalette = ["#6c7f92", "#9d8578", "#8f9073", "#6f8f84", "#8a7694", "#aa8576", "#6f90a3", "#ad9672", "#7d836e"];
  const naturePalette = ["#5b7f78", "#c57b57", "#7e91b2", "#b56d72", "#8aa184", "#9e82b2", "#bea55c", "#6e95a2", "#a38772", "#d4ccc2"];
  const groupList = Array.from(new Set(rankOptions.flatMap((level) => (levels[level]?.sample_series || []).map((entry) => String(entry?.group || groups[entry?.sample] || "未分组")))));
  const groupColorMap = new Map(groupList.map((group, index) => [group, groupPalette[index % groupPalette.length]]));
  const state = {
    level: rankOptions.includes("属") ? "属" : rankOptions[0],
    topN: 8,
  };
  const isAnnotatedCommunityTaxon = (label) => {
    const normalized = String(label || "").trim();
    if (!normalized) return false;
    const canonical = normalized.replace(/^[^:]+:\s*/, "").trim();
    const compact = canonical.toLowerCase();
    return !(
      canonical === "未注释"
      || canonical === "未分类"
      || canonical === "其他"
      || compact === "unclassified"
      || compact === "unassigned"
      || compact === "unknown"
      || compact === "norank"
      || compact === "uncultured"
      || compact === "ambiguous_taxa"
      || compact.startsWith("unclassified ")
      || compact.startsWith("unknown ")
    );
  };

  const renderPlot = () => {
    const levelBlock = levels?.[state.level] || {};
    const sampleSeries = Array.isArray(levelBlock?.sample_series) ? levelBlock.sample_series : [];
    const seriesMap = new Map(sampleSeries.map((item) => [String(item?.sample || "").trim(), item]));
    const orderedSeries = orderedSamples.map((sample) => seriesMap.get(sample)).filter(Boolean);
    if (!orderedSeries.length) {
      if (stageNode) {
        stageNode.innerHTML = `
          <div class="empty-table-state">
            <strong>当前没有可展示的聚类组成图</strong>
            <p class="empty-copy">样本组成数据与距离矩阵未成功对齐。</p>
          </div>
        `;
      }
      if (legendNode) legendNode.innerHTML = "";
      return;
    }
    const rowCount = orderedSeries.length;
    const rowHeight = Math.max(12, Math.min(24, Math.floor(860 / Math.max(rowCount, 1))));
    const barWidth = Math.max(420, availableWidth - treeWidth - stripWidth - 10);
    const plotHeight = rowCount * rowHeight;
    const leafCenters = Object.fromEntries(orderedIndices.map((sampleIndex, position) => [sampleIndex, position * rowHeight + rowHeight / 2]));
    const leftTree = layoutCommunityCluster(tree, leafCenters, maxTreeHeight, "left");
    const leftTreeSvg = renderCommunityClusterBranches(leftTree, treeWidth, plotHeight, "left");
    const taxonTotals = new Map();
    orderedSeries.forEach((item) => {
      const segments = Array.isArray(item?.segments) ? item.segments : [];
      segments.forEach((segment) => {
        const label = String(segment?.label || "").trim();
        if (!isAnnotatedCommunityTaxon(label)) return;
        taxonTotals.set(label, (taxonTotals.get(label) || 0) + (Number(segment?.ratio) || 0));
      });
    });
    const rankedLabels = Array.from(taxonTotals.entries()).sort((a, b) => b[1] - a[1]).map(([label]) => label);
    const topLabels = rankedLabels.slice(0, Math.max(1, state.topN));
    const colorMap = new Map(topLabels.map((label, index) => [
      label,
      naturePalette[index % naturePalette.length],
    ]));
    const normalizedSeries = orderedSeries.map((item) => {
      const rawSegments = Array.isArray(item?.segments) ? item.segments : [];
      const bucket = new Map();
      rawSegments.forEach((segment) => {
        const label = String(segment?.label || "").trim();
        if (!isAnnotatedCommunityTaxon(label)) return;
        if (!topLabels.includes(label)) return;
        bucket.set(label, (bucket.get(label) || 0) + (Number(segment?.ratio) || 0));
      });
      const segments = topLabels
        .map((label) => ({
          label,
          ratio: Number(bucket.get(label) || 0),
        }))
        .filter((segment) => segment.ratio > 0.01);
      return {
        sample: String(item?.sample || ""),
        group: String(item?.group || groups[item?.sample] || "未分组"),
        segments,
      };
    });
    if (tagNode) tagNode.textContent = `Bray + ${state.level} + Top ${state.topN}`;
    if (stageNode) {
      stageNode.innerHTML = `
        <div class="community-beta-composition-plot" style="grid-template-columns:${treeWidth}px ${stripWidth}px ${barWidth}px; grid-template-rows:repeat(${rowCount}, ${rowHeight}px);">
          <div class="community-beta-composition-tree" style="grid-column:1; grid-row:1 / span ${rowCount}; width:${treeWidth}px; height:${plotHeight}px;">
            <svg viewBox="0 0 ${treeWidth} ${plotHeight}" preserveAspectRatio="none" role="img" aria-label="样本层次聚类树">
              ${leftTreeSvg}
            </svg>
          </div>
          ${normalizedSeries.map((item, rowIndex) => `
            <div
              class="community-beta-composition-group-strip"
              style="grid-column:2; grid-row:${rowIndex + 1}; background:${escapeHtml(groupColorMap.get(item.group) || "#b6b3ae")};"
              title="${escapeHtml(`${item.sample} | ${item.group}`)}"
              data-community-beta-sample="${escapeHtml(item.sample)}"
              data-community-beta-group="${escapeHtml(item.group)}"
            ></div>
            <div
              class="community-beta-composition-bar"
              style="grid-column:3; grid-row:${rowIndex + 1};"
              title="${escapeHtml(`${item.sample} | ${item.group}`)}"
              data-community-beta-sample="${escapeHtml(item.sample)}"
              data-community-beta-group="${escapeHtml(item.group)}"
            >
              ${item.segments.map((segment) => `
                <span
                  class="community-beta-composition-segment"
                  style="width:${Math.max(segment.ratio, 0)}%; --segment-color:${escapeHtml(colorMap.get(segment.label) || "#d8d2ca")};"
                  title="${escapeHtml(`${item.sample} | ${item.group} | ${segment.label} | ${formatRate(segment.ratio)}`)}"
                  data-community-beta-sample="${escapeHtml(item.sample)}"
                  data-community-beta-group="${escapeHtml(item.group)}"
                  data-community-beta-taxon="${escapeHtml(segment.label)}"
                  data-community-beta-ratio="${escapeHtml(formatRate(segment.ratio))}"
                ></span>
              `).join("")}
            </div>
          `).join("")}
        </div>
        <div class="chart-tooltip" data-community-beta-tooltip hidden></div>
      `;
      const tooltipNode = stageNode.querySelector("[data-community-beta-tooltip]");
      const positionTooltip = (event, html) => {
        if (!tooltipNode) return;
        tooltipNode.innerHTML = html;
        tooltipNode.hidden = false;
        tooltipNode.style.position = "fixed";
        tooltipNode.style.pointerEvents = "none";
        const tooltipWidth = 220;
        const viewportWidth = window.innerWidth || 1280;
        const viewportHeight = window.innerHeight || 800;
        const nextLeft = Math.min(event.clientX + 16, viewportWidth - tooltipWidth - 20);
        const nextTop = Math.min(Math.max(event.clientY - 18, 12), viewportHeight - 110);
        tooltipNode.style.left = `${Math.max(12, nextLeft)}px`;
        tooltipNode.style.top = `${Math.max(12, nextTop)}px`;
        tooltipNode.style.right = "auto";
      };
      const hideTooltip = () => {
        if (!tooltipNode) return;
        tooltipNode.hidden = true;
      };
      stageNode.querySelectorAll("[data-community-beta-sample]").forEach((node) => {
        node.addEventListener("mouseenter", (event) => {
          const sample = node.getAttribute("data-community-beta-sample") || "--";
          const group = node.getAttribute("data-community-beta-group") || "未分组";
          const taxon = node.getAttribute("data-community-beta-taxon");
          const ratio = node.getAttribute("data-community-beta-ratio");
          if (taxon) {
            positionTooltip(event, `<strong>${escapeHtml(taxon)}</strong><div>样本: ${escapeHtml(sample)}</div><div>分组: ${escapeHtml(group)}</div><div>${escapeHtml(state.level)}水平占比: ${escapeHtml(ratio || "--")}</div>`);
            return;
          }
          positionTooltip(event, `<strong>${escapeHtml(sample)}</strong><div>分组: ${escapeHtml(group)}</div>`);
        });
        node.addEventListener("mousemove", (event) => {
          const sample = node.getAttribute("data-community-beta-sample") || "--";
          const group = node.getAttribute("data-community-beta-group") || "未分组";
          const taxon = node.getAttribute("data-community-beta-taxon");
          const ratio = node.getAttribute("data-community-beta-ratio");
          if (taxon) {
            positionTooltip(event, `<strong>${escapeHtml(taxon)}</strong><div>样本: ${escapeHtml(sample)}</div><div>分组: ${escapeHtml(group)}</div><div>${escapeHtml(state.level)}水平占比: ${escapeHtml(ratio || "--")}</div>`);
            return;
          }
          positionTooltip(event, `<strong>${escapeHtml(sample)}</strong><div>分组: ${escapeHtml(group)}</div>`);
        });
        node.addEventListener("mouseleave", hideTooltip);
      });
    }
    if (legendNode) {
      legendNode.innerHTML = topLabels.map((label) => `
        <span class="community-beta-composition-legend-item" title="${escapeHtml(label)}">
          <i style="background:${escapeHtml(colorMap.get(label) || "#d8d2ca")}"></i>
          <span>${escapeHtml(label)}</span>
        </span>
      `).join("");
    }
  };

  if (rankSelect) {
    rankSelect.value = state.level;
    rankSelect.addEventListener("change", () => {
      state.level = rankSelect.value;
      renderPlot();
    });
  }
  if (topnSelect) {
    topnSelect.value = String(state.topN);
    topnSelect.addEventListener("change", () => {
      state.topN = Number(topnSelect.value || 8);
      renderPlot();
    });
  }
  renderPlot();
}

function renderCommunityQcAssetCards(containerId, taskId, items, section) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const rows = normalizeCommunityAssetRows(items);
  if (!rows.length) {
    renderCommunityAssetLinks(containerId, taskId, items, "当前没有可展示的质控输出", "后续生成 demux / DADA2 结果后，这里会显示对应文件入口。");
    return;
  }
  const demuxCard = rows.find((item) => item.label === "demux.qzv");
  const denoiseCard = rows.find((item) => item.label === "denoising-stats-dada2.qzv");
  const cards = [
    demuxCard ? {
      ...demuxCard,
      title: "demux.qzv",
      copy: "读取每样本 reads 数和质量统计，用于判断原始测序数据分布是否均衡。",
      stats: [
        ["样本数", String(section?.summary?.demux_sample_count ?? "--")],
        ["总 forward", String(section?.demux_preview?.summary?.forward_total ?? "--")],
        ["总 reverse", String(section?.demux_preview?.summary?.reverse_total ?? "--")],
      ],
    } : null,
    denoiseCard ? {
      ...denoiseCard,
      title: "denoising-stats-dada2.qzv",
      copy: "读取 DADA2 去噪统计，直接展示过滤保留率、merged 和非嵌合保留情况。",
      stats: [
        ["过滤保留率", section?.denoise_preview?.summary?.avg_pass_filter_pct != null ? `${section.denoise_preview.summary.avg_pass_filter_pct}%` : "--"],
        ["非嵌合保留率", section?.denoise_preview?.summary?.avg_non_chimeric_pct != null ? `${section.denoise_preview.summary.avg_non_chimeric_pct}%` : "--"],
        ["merged 总数", String(section?.denoise_preview?.summary?.merged_total ?? "--")],
      ],
    } : null,
  ].filter(Boolean);
  container.innerHTML = `
    <div class="community-asset-grid">
      ${cards.map((item) => `
        <a class="community-asset-link community-asset-link-rich" href="/api/tasks/${encodeURIComponent(taskId)}/report-asset/${encodeURIComponent(item.relativePath)}" target="_blank" rel="noopener noreferrer">
          <span class="community-asset-name">${escapeHtml(item.title)}</span>
          <span class="community-asset-copy">${escapeHtml(item.copy)}</span>
          <div class="community-asset-stat-grid">
            ${item.stats.map((stat) => `
              <span class="community-asset-stat">
                <small>${escapeHtml(stat[0])}</small>
                <strong>${escapeHtml(stat[1])}</strong>
              </span>
            `).join("")}
          </div>
          <strong class="community-asset-status">${escapeHtml(item.status === "ready" ? "已生成" : item.statusLabel)}</strong>
        </a>
      `).join("")}
    </div>
  `;
}

function renderCommunityRankAbundance(containerId, levelData) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const sampleSeries = Array.isArray(levelData?.sample_series) ? levelData.sample_series : [];
  if (!sampleSeries.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前分类水平暂无可展示结果</strong>
        <p class="empty-copy">请切换到其他分类水平，或确认 taxa barplot 结果是否已生成。</p>
      </div>
    `;
    return;
  }
  const topnSelect = document.getElementById("community-topn-select");
  const sortSelect = document.getElementById("community-sort-select");
  const viewSelect = document.getElementById("community-view-select");
  const requestedTopN = Math.max(1, Number(topnSelect?.value || 10));
  const sortMode = String(sortSelect?.value || "group");
  const viewMode = String(viewSelect?.value || "sample");
  const palette = ["#4e6177", "#8a6654", "#7a7158", "#5d7c83", "#6d6481", "#9b7a3f", "#4f7f6b", "#8a4d47", "#6f8096", "#557a95", "#8a7a4f", "#7c5d83", "#5a897e", "#9a6b62", "#697789"];
  const topRows = Array.isArray(levelData?.rows) ? levelData.rows : [];
  const topLabels = topRows.slice(0, requestedTopN).map((row) => String(row?.[0] || ""));
  const preparedSeries = sampleSeries.map((sample) => {
    const segments = Array.isArray(sample?.segments) ? sample.segments : [];
    const kept = [];
    let otherRatio = 0;
    segments.forEach((segment) => {
      const label = String(segment?.label || "");
      const ratio = Math.max(0, Number(segment?.ratio) || 0);
      if (label === "其他") {
        otherRatio += ratio;
        return;
      }
      if (topLabels.includes(label)) {
        kept.push({ label, ratio });
      } else {
        otherRatio += ratio;
      }
    });
    if (otherRatio > 0.01) {
      kept.push({ label: "其他", ratio: Number(otherRatio.toFixed(2)) });
    }
    return {
      sample: String(sample?.sample || "--"),
      group: String(sample?.group || "未分组"),
      segments: kept,
    };
  });
  const displaySeries = (viewMode === "group"
    ? (() => {
      const grouped = new Map();
      preparedSeries.forEach((sample) => {
        const groupName = String(sample.group || "未分组");
        if (!grouped.has(groupName)) {
          grouped.set(groupName, []);
        }
        grouped.get(groupName).push(sample);
      });
      return Array.from(grouped.entries()).map(([groupName, samples]) => {
        const segmentTotals = new Map();
        samples.forEach((sample) => {
          (sample.segments || []).forEach((segment) => {
            const label = String(segment.label || "");
            const ratio = Math.max(0, Number(segment.ratio) || 0);
            segmentTotals.set(label, (segmentTotals.get(label) || 0) + ratio);
          });
        });
        const averagedSegments = Array.from(segmentTotals.entries())
          .map(([label, total]) => ({
            label,
            ratio: Number((total / Math.max(samples.length, 1)).toFixed(2)),
          }))
          .filter((segment) => segment.ratio > 0.01)
          .sort((left, right) => right.ratio - left.ratio);
        return {
          sample: groupName,
          group: groupName,
          sampleCount: samples.length,
          segments: averagedSegments,
        };
      });
    })()
    : preparedSeries
  ).sort((left, right) => {
    if (sortMode === "sample") {
      return String(left.sample).localeCompare(String(right.sample), "zh-CN");
    }
    const groupDelta = String(left.group).localeCompare(String(right.group), "zh-CN");
    if (groupDelta !== 0) return groupDelta;
    return String(left.sample).localeCompare(String(right.sample), "zh-CN");
  });
  const colorMap = new Map();
  let colorIndex = 0;
  displaySeries.forEach((sample) => {
    (sample.segments || []).forEach((segment) => {
      if (!colorMap.has(segment.label)) {
        colorMap.set(segment.label, palette[colorIndex % palette.length]);
        colorIndex += 1;
      }
    });
  });
  const longestLabelLength = displaySeries.reduce((maxLength, sample) => Math.max(maxLength, String(sample?.sample || "").length), 0);
  const labelReserve = Math.min(170, Math.max(96, longestLabelLength * 9));
  const width = Math.max(900, displaySeries.length * 52 + 180);
  const height = 280 + labelReserve;
  const padLeft = 70;
  const padRight = 24;
  const padTop = 28;
  const padBottom = labelReserve;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const barWidth = Math.max(14, Math.min(30, innerWidth / Math.max(displaySeries.length, 1) - 6));
  const step = innerWidth / Math.max(displaySeries.length, 1);
  const yTicks = [0, 25, 50, 75, 100];
  const axisAndBars = displaySeries.map((sample, index) => {
    const x = padLeft + index * step + Math.max((step - barWidth) / 2, 2);
    let currentTop = padTop + innerHeight;
    const segments = (sample.segments || []).map((segment) => {
      const ratio = Math.max(0, Number(segment.ratio) || 0);
      const segmentHeight = innerHeight * (ratio / 100);
      const y = currentTop - segmentHeight;
      currentTop = y;
      return `
        <rect
          x="${x}"
          y="${y}"
          width="${barWidth}"
          height="${Math.max(segmentHeight, 0)}"
          rx="2"
          fill="${escapeHtml(colorMap.get(segment.label) || '#4e6177')}"
        >
          <title>${escapeHtml(`${sample.sample} | ${segment.label} | ${segment.ratio}%`)}</title>
        </rect>
      `;
    }).join("");
    return `
      <g class="community-stack-column">
        ${segments}
        <text class="chart-axis-label x-axis community-stack-label" x="${x + barWidth / 2}" y="${padTop + innerHeight + labelReserve - 22}" transform="rotate(-38 ${x + barWidth / 2} ${padTop + innerHeight + labelReserve - 22})">${escapeHtml(String(sample.sample || '--'))}</text>
      </g>
    `;
  }).join("");
  container.innerHTML = `
    <div class="mini-chart-card community-abundance-card">
      <span class="mini-chart-title">${escapeHtml(String(levelData?.level || '当前层级'))}水平样本组成堆积柱状图</span>
      <div class="community-stacked-chart-wrap">
        <svg class="sparkline-svg community-stacked-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMinYMin meet" role="img" aria-label="${escapeHtml(String(levelData?.level || '当前层级'))}样本组成堆积柱状图">
          ${yTicks.map((tick) => {
            const y = padTop + innerHeight - innerHeight * (tick / 100);
            return `
              <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
              <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${tick}%</text>
            `;
          }).join('')}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${padTop + innerHeight}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop + innerHeight}" x2="${width - padRight}" y2="${padTop + innerHeight}"></line>
          ${axisAndBars}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 8}">${viewMode === "group" ? "分组" : "样本"}</text>
          <text class="chart-axis-title chart-axis-title-y" x="22" y="${height / 2}">相对丰度 (%)</text>
        </svg>
      </div>
      <div class="community-chart-caption">当前显示前 ${requestedTopN} 个优势分类单元，其余合并为 Other；${viewMode === "group" ? "按分组展示组内样本平均相对丰度" : "按样本展示单样本组成"}；${viewMode === "group" ? `分组按${sortMode === "group" ? "分组名" : "分组名"}排序` : `样本按${sortMode === "group" ? "分组" : "样本名"}排序`}。</div>
      <div class="community-abundance-legend">
        ${Array.from(colorMap.entries()).map(([label, color]) => `
          <span class="community-abundance-legend-item">
            <i style="background:${escapeHtml(color)}"></i>
            <span>${escapeHtml(label)}</span>
          </span>
        `).join('')}
      </div>
    </div>
  `;
}

function renderCommunityAlphaWorkspace(alphaSection) {
  const chartContainer = document.getElementById("community-alpha-chart");
  const detailContainer = document.getElementById("community-alpha-detail");
  const pairwiseContainerId = "community-alpha-pairwise-table";
  const tabGroup = document.querySelector(".community-alpha-metric-tabs");
  if (!chartContainer || !detailContainer || !tabGroup) return;
  const preferredOrder = ["shannon", "observed_features", "pielou_evenness", "simpson", "faith_pd", "chao1", "goods_coverage"];
  const metricKeys = preferredOrder.filter((key) => alphaSection?.boxplots?.[key]).concat(
    Object.keys(alphaSection?.boxplots || {}).filter((key) => !preferredOrder.includes(key)),
  );
  const metrics = Object.fromEntries(metricKeys.map((key) => {
    const metric = alphaSection?.boxplots?.[key] || {};
    return [key, {
      metric,
      pairwise: alphaSection?.pairwise?.[key] || {},
      tabLabel: String(metric?.tab_label || metric?.y_label || key),
    }];
  }));
  tabGroup.innerHTML = metricKeys.map((key, index) => `
    <button class="report-tab-button ${index === 0 ? "active" : ""}" type="button" data-alpha-metric="${escapeHtml(key)}">${escapeHtml(metrics[key].tabLabel)}</button>
  `).join("");
  const buttons = Array.from(tabGroup.querySelectorAll("[data-alpha-metric]"));
  const firstMetric = buttons.find((button) => Array.isArray(metrics?.[button.dataset.alphaMetric]?.metric?.groups) && metrics[button.dataset.alphaMetric].metric.groups.length)?.dataset.alphaMetric || metricKeys[0] || "shannon";
  const renderMetric = (metricKey, activeGroupLabel = "") => {
    buttons.forEach((button) => button.classList.toggle("active", button.dataset.alphaMetric === metricKey));
    const config = metrics[metricKey] || { metric: {}, pairwise: {}, tabLabel: metricKey };
    renderCommunityAlphaBoxplotCard(config.metric, chartContainer.id, detailContainer.id, config.tabLabel, activeGroupLabel);
    buildTableCard(pairwiseContainerId, `${config.tabLabel} 两两比较 P 值`, config.pairwise?.columns || [], config.pairwise?.rows || []);
  };
  buttons.forEach((button) => {
    const metricKey = String(button.dataset.alphaMetric || "");
    button.disabled = !Array.isArray(metrics?.[metricKey]?.metric?.groups) || metrics[metricKey].metric.groups.length === 0;
    button.classList.toggle("hidden", button.disabled);
    button.addEventListener("click", () => renderMetric(metricKey));
  });
  renderMetric(firstMetric);
}

function renderCommunityQcRarefaction(containerId, alphaSection) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const rarefaction = alphaSection?.rarefaction || {};
  const curves = Array.isArray(rarefaction?.curves) ? rarefaction.curves : [];
  if (!curves.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有可展示的稀释曲线</strong>
        <p class="empty-copy">请先生成 alpha-rarefaction.qzv 后再查看质控稀释曲线。</p>
      </div>
    `;
    return;
  }
  container.classList.remove("empty-box");
  const width = 760;
  const height = 420;
  const padLeft = 72;
  const padRight = 24;
  const padTop = 26;
  const padBottom = 58;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const allPoints = curves.flatMap((curve) => Array.isArray(curve?.points) ? curve.points : []);
  const xValues = allPoints.map((point) => Number(point.x) || 0);
  const yValues = allPoints.map((point) => Number(point.y) || 0);
  const minX = Math.min(...xValues, 0);
  const maxX = Math.max(...xValues, 1);
  const minY = 0;
  const maxY = Math.max(...yValues, 1);
  const xScale = (value) => padLeft + ((value - minX) / Math.max(maxX - minX, 1e-6)) * innerWidth;
  const yScale = (value) => padTop + innerHeight - ((value - minY) / Math.max(maxY - minY, 1e-6)) * innerHeight;
  const palette = ["#5f7c61", "#c57b57", "#7e91b2", "#b56d72", "#8aa184", "#9e82b2", "#bea55c", "#6e95a2"];
  const colorMap = new Map(curves.map((curve, index) => [String(curve?.group || `group-${index}`), palette[index % palette.length]]));
  const xTicks = 4;
  const yTicks = 4;
  const gridX = Array.from({ length: xTicks + 1 }, (_, index) => {
    const ratio = index / xTicks;
    const value = minX + (maxX - minX) * ratio;
    const x = padLeft + innerWidth * ratio;
    return `
      <line class="chart-grid-line" x1="${x}" y1="${padTop}" x2="${x}" y2="${height - padBottom}"></line>
      <text class="chart-axis-label x-axis" x="${x}" y="${height - 20}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const gridY = Array.from({ length: yTicks + 1 }, (_, index) => {
    const ratio = index / yTicks;
    const value = minY + (maxY - minY) * ratio;
    const y = padTop + innerHeight - innerHeight * ratio;
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const suggestedDepth = Number(rarefaction?.suggested_depth) || 0;
  const suggestedLine = suggestedDepth > 0 ? `<line class="community-rarefaction-depth-line" x1="${xScale(suggestedDepth)}" y1="${padTop}" x2="${xScale(suggestedDepth)}" y2="${height - padBottom}"></line>` : "";
  const suggestedLabel = suggestedDepth > 0 ? `<text class="community-rarefaction-depth-label" x="${xScale(suggestedDepth) + 6}" y="${padTop + 14}">建议深度 ${escapeHtml(String(suggestedDepth))}</text>` : "";
  const paths = curves.map((curve) => {
    const pointsText = (Array.isArray(curve?.points) ? curve.points : []).map((point, index) => {
      const x = xScale(Number(point.x) || 0);
      const y = yScale(Number(point.y) || 0);
      return `${index === 0 ? "M" : "L"} ${x} ${y}`;
    }).join(" ");
    const group = String(curve?.group || "未分组");
    const n = Number(curve?.n || 0);
    return `
      <path
        class="community-rarefaction-line"
        d="${pointsText}"
        stroke="${escapeHtml(colorMap.get(group) || "#5f7c61")}"
        data-group="${escapeHtml(group)}"
      >
        <title>${escapeHtml(`${group} | n=${n}`)}</title>
      </path>
    `;
  }).join("");
  container.innerHTML = `
    <div class="mini-chart-card community-rarefaction-card">
      <span class="mini-chart-title">${escapeHtml(String(rarefaction?.label || "Observed Features 稀释曲线"))}</span>
      ${buildChartInsight(`按分组汇总展示 observed features 随测序深度的变化趋势；曲线越早趋于平缓，通常说明当前测序深度越接近饱和。`)}
      <div class="chart-canvas community-rarefaction-canvas" style="--chart-height:${height}px">
        <svg class="sparkline-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="群落分析 observed features 稀释曲线">
          ${gridX}
          ${gridY}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${height - padBottom}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${height - padBottom}" x2="${width - padRight}" y2="${height - padBottom}"></line>
          ${suggestedLine}
          ${paths}
          ${suggestedLabel}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 4}">${escapeHtml(String(rarefaction?.x_label || "测序深度"))}</text>
          <text class="chart-axis-title chart-axis-title-y" x="20" y="${height / 2}">${escapeHtml(String(rarefaction?.y_label || "观察到的特征数"))}</text>
        </svg>
      </div>
      <div class="community-pcoa-legend">
        ${curves.map((curve) => `
          <span class="community-pcoa-legend-item">
            <i style="background:${escapeHtml(colorMap.get(String(curve?.group || "未分组")) || "#5f7c61")}"></i>
            <span>${escapeHtml(`${String(curve?.group || "未分组")} (n=${Number(curve?.n || 0)})`)}</span>
          </span>
        `).join("")}
      </div>
    </div>
  `;
}

function renderCommunityAlphaBoxplotCard(metric, containerId, detailContainerId, tabLabel, activeGroupLabel = "") {
  const container = document.getElementById(containerId);
  const detailContainer = document.getElementById(detailContainerId);
  if (!container || !detailContainer) return;
  const groups = Array.isArray(metric?.groups) ? metric.groups : [];
  if (!groups.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有可展示的 ${escapeHtml(tabLabel || "Alpha")} 图形结果</strong>
        <p class="empty-copy">后续读到 alpha-rarefaction 结果后，这里会显示当前指标的组间分布。</p>
      </div>
    `;
    detailContainer.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有可展示的统计摘要</strong>
        <p class="empty-copy">请先选择有结果的 Alpha 指标。</p>
      </div>
    `;
    return;
  }
  const longestLabelLength = groups.reduce((maxLength, group) => Math.max(maxLength, String(group?.label || "").length), 0);
  const labelReserve = Math.min(132, Math.max(72, longestLabelLength * 6.4));
  const width = Math.max(920, groups.length * 84 + 180);
  const height = 270 + labelReserve;
  const padLeft = 74;
  const padRight = 28;
  const padTop = 26;
  const padBottom = labelReserve;
  const innerWidth = width - padLeft - padRight;
  const innerHeight = height - padTop - padBottom;
  const maxValue = Math.max(...groups.map((group) => Number(group?.max) || 0), 1);
  const minValue = Math.min(...groups.map((group) => Number(group?.min) || 0), 0);
  const range = Math.max(maxValue - minValue, 1);
  const yPosition = (value) => padTop + innerHeight - ((Number(value) - minValue) / range) * innerHeight;
  const step = innerWidth / Math.max(groups.length, 1);
  const boxWidth = Math.max(24, Math.min(46, step * 0.46));
  const ticks = 4;
  const yGrid = Array.from({ length: ticks + 1 }, (_, index) => {
    const ratio = index / ticks;
    const value = minValue + range * ratio;
    const y = padTop + innerHeight - innerHeight * ratio;
    return `
      <line class="chart-grid-line" x1="${padLeft}" y1="${y}" x2="${width - padRight}" y2="${y}"></line>
      <text class="chart-axis-label y-axis" x="${padLeft - 10}" y="${y + 4}">${escapeHtml(formatChartValue(value))}</text>
    `;
  }).join("");
  const currentGroup = groups.find((group) => String(group?.label || "") === activeGroupLabel) || groups[0];
  const boxes = groups.map((group, index) => {
    const xCenter = padLeft + index * step + step / 2;
    const minY = yPosition(group.min);
    const q1Y = yPosition(group.q1);
    const medianY = yPosition(group.median);
    const q3Y = yPosition(group.q3);
    const maxY = yPosition(group.max);
    const boxX = xCenter - boxWidth / 2;
    const boxHeight = Math.max(q1Y - q3Y, 2);
    const labelY = padTop + innerHeight + Math.max(34, labelReserve - 20);
    const isActive = String(group?.label || "") === String(currentGroup?.label || "");
    return `
      <g class="community-boxplot-group ${isActive ? "is-active" : ""}" data-alpha-group="${escapeHtml(String(group.label || ""))}">
        <line class="community-boxplot-whisker" x1="${xCenter}" y1="${maxY}" x2="${xCenter}" y2="${q3Y}"></line>
        <line class="community-boxplot-whisker" x1="${xCenter}" y1="${q1Y}" x2="${xCenter}" y2="${minY}"></line>
        <line class="community-boxplot-cap" x1="${xCenter - boxWidth * 0.32}" y1="${maxY}" x2="${xCenter + boxWidth * 0.32}" y2="${maxY}"></line>
        <line class="community-boxplot-cap" x1="${xCenter - boxWidth * 0.32}" y1="${minY}" x2="${xCenter + boxWidth * 0.32}" y2="${minY}"></line>
        <line class="community-boxplot-mean" x1="${xCenter - boxWidth * 0.24}" y1="${yPosition(group.mean)}" x2="${xCenter + boxWidth * 0.24}" y2="${yPosition(group.mean)}"></line>
        <rect class="community-boxplot-box" x="${boxX}" y="${q3Y}" width="${boxWidth}" height="${boxHeight}" rx="8"></rect>
        <line class="community-boxplot-median" x1="${boxX}" y1="${medianY}" x2="${boxX + boxWidth}" y2="${medianY}"></line>
        <text class="chart-axis-label x-axis community-stack-label" x="${xCenter}" y="${labelY}" transform="rotate(-38 ${xCenter} ${labelY})">${escapeHtml(String(group.label || "--"))}</text>
        <title>${escapeHtml(`${group.label} | n=${group.n} | median=${group.median} | Q1=${group.q1} | Q3=${group.q3}`)}</title>
      </g>
    `;
  }).join("");
  const pValue = metric?.p_value;
  const significance = metric?.significant;
  const pText = pValue == null ? "p = --" : `p = ${pValue}`;
  const significanceText = pValue == null ? "样本量不足，未计算显著性" : (significance ? "存在显著性差异" : "未见显著性差异");
  container.innerHTML = `
    <div class="mini-chart-card community-alpha-boxplot-card" data-alpha-metric-card="${escapeHtml(tabLabel || '')}">
      <div class="chart-titlebar community-alpha-titlebar">
        <span class="mini-chart-title">${escapeHtml(metric.label || "Alpha 箱线图")}</span>
        <span class="community-significance-badge ${significance ? "is-significant" : "is-neutral"}">${escapeHtml(significanceText)}</span>
      </div>
      <div class="chart-insight" role="note" aria-label="${escapeHtml(metric.label || "Alpha 箱线图说明")}">
        <span class="chart-insight-label">统计说明</span>
        <p>${escapeHtml(metric.test || "置换 Kruskal-Wallis")}；${escapeHtml(pText)}。箱体表示四分位区间，中线表示中位数，短横线表示均值，须表示最小值和最大值。</p>
      </div>
      <div class="community-stacked-chart-wrap community-boxplot-wrap">
        <svg class="community-stacked-svg community-boxplot-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMinYMin meet" role="img" aria-label="${escapeHtml(metric.label || "Alpha 箱线图")}">
          ${yGrid}
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${padTop + innerHeight}"></line>
          <line class="chart-axis-line" x1="${padLeft}" y1="${padTop + innerHeight}" x2="${width - padRight}" y2="${padTop + innerHeight}"></line>
          ${boxes}
          <text class="chart-axis-title" x="${width / 2}" y="${height - 10}">${escapeHtml(metric.x_label || "分组")}</text>
          <text class="chart-axis-title chart-axis-title-y" x="24" y="${height / 2}">${escapeHtml(metric.y_label || "指标值")}</text>
        </svg>
      </div>
    </div>
  `;
  const renderDetail = (groupLabel) => {
    const detailGroup = groups.find((item) => String(item?.label || "") === String(groupLabel || "")) || currentGroup;
    if (!detailGroup) return;
    detailContainer.innerHTML = `
      <div class="card-head">
        <h3>${escapeHtml(tabLabel || "Alpha 指标")}判读</h3>
        <span class="card-tag">${escapeHtml(String(detailGroup.label || "--"))}</span>
      </div>
      <div class="chart-insight" role="note" aria-label="${escapeHtml(String(detailGroup.label || "--"))} 统计说明">
        <span class="chart-insight-label">当前聚焦分组</span>
        <p>悬停或点击箱体可以切换分组详情。这里显示当前分组的样本量、中位数、四分位区间和均值，帮助判断该组分布位置与离散程度。</p>
      </div>
      <div class="mini-stat-grid community-alpha-detail-grid">
        <div class="mini-stat-card"><span class="mini-stat-label">样本数</span><strong>${escapeHtml(String(detailGroup.n ?? "--"))}</strong></div>
        <div class="mini-stat-card"><span class="mini-stat-label">中位数</span><strong>${escapeHtml(String(detailGroup.median ?? "--"))}</strong></div>
        <div class="mini-stat-card"><span class="mini-stat-label">Q1 / Q3</span><strong>${escapeHtml(`${detailGroup.q1 ?? "--"} / ${detailGroup.q3 ?? "--"}`)}</strong></div>
        <div class="mini-stat-card"><span class="mini-stat-label">均值</span><strong>${escapeHtml(String(detailGroup.mean ?? "--"))}</strong></div>
        <div class="mini-stat-card"><span class="mini-stat-label">最小 / 最大</span><strong>${escapeHtml(`${detailGroup.min ?? "--"} / ${detailGroup.max ?? "--"}`)}</strong></div>
        <div class="mini-stat-card"><span class="mini-stat-label">整体显著性</span><strong>${escapeHtml(significanceText)}</strong></div>
      </div>
      <p class="community-chart-caption">全局检验：${escapeHtml(metric.test || "置换 Kruskal-Wallis")}，${escapeHtml(pText)}。</p>
    `;
    container.querySelectorAll(".community-boxplot-group").forEach((node) => {
      node.classList.toggle("is-active", node.getAttribute("data-alpha-group") === String(detailGroup.label || ""));
    });
  };
  container.querySelectorAll(".community-boxplot-group").forEach((node) => {
    const groupLabel = node.getAttribute("data-alpha-group") || "";
    node.addEventListener("mouseenter", () => renderDetail(groupLabel));
    node.addEventListener("click", () => renderDetail(groupLabel));
  });
  renderDetail(String(currentGroup?.label || ""));
}

function bindCommunityRankTabs(section) {
  const tabGroup = document.querySelector(".community-rank-tabs");
  if (!tabGroup) return;
  const buttons = Array.from(tabGroup.querySelectorAll("[data-community-rank]"));
  const levels = section?.taxa_abundance?.levels || {};
  const availableOrder = ["门", "科", "属"];
  const firstAvailable = availableOrder.find((key) => Array.isArray(levels?.[key]?.rows) && levels[key].rows.length) || Object.keys(levels)[0] || "门";
  const renderRank = (rank) => {
    buttons.forEach((button) => button.classList.toggle("active", button.dataset.communityRank === rank));
    const levelData = levels?.[rank] || { level: rank, columns: [], rows: [] };
    renderCommunityRankAbundance("community-taxonomy-chart", levelData);
    buildTableCard("community-taxonomy-table", `${rank}水平相对丰度汇总`, levelData.columns || [], levelData.rows || []);
  };
  buttons.forEach((button) => {
    const rank = String(button.dataset.communityRank || "");
    const enabled = Array.isArray(levels?.[rank]?.rows) && levels[rank].rows.length > 0;
    button.disabled = !enabled;
    button.classList.toggle("hidden", !enabled);
    button.addEventListener("click", () => renderRank(rank));
  });
  document.getElementById("community-topn-select")?.addEventListener("change", () => {
    const active = buttons.find((button) => button.classList.contains("active"));
    renderRank(String(active?.dataset.communityRank || firstAvailable));
  });
  document.getElementById("community-sort-select")?.addEventListener("change", () => {
    const active = buttons.find((button) => button.classList.contains("active"));
    renderRank(String(active?.dataset.communityRank || firstAvailable));
  });
  document.getElementById("community-view-select")?.addEventListener("change", () => {
    const active = buttons.find((button) => button.classList.contains("active"));
    renderRank(String(active?.dataset.communityRank || firstAvailable));
  });
  renderRank(firstAvailable);
}

function buildPathoSourceNav() {
  return `
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-interpretation">1. 判读路径</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-overview">2. 概览</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-cluster">3. 成簇信息</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-distance">4. SNP 距离分布</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-ani">5. ANI 距离</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-grapetree">6. GrapeTree</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-mlst">7. MLST</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-mutation">8. 突变信息</a></div>
    <div class="report-nav-group"><a class="report-nav-link" href="#section-patho-core-tree">9. 遗传进化树</a></div>
  `;
}

function buildPathoSourceLayout() {
  return `
    <section class="report-intro-card patho-intro-card">
      <div class="report-intro-main">
        <p class="report-kicker">PathoSource Report</p>
        <h2 id="report-sample-title">分子溯源分析结果</h2>
        <p id="report-sample-copy">围绕成簇关系、SNP 距离、ANI、MLST 以及多树拓扑组织结果，方便快速判读样本间亲缘关系。</p>
        <div id="report-sample-switcher" class="report-sample-switcher hidden" aria-label="样本列表"></div>
      </div>
      <dl class="report-meta-grid">
        <div><dt>任务归属</dt><dd id="meta-owner">-</dd></div>
        <div><dt>用户组</dt><dd id="meta-group">-</dd></div>
        <div><dt>比对方法</dt><dd id="meta-asm-type">-</dd></div>
        <div><dt>建树方法</dt><dd id="meta-method">-</dd></div>
        <div><dt>输入路径</dt><dd id="meta-input">-</dd></div>
        <div><dt>输出目录</dt><dd id="meta-output">-</dd></div>
      </dl>
    </section>

    <section id="section-patho-interpretation" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 1</p>
        <h2>判读路径</h2>
        <p>先回答是否成簇，再看簇内距离是否支持传播关联，最后指出下一步最值得追的样本与证据。</p>
      </div>
      <div id="patho-interpretation-band" class="patho-interpretation-band"></div>
    </section>

    <section id="section-overview" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 2</p>
        <h2>溯源概览</h2>
        <p>先回答三个问题：纳入多少样本、聚成几簇、哪几棵树已经具备判读条件。</p>
      </div>
      <div id="overview-metrics" class="metric-grid"></div>
    </section>

    <section id="section-patho-cluster" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 3</p>
        <h2>成簇信息整理</h2>
        <p>读取 Cluster.tsv，先看聚类规模，再看每个簇内部 SNP 差异范围。</p>
      </div>
      <div id="patho-cluster-summary" class="mini-stat-grid"></div>
      <div id="patho-cluster-table" class="report-table-card"></div>
    </section>

    <section id="section-patho-distance" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 4</p>
        <h2>SNP 距离分布与矩阵</h2>
        <p>距离分布决定了整体群体结构，矩阵预览则方便你快速定位近距离样本对。</p>
      </div>
      <div class="two-column">
        <article class="result-card">
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">分布</span>
              <h3>样本距离分布</h3>
            </div>
            <span class="card-tag">dis_bin.tsv</span>
          </div>
          <div id="patho-distance-bins" class="empty-box">
            <p>展示 0、1-10、10-100、100-1000、1000+ 五个距离区间内的样本对数量。</p>
          </div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">摘要</span>
              <h3>SNP 距离矩阵摘要</h3>
            </div>
            <span class="card-tag">dis.mat.txt</span>
          </div>
          <div id="patho-snp-summary" class="mini-stat-grid"></div>
        </article>
      </div>
      <div id="patho-snp-matrix-table" class="report-table-card"></div>
    </section>

    <section id="section-patho-ani" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 5</p>
        <h2>ANI 距离</h2>
        <p>从全基因组相似性补充验证样本间亲缘关系，适合和 SNP 拓扑交叉判断。</p>
      </div>
      <div class="two-column">
        <article class="result-card">
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">摘要</span>
              <h3>ANI 概览</h3>
            </div>
            <span class="card-tag">Full_ANI.txt</span>
          </div>
          <div id="patho-ani-summary" class="mini-stat-grid"></div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">高相似样本</span>
              <h3>ANI 最高样本对</h3>
            </div>
            <span class="card-tag">Top Pairs</span>
          </div>
          <div id="patho-ani-top-pairs" class="report-table-card report-table-card-embedded"></div>
        </article>
      </div>
      <div id="patho-ani-preview-table" class="report-table-card"></div>
    </section>

    <section id="section-patho-grapetree" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 6</p>
        <h2>GrapeTree 拓扑</h2>
        <p>根据 grapetree.nwk 渲染主树视图，用于快速查看主簇、分支和孤立样本。</p>
      </div>
      <div id="patho-grapetree-card"></div>
    </section>

    <section id="section-patho-mlst" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 7</p>
        <h2>MLST 结果与树</h2>
        <p>保留每个样本的 MLST 结果，同时单独展示 mlst.nwk 对应的拓扑信息。</p>
      </div>
      <div class="two-column">
        <article class="result-card">
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">统计</span>
              <h3>MLST 汇总</h3>
            </div>
            <span class="card-tag">mlst.txt</span>
          </div>
          <div id="patho-mlst-summary" class="mini-stat-grid"></div>
        </article>
        <article class="result-card">
          <div class="card-head">
            <div class="card-title-stack">
              <span class="section-chip">Tree</span>
              <h3>MLST 进化树</h3>
            </div>
            <span class="card-tag">mlst.nwk</span>
          </div>
          <div id="patho-mlst-tree-card"></div>
        </article>
      </div>
      <div id="patho-mlst-table" class="report-table-card"></div>
    </section>

    <section id="section-patho-mutation" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 8</p>
        <h2>突变信息</h2>
        <p>读取 Mutate.tsv，汇总发生变异的位点、受影响样本以及样本突变负荷，方便定位关键差异位点。</p>
      </div>
      <article class="result-card">
        <div class="card-head">
          <div class="card-title-stack">
            <span class="section-chip">摘要</span>
            <h3>突变位点概览</h3>
          </div>
          <span class="card-tag">Mutate.tsv</span>
        </div>
        <div id="patho-mutation-summary" class="mini-stat-grid"></div>
      </article>
      <div id="patho-mutation-table" class="report-table-card"></div>
    </section>

    <section id="section-patho-core-tree" class="report-section">
      <div class="section-heading">
        <p class="report-kicker">Section 9</p>
        <h2>遗传进化树</h2>
        <p>读取 rmref.core.aln.contree，展示基于核心比对结果构建的系统发育关系。</p>
      </div>
      <div id="patho-core-tree-card"></div>
    </section>
  `;
}

function renderMiniStatGrid(containerId, items) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (!Array.isArray(items) || !items.length) {
    container.innerHTML = `<p class="empty-copy">暂无可展示的数据。</p>`;
    return;
  }
  container.innerHTML = items.map((item) => `
    <article class="mini-stat-card">
      <span>${escapeHtml(item.label)}</span>
      <strong>${escapeHtml(item.value)}</strong>
    </article>
  `).join("");
}

function encodeReportAssetPath(assetName) {
  return String(assetName || "")
    .split("/")
    .filter((part) => part !== "")
    .map((part) => encodeURIComponent(part))
    .join("/");
}

function buildReportAssetUrl(taskId, assetName, options = {}) {
  if (!taskId || !assetName) {
    throw new Error("缺少结果资产信息");
  }
  const shell = document.querySelector(".report-shell");
  const endpoint = new URL(`/api/tasks/${encodeURIComponent(taskId)}/report-asset/${encodeReportAssetPath(assetName)}`, window.location.origin);
  const currentSample = shell?.dataset?.selectedSample || new URLSearchParams(window.location.search).get("sample") || "";
  if (currentSample) {
    endpoint.searchParams.set("sample", currentSample);
  }
  if (options?.renderAs) {
    endpoint.searchParams.set("render_as", options.renderAs);
  }
  return endpoint.toString();
}

async function fetchReportAsset(taskId, assetName, responseType = "json", options = {}) {
  const response = await fetch(buildReportAssetUrl(taskId, assetName, options), { credentials: "same-origin" });
  if (!response.ok) {
    throw new Error(`读取 CGView 资产失败：${assetName}`);
  }
  if (responseType === "text") {
    return response.text();
  }
  return response.json();
}

let hivBootscanStyleMounted = false;

function ensureHivBootscanStyles() {
  if (hivBootscanStyleMounted) return;
  const style = document.createElement("style");
  style.textContent = `
    .hiv-bootscan-card {
      display: grid;
      gap: 18px;
    }
    .hiv-bootscan-toolbar {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
    }
    .hiv-bootscan-segmented,
    .hiv-bootscan-actions {
      display: inline-flex;
      align-items: center;
      flex-wrap: wrap;
      gap: 8px;
    }
    .hiv-bootscan-segmented button,
    .hiv-bootscan-actions button,
    .hiv-bootscan-legend button,
    .hiv-bootscan-downloads a {
      appearance: none;
      border: 1px solid rgba(74, 86, 105, 0.18);
      background: rgba(255, 255, 255, 0.9);
      color: #263247;
      border-radius: 999px;
      padding: 8px 12px;
      font: inherit;
      font-size: 0.84rem;
      line-height: 1;
      cursor: pointer;
      text-decoration: none;
      transition: background-color .16s ease, border-color .16s ease, color .16s ease, transform .16s ease;
    }
    .hiv-bootscan-segmented button.active,
    .hiv-bootscan-actions button.active,
    .hiv-bootscan-downloads a:hover,
    .hiv-bootscan-legend button.active {
      background: #243b67;
      border-color: #243b67;
      color: #f7f8fb;
    }
    .hiv-bootscan-segmented button:hover,
    .hiv-bootscan-actions button:hover,
    .hiv-bootscan-legend button:hover,
    .hiv-bootscan-downloads a:hover {
      transform: translateY(-1px);
    }
    .hiv-bootscan-workspace {
      display: grid;
      gap: 14px;
      padding: 18px;
      border: 1px solid rgba(42, 56, 84, 0.08);
      border-radius: 20px;
      background:
        linear-gradient(180deg, rgba(250, 248, 242, 0.95), rgba(255, 255, 255, 0.98)),
        rgba(255,255,255,0.96);
      box-shadow: inset 0 1px 0 rgba(255,255,255,0.85);
    }
    .hiv-bootscan-headline {
      display: flex;
      flex-wrap: wrap;
      align-items: baseline;
      justify-content: space-between;
      gap: 10px;
    }
    .hiv-bootscan-headline strong {
      font-size: 1rem;
      color: #1f2a3d;
    }
    .hiv-bootscan-headline span {
      font-size: 0.82rem;
      color: #5f6f86;
    }
    .hiv-bootscan-stage {
      position: relative;
      overflow: hidden;
      border-radius: 18px;
      border: 1px solid rgba(54, 68, 94, 0.12);
      background: #ffffff;
      padding: 14px 14px 10px;
    }
    .hiv-bootscan-svg {
      width: 100%;
      height: auto;
      display: block;
    }
    .hiv-bootscan-grid-line {
      stroke: rgba(101, 114, 136, 0.14);
      stroke-width: 1;
    }
    .hiv-bootscan-axis-line {
      stroke: rgba(42, 55, 77, 0.58);
      stroke-width: 1.35;
    }
    .hiv-bootscan-threshold {
      stroke: rgba(220, 74, 74, 0.72);
      stroke-width: 1.6;
      stroke-dasharray: 7 6;
    }
    .hiv-bootscan-path {
      fill: none;
      stroke-width: 2.35;
      stroke-linecap: round;
      stroke-linejoin: round;
      transition: opacity .14s ease;
    }
    .hiv-bootscan-path.is-muted {
      opacity: 0.1;
    }
    .hiv-bootscan-focus-line {
      stroke: rgba(26, 35, 56, 0.28);
      stroke-width: 1.4;
      stroke-dasharray: 4 5;
    }
    .hiv-bootscan-focus-dot {
      stroke: #fff;
      stroke-width: 1.6;
    }
    .hiv-bootscan-axis-text,
    .hiv-bootscan-caption {
      fill: #58677d;
      font-size: 12px;
    }
    .hiv-bootscan-axis-label {
      fill: #2a354d;
      font-size: 13px;
      font-weight: 600;
    }
    .hiv-bootscan-chart-title {
      fill: #24324a;
      font-size: 15px;
      font-weight: 700;
    }
    .hiv-bootscan-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    .hiv-bootscan-legend button {
      display: inline-flex;
      align-items: center;
      gap: 8px;
    }
    .hiv-bootscan-swatch {
      width: 10px;
      height: 10px;
      border-radius: 999px;
      flex: 0 0 auto;
    }
    .hiv-bootscan-meta {
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
    }
    .hiv-bootscan-meta-card {
      padding: 12px 14px;
      border-radius: 14px;
      background: rgba(248, 244, 236, 0.85);
      border: 1px solid rgba(57, 69, 90, 0.08);
      min-height: 74px;
    }
    .hiv-bootscan-meta-card span {
      display: block;
      font-size: 0.75rem;
      color: #69798f;
      margin-bottom: 6px;
    }
    .hiv-bootscan-meta-card strong {
      display: block;
      color: #223047;
      font-size: 0.98rem;
      line-height: 1.35;
      word-break: break-word;
    }
    .hiv-bootscan-tooltip {
      position: absolute;
      pointer-events: none;
      min-width: 190px;
      max-width: 240px;
      padding: 12px 14px;
      border-radius: 14px;
      background: rgba(24, 31, 46, 0.94);
      color: #f7f8fb;
      box-shadow: 0 14px 36px rgba(10, 16, 28, 0.24);
      font-size: 0.8rem;
      line-height: 1.45;
      opacity: 0;
      transform: translateY(6px);
      transition: opacity .12s ease, transform .12s ease;
      z-index: 4;
    }
    .hiv-bootscan-tooltip.visible {
      opacity: 1;
      transform: translateY(0);
    }
    .hiv-bootscan-tooltip strong {
      display: block;
      margin-bottom: 4px;
      font-size: 0.86rem;
    }
    .hiv-bootscan-empty {
      padding: 24px;
      border-radius: 18px;
      border: 1px dashed rgba(95, 109, 134, 0.28);
      color: #5e6d83;
      background: rgba(255,255,255,0.72);
    }
    .hiv-bootscan-downloads {
      display: inline-flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    @media (max-width: 900px) {
      .hiv-bootscan-meta {
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }
    }
    @media (max-width: 640px) {
      .hiv-bootscan-toolbar,
      .hiv-bootscan-headline {
        flex-direction: column;
        align-items: flex-start;
      }
      .hiv-bootscan-meta {
        grid-template-columns: minmax(0, 1fr);
      }
    }
  `;
  document.head.appendChild(style);
  hivBootscanStyleMounted = true;
}

function parseHivBootscanCsv(text) {
  const lines = String(text || "").trim().split(/\r?\n/).filter(Boolean);
  if (lines.length < 2) return { rows: [], groups: [], maxMidpoint: 0 };
  const headers = lines[0].split(",").map((item) => item.trim());
  const groups = headers.slice(3);
  const rows = lines.slice(1).map((line) => {
    const parts = line.split(",").map((item) => item.trim());
    const row = {
      midpoint: Number(parts[0] || 0),
      start: Number(parts[1] || 0),
      end: Number(parts[2] || 0),
      values: {},
    };
    groups.forEach((group, index) => {
      row.values[group] = Number(parts[index + 3] || 0);
    });
    return row;
  }).filter((row) => Number.isFinite(row.midpoint));
  return {
    rows,
    groups,
    maxMidpoint: rows.reduce((max, row) => Math.max(max, row.end || row.midpoint || 0), 0),
  };
}

function getHivBootscanPalette(groups = []) {
  const fixed = {
    B: "#3b93ff",
    A1: "#f05454",
    A2: "#ef8cc5",
    C: "#9d6a1f",
    D: "#d8b3df",
    F1: "#9fe000",
    F2: "#65c03f",
    G: "#55c878",
    H: "#ffbf1f",
    J: "#37cce3",
    K: "#8b5cf6",
  };
  const fallback = ["#3b93ff", "#f05454", "#ef8cc5", "#9d6a1f", "#9fe000", "#55c878", "#ffbf1f", "#37cce3", "#8b5cf6", "#7386a6"];
  return groups.reduce((acc, group, index) => {
    acc[group] = fixed[group] || fallback[index % fallback.length];
    return acc;
  }, {});
}

function getDefaultHivBootscanGroups(dataset) {
  const peaks = (dataset?.groups || []).map((group) => ({
    group,
    peak: Math.max(...(dataset.rows || []).map((row) => Number(row.values?.[group] || 0)), 0),
  }));
  return peaks
    .filter((item) => item.peak >= 1)
    .sort((a, b) => b.peak - a.peak)
    .slice(0, 5)
    .map((item) => item.group);
}

function formatBootscanSupport(value) {
  const numeric = Number(value || 0);
  if (!Number.isFinite(numeric)) return "--";
  if (numeric >= 99.995) return "100.0%";
  if (numeric >= 10) return `${numeric.toFixed(1)}%`;
  if (numeric >= 1) return `${numeric.toFixed(2)}%`;
  return `${numeric.toFixed(3)}%`;
}

async function initializeHivBootscanExplorer({
  containerId,
  taskId,
  overallCsvAsset,
  pureCsvAsset,
  embeddedData = {},
  summary = {},
}) {
  const container = document.getElementById(containerId);
  if (!(container instanceof HTMLElement) || !taskId || (!overallCsvAsset && !pureCsvAsset)) return;
  ensureHivBootscanStyles();
  container.innerHTML = `<div class="hiv-bootscan-empty">正在加载 bootscan 数据…</div>`;
  try {
    const overallEmbedded = String(embeddedData?.overall_csv_text || "").trim();
    const pureEmbedded = String(embeddedData?.pure_csv_text || "").trim();
    const [overallText, pureText] = await Promise.all([
      overallEmbedded ? Promise.resolve(overallEmbedded) : (overallCsvAsset ? fetchReportAsset(taskId, overallCsvAsset, "text") : Promise.resolve("")),
      pureEmbedded ? Promise.resolve(pureEmbedded) : (pureCsvAsset ? fetchReportAsset(taskId, pureCsvAsset, "text") : Promise.resolve("")),
    ]);
    const datasets = {
      overall: parseHivBootscanCsv(overallText),
      pure: parseHivBootscanCsv(pureText),
    };
    const availableModes = ["overall", "pure"].filter((mode) => Array.isArray(datasets[mode]?.rows) && datasets[mode].rows.length);
    if (!availableModes.length) {
      container.innerHTML = `<div class="hiv-bootscan-empty">当前样本没有可用于交互绘制的 bootscan 窗口数据。</div>`;
      return;
    }
    const defaultMode = availableModes.includes("overall") ? "overall" : availableModes[0];
    const palette = getHivBootscanPalette(Array.from(new Set([...datasets.overall.groups, ...datasets.pure.groups])));
    const state = {
      mode: defaultMode,
      enabled: new Set(getDefaultHivBootscanGroups(datasets[defaultMode])),
      focusIndex: -1,
    };
    if (!state.enabled.size) {
      state.enabled = new Set((datasets[defaultMode]?.groups || []).slice(0, 4));
    }
    container.innerHTML = `
      <div class="hiv-bootscan-card">
        <div class="hiv-bootscan-toolbar">
          <div class="hiv-bootscan-segmented" role="tablist" aria-label="Bootscan 视角切换">
            ${availableModes.map((mode) => `<button type="button" data-hiv-bootscan-mode="${mode}" class="${mode === state.mode ? "active" : ""}">${mode === "overall" ? "Overall" : "Pure subtype"}</button>`).join("")}
          </div>
          <div class="hiv-bootscan-actions">
            <button type="button" data-hiv-bootscan-action="focus">仅主信号</button>
            <button type="button" data-hiv-bootscan-action="all">全部亚型</button>
            <button type="button" data-hiv-bootscan-action="reset">重置筛选</button>
          </div>
        </div>
        <div class="hiv-bootscan-workspace">
          <div class="hiv-bootscan-headline">
            <strong data-hiv-bootscan-title>Bootscan Explorer</strong>
            <span data-hiv-bootscan-caption></span>
          </div>
          <div class="hiv-bootscan-stage">
            <svg class="hiv-bootscan-svg" viewBox="0 0 1120 440" preserveAspectRatio="xMidYMid meet" role="img" aria-label="HIV bootscan interactive plot"></svg>
            <div class="hiv-bootscan-tooltip" hidden></div>
          </div>
          <div class="hiv-bootscan-legend"></div>
          <div class="hiv-bootscan-meta">
            <article class="hiv-bootscan-meta-card"><span>Window</span><strong data-hiv-meta-window>--</strong></article>
            <article class="hiv-bootscan-meta-card"><span>Dominant lineage</span><strong data-hiv-meta-dominant>--</strong></article>
            <article class="hiv-bootscan-meta-card"><span>Peak support</span><strong data-hiv-meta-support>--</strong></article>
            <article class="hiv-bootscan-meta-card"><span>Tree context</span><strong data-hiv-meta-tree>--</strong></article>
          </div>
          <div class="hiv-bootscan-downloads">
            ${overallCsvAsset ? `<a href="${escapeHtml(buildReportAssetUrl(taskId, overallCsvAsset))}" target="_blank" rel="noreferrer">Overall CSV</a>` : ""}
            ${pureCsvAsset ? `<a href="${escapeHtml(buildReportAssetUrl(taskId, pureCsvAsset))}" target="_blank" rel="noreferrer">Pure CSV</a>` : ""}
          </div>
        </div>
      </div>
    `;
    const svg = container.querySelector(".hiv-bootscan-svg");
    const tooltip = container.querySelector(".hiv-bootscan-tooltip");
    const legend = container.querySelector(".hiv-bootscan-legend");
    const titleNode = container.querySelector("[data-hiv-bootscan-title]");
    const captionNode = container.querySelector("[data-hiv-bootscan-caption]");
    const windowNode = container.querySelector("[data-hiv-meta-window]");
    const dominantNode = container.querySelector("[data-hiv-meta-dominant]");
    const supportNode = container.querySelector("[data-hiv-meta-support]");
    const treeNode = container.querySelector("[data-hiv-meta-tree]");
    if (!(svg instanceof SVGElement) || !(tooltip instanceof HTMLElement) || !(legend instanceof HTMLElement)) return;

    const width = 1120;
    const height = 440;
    const margin = { top: 26, right: 34, bottom: 46, left: 58 };
    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;

    const treeSummaryText = () => {
      const key = state.mode === "overall" ? "overall_tree_support" : "pure_tree_support";
      const label = state.mode === "overall" ? "全参考树支持度" : "纯亚型树支持度";
      return `${label} ${String(summary?.[key] || "--")}`;
    };

    const linePathForGroup = (rows, group, xScale, yScale) => {
      return rows.map((row, index) => {
        const x = xScale(row.midpoint);
        const y = yScale(Number(row.values?.[group] || 0));
        return `${index === 0 ? "M" : "L"}${x.toFixed(2)} ${y.toFixed(2)}`;
      }).join(" ");
    };

    const renderLegend = (dataset) => {
      legend.innerHTML = dataset.groups.map((group) => `
        <button type="button" class="${state.enabled.has(group) ? "active" : ""}" data-hiv-group="${escapeHtml(group)}">
          <span class="hiv-bootscan-swatch" style="background:${palette[group]};"></span>
          <span>${escapeHtml(group)}</span>
        </button>
      `).join("");
      legend.querySelectorAll("[data-hiv-group]").forEach((button) => {
        button.addEventListener("click", () => {
          const group = String(button.getAttribute("data-hiv-group") || "").trim();
          if (!group) return;
          if (state.enabled.has(group) && state.enabled.size > 1) state.enabled.delete(group);
          else state.enabled.add(group);
          render();
        });
      });
    };

    const renderFocus = (dataset, index) => {
      if (index < 0 || index >= dataset.rows.length) {
        windowNode.textContent = "--";
        dominantNode.textContent = "--";
        supportNode.textContent = "--";
        treeNode.textContent = treeSummaryText();
        tooltip.hidden = true;
        tooltip.classList.remove("visible");
        return;
      }
      const row = dataset.rows[index];
      const ranked = dataset.groups
        .map((group) => ({ group, value: Number(row.values?.[group] || 0) }))
        .sort((a, b) => b.value - a.value);
      const dominant = ranked[0] || { group: "--", value: 0 };
      const runnerUp = ranked[1] || null;
      windowNode.textContent = `${row.start}-${row.end} nt`;
      dominantNode.textContent = runnerUp && runnerUp.value > 0.5
        ? `${dominant.group} > ${runnerUp.group}`
        : dominant.group;
      supportNode.textContent = formatBootscanSupport(dominant.value);
      treeNode.textContent = treeSummaryText();
      tooltip.innerHTML = `
        <strong>主导谱系：${escapeHtml(dominant.group)}</strong>
        中点 ${escapeHtml(String(row.midpoint))} nt<br>
        窗口 ${escapeHtml(String(row.start))}-${escapeHtml(String(row.end))} nt<br>
        支持度 ${escapeHtml(formatBootscanSupport(dominant.value))}
        ${runnerUp && runnerUp.value > 0.5 ? `<br>次优势谱系 ${escapeHtml(runnerUp.group)} · ${escapeHtml(formatBootscanSupport(runnerUp.value))}` : ""}
      `;
    };

    const render = () => {
      const dataset = datasets[state.mode];
      if (!dataset.rows.length) return;
      if (![...state.enabled].some((group) => dataset.groups.includes(group))) {
        state.enabled = new Set(getDefaultHivBootscanGroups(dataset));
      }
      const activeGroups = dataset.groups.filter((group) => state.enabled.has(group));
      const maxX = dataset.maxMidpoint || Math.max(...dataset.rows.map((row) => row.end || row.midpoint), 1);
      const xScale = (value) => margin.left + (Math.max(0, value) / maxX) * innerWidth;
      const yScale = (value) => margin.top + innerHeight - (Math.max(0, Math.min(100, value)) / 100) * innerHeight;
      const yTicks = [0, 20, 40, 60, 80, 100];
      const xTicks = Array.from({ length: 6 }, (_, idx) => Math.round((maxX / 5) * idx));
      titleNode.textContent = state.mode === "overall" ? "整体 bootscan" : "纯亚型 bootscan";
      captionNode.textContent = `窗口 400 bp · 步长 40 bp · 当前显示 ${activeGroups.length}/${dataset.groups.length} 条线`;
      renderLegend(dataset);

      const gridMarkup = yTicks.map((tick) => `
        <g>
          <line class="hiv-bootscan-grid-line" x1="${margin.left}" y1="${yScale(tick)}" x2="${width - margin.right}" y2="${yScale(tick)}"></line>
          <text class="hiv-bootscan-axis-text" x="${margin.left - 10}" y="${yScale(tick) + 4}" text-anchor="end">${tick}</text>
        </g>
      `).join("") + xTicks.map((tick) => `
        <g>
          <line class="hiv-bootscan-grid-line" x1="${xScale(tick)}" y1="${margin.top}" x2="${xScale(tick)}" y2="${margin.top + innerHeight}"></line>
          <text class="hiv-bootscan-axis-text" x="${xScale(tick)}" y="${height - 12}" text-anchor="middle">${tick}</text>
        </g>
      `).join("");

      const pathMarkup = dataset.groups.map((group) => `
        <path class="hiv-bootscan-path ${state.enabled.has(group) ? "" : "is-muted"}" d="${linePathForGroup(dataset.rows, group, xScale, yScale)}" stroke="${palette[group]}" data-group="${escapeHtml(group)}"></path>
      `).join("");

      svg.innerHTML = `
        <text class="hiv-bootscan-chart-title" x="${margin.left}" y="16">不同基因组窗口的 bootscan 支持度</text>
        <text class="hiv-bootscan-axis-label" x="${margin.left}" y="${margin.top - 6}">支持度</text>
        <text class="hiv-bootscan-axis-label" x="${margin.left + innerWidth / 2}" y="${height - 2}" text-anchor="middle">核苷酸位置</text>
        ${gridMarkup}
        <line class="hiv-bootscan-threshold" x1="${margin.left}" y1="${yScale(70)}" x2="${width - margin.right}" y2="${yScale(70)}"></line>
        <line class="hiv-bootscan-axis-line" x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top + innerHeight}"></line>
        <line class="hiv-bootscan-axis-line" x1="${margin.left}" y1="${margin.top + innerHeight}" x2="${width - margin.right}" y2="${margin.top + innerHeight}"></line>
        ${pathMarkup}
        <line class="hiv-bootscan-focus-line" x1="0" y1="${margin.top}" x2="0" y2="${margin.top + innerHeight}" visibility="hidden"></line>
        <g class="hiv-bootscan-focus-points"></g>
      `;

      const focusLine = svg.querySelector(".hiv-bootscan-focus-line");
      const focusPoints = svg.querySelector(".hiv-bootscan-focus-points");
      const updateFocus = (clientX) => {
        const rect = svg.getBoundingClientRect();
        const svgX = ((clientX - rect.left) / Math.max(rect.width, 1)) * width;
        const relativeX = Math.min(Math.max(svgX, margin.left), width - margin.right);
        let bestIndex = 0;
        let bestDistance = Infinity;
        dataset.rows.forEach((row, index) => {
          const distance = Math.abs(xScale(row.midpoint) - relativeX);
          if (distance < bestDistance) {
            bestDistance = distance;
            bestIndex = index;
          }
        });
        state.focusIndex = bestIndex;
        const row = dataset.rows[bestIndex];
        renderFocus(dataset, bestIndex);
        if (focusLine) {
          const focusX = xScale(row.midpoint);
          focusLine.setAttribute("x1", String(focusX));
          focusLine.setAttribute("x2", String(focusX));
          focusLine.setAttribute("visibility", "visible");
          if (focusPoints) {
            focusPoints.innerHTML = activeGroups.map((group) => `
              <circle class="hiv-bootscan-focus-dot" cx="${focusX}" cy="${yScale(Number(row.values?.[group] || 0))}" r="4.2" fill="${palette[group]}"></circle>
            `).join("");
          }
          const ranked = dataset.groups
            .map((group) => ({ group, value: Number(row.values?.[group] || 0) }))
            .sort((a, b) => b.value - a.value);
          tooltip.hidden = false;
          tooltip.classList.add("visible");
          const left = Math.min(Math.max(clientX - rect.left + 12, 8), rect.width - 220);
          const top = Math.max(8, ((yScale(ranked[0]?.value || 0) / height) * rect.height) - 22);
          tooltip.style.left = `${left}px`;
          tooltip.style.top = `${top}px`;
        }
      };

      renderFocus(dataset, state.focusIndex >= 0 ? state.focusIndex : 0);
      svg.onmousemove = (event) => updateFocus(event.clientX);
      svg.onmouseenter = (event) => updateFocus(event.clientX);
      svg.onmouseleave = () => {
        if (focusLine) focusLine.setAttribute("visibility", "hidden");
        if (focusPoints) focusPoints.innerHTML = "";
        renderFocus(dataset, -1);
      };
      container.querySelectorAll("[data-hiv-bootscan-mode]").forEach((button) => {
        button.classList.toggle("active", button.getAttribute("data-hiv-bootscan-mode") === state.mode);
        button.onclick = () => {
          const nextMode = String(button.getAttribute("data-hiv-bootscan-mode") || "");
          if (!nextMode || nextMode === state.mode) return;
          state.mode = nextMode;
          state.enabled = new Set(getDefaultHivBootscanGroups(datasets[state.mode]));
          state.focusIndex = -1;
          render();
        };
      });
      container.querySelectorAll("[data-hiv-bootscan-action]").forEach((button) => {
        button.onclick = () => {
          const action = String(button.getAttribute("data-hiv-bootscan-action") || "");
          const currentDataset = datasets[state.mode];
          if (action === "all") state.enabled = new Set(currentDataset.groups);
          else if (action === "focus" || action === "reset") state.enabled = new Set(getDefaultHivBootscanGroups(currentDataset));
          render();
        };
      });
    };
    render();
  } catch (error) {
    container.innerHTML = `<div class="hiv-bootscan-empty">Bootscan 数据加载失败：${escapeHtml(error instanceof Error ? error.message : "未知错误")}</div>`;
  }
}

let hivResistanceStyleMounted = false;
let hivMutationMapStyleMounted = false;

function ensureHivResistanceStyles() {
  if (hivResistanceStyleMounted) return;
  const style = document.createElement("style");
  style.textContent = `
    .hiv-resistance-shell {
      display: grid;
      gap: 18px;
    }
    .hiv-resistance-overview {
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 12px;
    }
    .hiv-resistance-overview-card,
    .hiv-resistance-gene-card,
    .hiv-resistance-panel,
    .hiv-resistance-score-card,
    .hiv-resistance-alert-card {
      border: 1px solid rgba(46, 60, 86, 0.08);
      background:
        linear-gradient(180deg, rgba(250, 248, 242, 0.92), rgba(255, 255, 255, 0.98)),
        rgba(255,255,255,0.96);
      border-radius: 18px;
      box-shadow: inset 0 1px 0 rgba(255,255,255,0.85);
    }
    .hiv-resistance-overview-card {
      padding: 14px 16px;
      min-height: 110px;
      display: grid;
      gap: 8px;
    }
    .hiv-resistance-overview-card span {
      font-size: 0.76rem;
      color: #63738a;
      letter-spacing: 0.01em;
      text-transform: uppercase;
    }
    .hiv-resistance-overview-card strong {
      font-size: 1rem;
      line-height: 1.35;
      color: #1f2a3d;
    }
    .hiv-resistance-overview-card p {
      margin: 0;
      font-size: 0.83rem;
      color: #667892;
      line-height: 1.45;
    }
    .hiv-resistance-tabbar {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    .hiv-resistance-tabbar button {
      appearance: none;
      border: 1px solid rgba(74, 86, 105, 0.18);
      background: rgba(255, 255, 255, 0.92);
      color: #24324a;
      border-radius: 999px;
      padding: 8px 12px;
      font: inherit;
      font-size: 0.84rem;
      line-height: 1;
      cursor: pointer;
      transition: background-color .16s ease, border-color .16s ease, color .16s ease, transform .16s ease;
    }
    .hiv-resistance-tabbar button:hover {
      transform: translateY(-1px);
    }
    .hiv-resistance-tabbar button.active {
      color: #fcfdff;
      border-color: #223552;
      background: #223552;
    }
    .hiv-resistance-panel {
      display: grid;
      gap: 18px;
      padding: 18px;
    }
    .hiv-resistance-panel-head {
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      gap: 12px;
      align-items: flex-start;
    }
    .hiv-resistance-panel-head h4,
    .hiv-resistance-gene-card h4 {
      margin: 0;
      font-size: 1.04rem;
      line-height: 1.25;
      color: #1d2940;
    }
    .hiv-resistance-panel-head p,
    .hiv-resistance-gene-card p,
    .hiv-resistance-alert-card p {
      margin: 4px 0 0;
      font-size: 0.84rem;
      line-height: 1.52;
      color: #66788e;
    }
    .hiv-resistance-panel-grid {
      display: grid;
      grid-template-columns: minmax(0, 1.5fr) minmax(320px, 0.95fr);
      gap: 16px;
      align-items: start;
    }
    .hiv-resistance-drug-list {
      display: grid;
      gap: 10px;
    }
    .hiv-resistance-drug-row {
      border-radius: 16px;
      border: 1px solid rgba(60, 74, 101, 0.08);
      background: rgba(255,255,255,0.94);
      padding: 14px 16px;
      display: grid;
      gap: 10px;
    }
    .hiv-resistance-drug-grid {
      display: grid;
      grid-template-columns: minmax(0, 1.25fr) 92px minmax(170px, 0.9fr);
      gap: 12px;
      align-items: center;
    }
    .hiv-resistance-drug-name strong {
      display: block;
      color: #1f2a3d;
      font-size: 0.98rem;
      line-height: 1.3;
    }
    .hiv-resistance-drug-name span,
    .hiv-resistance-drug-score span,
    .hiv-resistance-rule-summary span,
    .hiv-resistance-gene-meta span,
    .hiv-resistance-panel-meta span {
      display: block;
      color: #6a7b91;
      font-size: 0.76rem;
      margin-top: 4px;
    }
    .hiv-resistance-drug-score strong {
      display: block;
      color: #223047;
      font-size: 1rem;
      line-height: 1.2;
    }
    .hiv-resistance-level {
      display: inline-flex;
      align-items: center;
      gap: 8px;
      border-radius: 999px;
      padding: 7px 11px;
      font-size: 0.8rem;
      line-height: 1;
      font-weight: 600;
      white-space: nowrap;
    }
    .hiv-resistance-level::before {
      content: "";
      width: 8px;
      height: 8px;
      border-radius: 999px;
      background: currentColor;
      opacity: 0.8;
    }
    .hiv-resistance-level.is-susceptible {
      color: #22614a;
      background: rgba(220, 243, 233, 0.95);
      border: 1px solid rgba(34, 97, 74, 0.12);
    }
    .hiv-resistance-level.is-potential {
      color: #94611d;
      background: rgba(255, 242, 211, 0.96);
      border: 1px solid rgba(148, 97, 29, 0.12);
    }
    .hiv-resistance-level.is-low {
      color: #9d5b14;
      background: rgba(255, 232, 203, 0.96);
      border: 1px solid rgba(157, 91, 20, 0.12);
    }
    .hiv-resistance-level.is-intermediate {
      color: #b24d1f;
      background: rgba(255, 221, 207, 0.96);
      border: 1px solid rgba(178, 77, 31, 0.12);
    }
    .hiv-resistance-level.is-high {
      color: #b13030;
      background: rgba(255, 219, 219, 0.96);
      border: 1px solid rgba(177, 48, 48, 0.14);
    }
    .hiv-resistance-drug-details {
      display: grid;
      gap: 10px;
      padding-top: 2px;
    }
    .hiv-resistance-drug-details details {
      border-radius: 14px;
      border: 1px solid rgba(63, 78, 106, 0.08);
      background: rgba(247, 244, 236, 0.78);
      padding: 10px 12px;
    }
    .hiv-resistance-drug-details summary,
    .hiv-resistance-gene-card summary {
      cursor: pointer;
      color: #23324a;
      font-size: 0.84rem;
      font-weight: 600;
    }
    .hiv-resistance-detail-list,
    .hiv-resistance-comment-list {
      margin: 10px 0 0;
      padding-left: 18px;
      display: grid;
      gap: 8px;
      color: #304057;
      font-size: 0.83rem;
      line-height: 1.5;
    }
    .hiv-resistance-side {
      display: grid;
      gap: 12px;
    }
    .hiv-resistance-score-card {
      padding: 15px 16px;
      display: grid;
      gap: 12px;
    }
    .hiv-resistance-score-card h5,
    .hiv-resistance-alert-card h5 {
      margin: 0;
      color: #203047;
      font-size: 0.9rem;
    }
    .hiv-resistance-score-list {
      display: grid;
      gap: 10px;
    }
    .hiv-resistance-score-item {
      padding: 10px 12px;
      border-radius: 14px;
      background: rgba(255,255,255,0.92);
      border: 1px solid rgba(56, 70, 97, 0.08);
    }
    .hiv-resistance-score-item strong {
      display: block;
      color: #203047;
      font-size: 0.84rem;
      line-height: 1.45;
    }
    .hiv-resistance-score-item p {
      margin: 6px 0 0;
      color: #657791;
      font-size: 0.79rem;
      line-height: 1.5;
      word-break: break-word;
    }
    .hiv-resistance-gene-grid {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 14px;
    }
    .hiv-resistance-gene-card,
    .hiv-resistance-alert-card {
      padding: 16px;
      display: grid;
      gap: 12px;
      align-content: start;
    }
    .hiv-resistance-mutations {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }
    .hiv-resistance-mutation-pill {
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      padding: 7px 10px;
      background: rgba(36, 59, 103, 0.08);
      color: #233d67;
      border: 1px solid rgba(36, 59, 103, 0.12);
      font-size: 0.8rem;
      line-height: 1;
      white-space: nowrap;
    }
    .hiv-resistance-gene-meta {
      display: grid;
      gap: 8px;
    }
    .hiv-resistance-gene-meta strong {
      color: #203047;
      font-size: 0.92rem;
      line-height: 1.4;
    }
    .hiv-resistance-panel-meta {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      align-items: center;
    }
    .hiv-resistance-panel-meta strong {
      color: #203047;
      font-size: 0.84rem;
    }
    .hiv-resistance-empty {
      padding: 18px;
      border-radius: 18px;
      border: 1px dashed rgba(95, 109, 134, 0.28);
      color: #5e6d83;
      background: rgba(255,255,255,0.72);
    }
    @media (max-width: 1180px) {
      .hiv-resistance-overview {
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }
      .hiv-resistance-panel-grid,
      .hiv-resistance-gene-grid {
        grid-template-columns: minmax(0, 1fr);
      }
    }
    @media (max-width: 760px) {
      .hiv-resistance-overview {
        grid-template-columns: minmax(0, 1fr);
      }
      .hiv-resistance-drug-grid {
        grid-template-columns: minmax(0, 1fr);
      }
    }
  `;
  document.head.appendChild(style);
  hivResistanceStyleMounted = true;
}

function ensureHivMutationMapStyles() {
  if (hivMutationMapStyleMounted) return;
  const style = document.createElement("style");
  style.textContent = `
    .hiv-mutation-map-shell {
      display: grid;
      gap: 16px;
    }
    .hiv-mutation-map-browser {
      display: grid;
      gap: 16px;
    }
    .hiv-mutation-map-browser.is-collapsed {
      grid-template-columns: 220px minmax(0, 1fr);
      align-items: start;
    }
    .hiv-mutation-map-gene-list {
      display: none;
      border: 1px solid rgba(179, 41, 41, 0.2);
      background: rgba(255,255,255,0.96);
    }
    .hiv-mutation-map-browser.is-collapsed .hiv-mutation-map-gene-list {
      display: grid;
    }
    .hiv-mutation-map-gene-button {
      appearance: none;
      border: 0;
      border-top: 1px solid rgba(171, 145, 79, 0.18);
      background: rgba(255,255,255,0.95);
      color: #2d3138;
      text-align: left;
      padding: 18px 18px;
      font: inherit;
      font-size: 0.98rem;
      line-height: 1.3;
      cursor: pointer;
      transition: background-color .16s ease, color .16s ease, box-shadow .16s ease;
    }
    .hiv-mutation-map-gene-button:first-child {
      border-top: 0;
    }
    .hiv-mutation-map-gene-button:hover {
      background: rgba(248, 243, 232, 0.92);
    }
    .hiv-mutation-map-gene-button.active {
      background: rgba(255,255,255,1);
      box-shadow: inset 6px 0 0 #a12525;
      color: #1f2a3d;
      font-weight: 700;
    }
    .hiv-mutation-map-toolbar {
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
      gap: 12px;
      padding: 14px 16px;
      border: 1px solid rgba(201, 163, 48, 0.65);
      border-radius: 0;
      background: rgba(255, 252, 241, 0.74);
    }
    .hiv-mutation-map-toolbar strong {
      color: #2d3440;
      font-size: 0.98rem;
    }
    .hiv-mutation-map-actions {
      display: inline-flex;
      align-items: center;
      flex-wrap: wrap;
      gap: 12px;
    }
    .hiv-mutation-map-checkbox {
      display: inline-flex;
      align-items: center;
      gap: 8px;
      color: #32445f;
      font-size: 0.88rem;
      cursor: pointer;
    }
    .hiv-mutation-map-checkbox input {
      width: 16px;
      height: 16px;
      margin: 0;
      accent-color: #8e2f2f;
    }
    .hiv-mutation-map-download {
      appearance: none;
      border: 0;
      background: #982f2f;
      color: #fff8f8;
      border-radius: 4px;
      padding: 9px 14px;
      font: inherit;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
      transition: background-color .16s ease, transform .16s ease;
    }
    .hiv-mutation-map-download:hover {
      background: #812626;
      transform: translateY(-1px);
    }
    .hiv-mutation-map-stage {
      padding: 2px 10px 0 4px;
      overflow-x: auto;
    }
    .hiv-mutation-map-svg {
      width: 100%;
      min-width: 1280px;
      height: auto;
      display: block;
      overflow: visible;
    }
    .hiv-mutation-gene-title {
      fill: #2d2f34;
      font-size: 18px;
      font-weight: 700;
    }
    .hiv-mutation-axis {
      stroke: #1f232a;
      stroke-width: 2;
    }
    .hiv-mutation-axis-tick {
      stroke: #1f232a;
      stroke-width: 1.5;
    }
    .hiv-mutation-axis-label {
      fill: #2f3338;
      font-size: 11px;
      font-weight: 500;
    }
    .hiv-mutation-track {
      stroke: rgba(255,255,255,0.42);
      stroke-width: 1;
    }
    .hiv-mutation-track-label {
      font-size: 14px;
      font-weight: 600;
    }
    .hiv-mutation-stick {
      stroke: #30343b;
      stroke-width: 1.4;
      fill: none;
    }
    .hiv-mutation-stick.is-highlight {
      stroke: #1d8de3;
      stroke-width: 2.1;
    }
    .hiv-mutation-stick.is-alert {
      stroke: #c2362f;
      stroke-width: 1.9;
    }
    .hiv-mutation-cap {
      stroke: inherit;
      stroke-width: inherit;
    }
    .hiv-mutation-label {
      fill: #2d3138;
      font-size: 11px;
      font-weight: 500;
    }
    .hiv-mutation-label.is-highlight {
      fill: #1d8de3;
      font-weight: 700;
    }
    .hiv-mutation-label.is-alert {
      fill: #c2362f;
      font-weight: 700;
    }
    .hiv-mutation-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      color: #5a6980;
      font-size: 0.8rem;
      padding: 0 4px;
    }
    .hiv-mutation-legend span {
      display: inline-flex;
      align-items: center;
      gap: 7px;
    }
    .hiv-mutation-legend i {
      width: 12px;
      height: 2px;
      display: inline-block;
      background: #30343b;
    }
    .hiv-mutation-legend i.highlight {
      background: #1d8de3;
      height: 3px;
    }
    .hiv-mutation-legend i.alert {
      background: #c2362f;
      height: 3px;
    }
    .hiv-mutation-map-empty {
      padding: 18px;
      border-radius: 18px;
      border: 1px dashed rgba(95, 109, 134, 0.28);
      color: #5e6d83;
      background: rgba(255,255,255,0.72);
    }
    @media (max-width: 760px) {
      .hiv-mutation-map-toolbar {
        align-items: flex-start;
        flex-direction: column;
      }
      .hiv-mutation-map-browser.is-collapsed {
        grid-template-columns: minmax(0, 1fr);
      }
    }
  `;
  document.head.appendChild(style);
  hivMutationMapStyleMounted = true;
}

function normalizeHivResistanceLevelName(value) {
  const text = String(value || "").trim().toLowerCase();
  if (!text) return "susceptible";
  if (text.includes("high")) return "high";
  if (text.includes("intermediate")) return "intermediate";
  if (text.includes("low-level") && text.includes("potential")) return "potential";
  if (text.includes("low-level") || text.includes("low level")) return "low";
  return "susceptible";
}

function getHivResistanceLevelRank(value) {
  const normalized = normalizeHivResistanceLevelName(value);
  if (normalized === "high") return 5;
  if (normalized === "intermediate") return 4;
  if (normalized === "low") return 3;
  if (normalized === "potential") return 2;
  return 1;
}

function localizeHivAssignmentLabel(value) {
  const text = String(value || "").trim();
  const upper = text.toUpperCase();
  if (!text) return "--";
  if (upper === "PURE") return "纯亚型";
  if (upper === "PURE-LIKE") return "近似纯亚型";
  if (upper === "CRF") return "循环重组型（CRF）";
  if (upper === "CRF-LIKE") return "近似循环重组型";
  if (upper === "PURE RECOMBINANT") return "纯亚型重组";
  if (upper === "CRF PURE RECOMBINANT") return "CRF/纯亚型重组";
  if (upper === "POTENTIAL RECOMBINANT") return "疑似重组";
  if (upper === "NOT ASSIGNED") return "未定型";
  return text;
}

function parseHivResistanceList(value, splitOnComma = true) {
  if (Array.isArray(value)) {
    return value.map((item) => String(item || "").trim()).filter(Boolean);
  }
  const text = String(value || "").trim();
  if (!text) return [];
  if ((text.startsWith("[") && text.endsWith("]")) || (text.startsWith("{") && text.endsWith("}"))) {
    try {
      const parsed = JSON.parse(text);
      return parseHivResistanceList(parsed, splitOnComma);
    } catch (error) {
      // ignore parse failures and fall back to plain text splitting
    }
  }
  const pattern = splitOnComma ? /(?:\r?\n|；|;|\s+\|\s+|,\s*)+/ : /(?:\r?\n|；|;|\s+\|\s+)+/;
  return text.split(pattern).map((item) => item.trim()).filter(Boolean);
}

function buildHivResistanceData(table, payload, mutationPanels) {
  const columns = Array.isArray(table?.columns) ? table.columns : [];
  const rows = Array.isArray(table?.rows) ? table.rows : [];
  const payloadRoot = payload && typeof payload === "object" ? payload : {};
  const payloadSample = payloadRoot?.sample && typeof payloadRoot.sample === "object" ? payloadRoot.sample : {};
  let drugResults = Array.isArray(payloadSample?.drug_results) ? payloadSample.drug_results : [];
  if (!drugResults.length && columns.length && rows.length) {
    drugResults = rows.map((row) => {
      const valueOf = (name) => {
        const index = columns.indexOf(name);
        return index >= 0 && Array.isArray(row) ? row[index] : "";
      };
      return {
        drug_class: valueOf("drug_class"),
        drug: valueOf("drug"),
        fullname: valueOf("fullname"),
        score: valueOf("score"),
        level: valueOf("level"),
        level_name: valueOf("level_name"),
        sir: valueOf("sir"),
        triggered_rules: parseHivResistanceList(valueOf("triggered_rules"), false),
        result_comments: parseHivResistanceList(valueOf("result_comments"), false),
      };
    });
  }
  const grouped = new Map();
  drugResults.forEach((entry) => {
    if (!entry || typeof entry !== "object") return;
    const className = String(entry.drug_class || "").trim() || "Other";
    if (!grouped.has(className)) grouped.set(className, []);
    grouped.get(className).push({
      drug_class: className,
      drug: String(entry.drug || "").trim() || "--",
      fullname: String(entry.fullname || "").trim() || "",
      score: String(entry.score ?? "").trim(),
      level_name: String(entry.level_name || "Susceptible").trim() || "Susceptible",
      sir: String(entry.sir || "").trim(),
      triggered_rules: parseHivResistanceList(entry.triggered_rules, false),
      result_comments: parseHivResistanceList(entry.result_comments, false),
    });
  });
  const panelMutations = {};
  (Array.isArray(mutationPanels) ? mutationPanels : []).forEach((item) => {
    const label = String(item?.label || "").trim();
    const value = String(item?.value || "").trim();
    if (!label || !value || value === "-") return;
    const gene = label.split(/\s+/)[0].toUpperCase();
    panelMutations[gene] = parseHivResistanceList(value);
  });
  const inputMutations = payloadSample?.input_mutations && typeof payloadSample.input_mutations === "object"
    ? payloadSample.input_mutations
    : {};
  const comments = Array.isArray(payloadSample?.mutation_comments) ? payloadSample.mutation_comments : [];
  const sequenceAlerts = Array.isArray(payloadSample?.sequence_alerts) ? payloadSample.sequence_alerts : [];
  const classOrder = ["NRTI", "NNRTI", "PI", "INSTI", "CAI"];
  const classSummaries = classOrder
    .filter((className) => grouped.has(className))
    .map((className) => {
      const items = grouped.get(className) || [];
      const top = items.reduce((best, item) => {
        return getHivResistanceLevelRank(item.level_name) > getHivResistanceLevelRank(best?.level_name)
          ? item
          : best;
      }, items[0] || null);
      const resistantCount = items.filter((item) => getHivResistanceLevelRank(item.level_name) >= 3).length;
      return {
        className,
        items,
        topLevel: top?.level_name || "Susceptible",
        topDrug: top?.drug || "--",
        resistantCount,
      };
    });
  const genes = ["PR", "RT", "IN", "CA"]
    .map((gene) => {
      const observed = Array.isArray(inputMutations?.[gene]) ? inputMutations[gene] : (panelMutations[gene] || []);
      const geneComments = comments.filter((entry) => String(entry?.gene || "").trim().toUpperCase() === gene);
      return {
        gene,
        mutations: observed.map((item) => String(item || "").trim()).filter(Boolean),
        comments: geneComments,
      };
    })
    .filter((entry) => entry.mutations.length || entry.comments.length);
  return {
    algorithm: payloadRoot?.algorithm && typeof payloadRoot.algorithm === "object" ? payloadRoot.algorithm : {},
    classOrder,
    classSummaries,
    grouped,
    genes,
    sequenceAlerts: sequenceAlerts.map((item) => String(item || "").trim()).filter(Boolean),
  };
}

const HIV_MUTATION_TRACKS = [
  {
    gene: "PR",
    title: "蛋白酶（PR）",
    shortLabel: "蛋白酶",
    length: 99,
    color: "#a9d97f",
    labelColor: "#3d9830",
    ticks: [1, 5, 10, 15, 25, 30, 35, 40, 45, 50, 55, 65, 70, 75, 80, 85, 90, 99],
  },
  {
    gene: "RT",
    title: "逆转录酶（RT）",
    shortLabel: "逆转录酶",
    length: 560,
    color: "#9cc5df",
    labelColor: "#2474b0",
    ticks: [1, 20, 40, 65, 85, 110, 130, 155, 175, 200, 220, 245, 270, 290, 315, 335, 390, 560],
    breakAfter: 335,
  },
  {
    gene: "IN",
    title: "整合酶（IN）",
    shortLabel: "整合酶",
    length: 288,
    color: "#cab5df",
    labelColor: "#7f57b8",
    ticks: [1, 15, 30, 50, 65, 80, 100, 115, 130, 150, 165, 180, 200, 215, 230, 250, 265, 288],
  },
];

function parseHivMutationToken(token) {
  const text = String(token || "").trim();
  if (!text) return null;
  const match = text.match(/^([A-Za-z*]+)?(\d+)([A-Za-z*]+)$/);
  if (!match) return null;
  const ref = String(match[1] || "");
  const position = Number(match[2] || 0);
  const alt = String(match[3] || "");
  if (!Number.isFinite(position) || position <= 0) return null;
  return {
    raw: text,
    ref,
    alt,
    position,
    isAlert: alt.includes("*") || alt.includes("d") || alt.includes("i"),
  };
}

function buildHivMutationTracks(data) {
  const commentIndex = new Map();
  (data?.genes || []).forEach((geneEntry) => {
    const gene = String(geneEntry?.gene || "").trim().toUpperCase();
    (geneEntry?.comments || []).forEach((comment) => {
      const key = `${gene}:${String(comment?.condition || "").trim()}`;
      if (key) commentIndex.set(key, true);
    });
  });
  return HIV_MUTATION_TRACKS.map((track) => {
    const geneEntry = (data?.genes || []).find((item) => String(item?.gene || "").trim().toUpperCase() === track.gene);
    const mutations = Array.isArray(geneEntry?.mutations) ? geneEntry.mutations : [];
    const parsed = mutations
      .map((item) => parseHivMutationToken(item))
      .filter(Boolean)
      .map((item) => ({
        ...item,
        isHighlight: commentIndex.has(`${track.gene}:${item.raw}`),
      }))
      .sort((a, b) => a.position - b.position || a.raw.localeCompare(b.raw));
    return {
      ...track,
      mutations: parsed,
    };
  }).filter((track) => track.mutations.length);
}

function escapeXml(value) {
  return escapeHtml(value);
}

function downloadSvgMarkup(filename, markup) {
  const blob = new Blob([markup], { type: "image/svg+xml;charset=utf-8" });
  const objectUrl = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = objectUrl;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(objectUrl);
}

function renderHivMutationMap({
  containerId,
  table,
  payload,
  mutationPanels,
}) {
  const container = document.getElementById(containerId);
  if (!(container instanceof HTMLElement)) return;
  ensureHivMutationMapStyles();
  const data = buildHivResistanceData(table, payload, mutationPanels);
  const tracks = buildHivMutationTracks(data);
  if (!tracks.length) {
    container.innerHTML = `<div class="hiv-mutation-map-empty">当前没有可绘制的 PR / RT / IN 突变位点。</div>`;
    return;
  }

  const state = { collapsed: false, activeGene: tracks[0]?.gene || "" };
  container.innerHTML = `
    <div class="hiv-mutation-map-shell">
      <div class="hiv-mutation-map-toolbar">
        <strong>序列质量评估</strong>
        <div class="hiv-mutation-map-actions">
          <label class="hiv-mutation-map-checkbox">
            <input type="checkbox" data-hiv-collapse-genes>
            <span>折叠基因</span>
          </label>
          <button type="button" class="hiv-mutation-map-download" data-hiv-download-svg>下载 SVG</button>
        </div>
      </div>
      <div class="hiv-mutation-map-browser">
        <div class="hiv-mutation-map-gene-list"></div>
        <div class="hiv-mutation-map-stage">
          <div class="hiv-mutation-map-svg-wrap"></div>
        </div>
      </div>
      <div class="hiv-mutation-legend">
        <span><i></i>观测突变</span>
        <span><i class="highlight"></i>HIVDB 注释位点</span>
        <span><i class="alert"></i>终止/插缺相关位点</span>
      </div>
      <details class="report-detail-block">
        <summary>查看原始突变清单</summary>
        <div class="mini-stat-grid" style="margin-top:12px;">
          ${tracks.map((track) => `
            <article class="mini-stat-card" style="align-items:flex-start;">
              <span>${escapeHtml(track.title)}</span>
              <strong style="font-size:0.95rem;line-height:1.55;white-space:normal;word-break:break-word;">${escapeHtml(track.mutations.map((item) => item.raw).join(", "))}</strong>
            </article>
          `).join("")}
          ${Array.isArray(mutationPanels) ? mutationPanels.filter((item) => String(item?.label || "").includes("候选父本")).map((item) => `
            <article class="mini-stat-card" style="align-items:flex-start;">
              <span>${escapeHtml(String(item?.label || "--"))}</span>
              <strong style="font-size:0.95rem;line-height:1.55;white-space:normal;word-break:break-word;">${escapeHtml(String(item?.value || "--"))}</strong>
            </article>
          `).join("") : ""}
        </div>
      </details>
    </div>
  `;

  const svgWrap = container.querySelector(".hiv-mutation-map-svg-wrap");
  const browserNode = container.querySelector(".hiv-mutation-map-browser");
  const geneListNode = container.querySelector(".hiv-mutation-map-gene-list");
  const collapseInput = container.querySelector("[data-hiv-collapse-genes]");
  const downloadButton = container.querySelector("[data-hiv-download-svg]");
  if (!(svgWrap instanceof HTMLElement)) return;

  const render = () => {
    const visibleTracks = state.collapsed
      ? tracks.filter((track) => track.gene === state.activeGene)
      : tracks;
    if (!visibleTracks.length) return;
    if (browserNode instanceof HTMLElement) {
      browserNode.classList.toggle("is-collapsed", state.collapsed);
    }
    if (geneListNode instanceof HTMLElement) {
      geneListNode.innerHTML = tracks.map((track) => `
        <button type="button" class="hiv-mutation-map-gene-button ${state.collapsed && track.gene === state.activeGene ? "active" : ""}" data-hiv-gene="${escapeHtml(track.gene)}">
          ${escapeHtml(track.title)}
        </button>
      `).join("");
      geneListNode.querySelectorAll("[data-hiv-gene]").forEach((button) => {
        button.addEventListener("click", () => {
          const gene = String(button.getAttribute("data-hiv-gene") || "").trim();
          if (!gene) return;
          state.activeGene = gene;
          state.collapsed = true;
          if (collapseInput instanceof HTMLInputElement) collapseInput.checked = true;
          render();
        });
      });
    }
    const width = 1880;
    const left = 138;
    const right = 220;
    const axisStart = left;
    const axisEnd = width - right;
    const usableWidth = axisEnd - axisStart;
    const baseGap = state.collapsed ? 170 : 248;
    const topPad = 40;
    const totalHeight = topPad + visibleTracks.length * baseGap + 24;
    const axisYStart = 84;
    const trackOffset = 46;
    const labelBaseOffset = 34;
    const fontFamily = "Source Sans Pro, system-ui, sans-serif";
    const rows = [];
    visibleTracks.forEach((track, index) => {
      const top = topPad + index * baseGap;
      const axisY = top + axisYStart;
      const barY = axisY + trackOffset;
      const barHeight = 26;
      const labelStartY = barY + labelBaseOffset;
      const posToX = (position) => {
        const ratio = (Math.max(1, Math.min(track.length, position)) - 1) / Math.max(track.length - 1, 1);
        return axisStart + ratio * usableWidth;
      };
      const tickMarkup = (track.ticks || []).map((tick) => {
        const x = posToX(tick);
        return `
          <line class="hiv-mutation-axis-tick" x1="${x}" y1="${axisY}" x2="${x}" y2="${axisY + 10}"></line>
          <text class="hiv-mutation-axis-label" x="${x}" y="${axisY - 6}" text-anchor="middle">${escapeXml(String(tick))}</text>
        `;
      }).join("");
      const breakMarkup = track.breakAfter ? `
        <path d="M ${posToX(track.breakAfter) + 28} ${axisY + 1} l 6 -7 l 0 14 l 6 -7" fill="none" stroke="#1f232a" stroke-width="2"></path>
      ` : "";
      const mutationMarkup = [];
      let denseRightIndex = 0;
      let prevX = -Infinity;
      track.mutations.forEach((mutation, mutationIndex) => {
        const x = posToX(mutation.position);
        const nextX = mutationIndex + 1 < track.mutations.length ? posToX(track.mutations[mutationIndex + 1].position) : x + 999;
        const nearRight = x > axisEnd - 220;
        const tightCluster = x - prevX < 18 || nextX - x < 18;
        const lane = tightCluster ? (mutationIndex % 3) : 0;
        let labelX = x;
        let labelY = labelStartY + lane * 18;
        let pathD = `M ${x} ${barY} L ${x} ${labelY - 10}`;
        if (nearRight) {
          denseRightIndex += 1;
          labelX = Math.max(axisStart + 160, axisEnd - 42 - denseRightIndex * 22);
          labelY = labelStartY + denseRightIndex * 13;
          pathD = `M ${x} ${barY} L ${x} ${labelY - 28} L ${labelX + 12} ${labelY - 28} L ${labelX + 12} ${labelY - 10}`;
        } else if (tightCluster) {
          labelX = x + lane * 8;
          pathD = `M ${x} ${barY} L ${x} ${labelY - 14} L ${labelX} ${labelY - 10}`;
        }
        const lineClass = mutation.isAlert ? "is-alert" : (mutation.isHighlight ? "is-highlight" : "");
        const textClass = mutation.isAlert ? "is-alert" : (mutation.isHighlight ? "is-highlight" : "");
        mutationMarkup.push(`
          <path class="hiv-mutation-stick ${lineClass}" d="${pathD}"></path>
          <line class="hiv-mutation-cap ${lineClass}" x1="${x - 5}" y1="${barY}" x2="${x + 5}" y2="${barY}"></line>
          <text class="hiv-mutation-label ${textClass}" x="${labelX}" y="${labelY}" transform="rotate(-63 ${labelX} ${labelY})" text-anchor="end">${escapeXml(mutation.raw)}</text>
        `);
        prevX = x;
      });
      rows.push(`
        <g>
          <text class="hiv-mutation-gene-title" x="18" y="${top + 8}">${escapeXml(track.title)}</text>
          <line class="hiv-mutation-axis" x1="${axisStart}" y1="${axisY}" x2="${axisEnd}" y2="${axisY}"></line>
          ${tickMarkup}
          ${breakMarkup}
          <rect class="hiv-mutation-track" x="${axisStart}" y="${barY}" width="${usableWidth}" height="${barHeight}" rx="4" ry="4" fill="${track.color}"></rect>
          <text class="hiv-mutation-track-label" x="${axisStart + 10}" y="${barY + 18}" fill="${track.labelColor}">${escapeXml(track.shortLabel)}</text>
          ${mutationMarkup.join("")}
        </g>
      `);
    });
    const svgMarkup = `
      <svg class="hiv-mutation-map-svg" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${totalHeight}" role="img" aria-label="HIV mutation map" style="font-family:${escapeXml(fontFamily)};">
        <rect x="0" y="0" width="${width}" height="${totalHeight}" fill="#ffffff"></rect>
        ${rows.join("")}
      </svg>
    `;
    svgWrap.innerHTML = svgMarkup;
    if (downloadButton instanceof HTMLButtonElement) {
      downloadButton.onclick = () => downloadSvgMarkup("hiv-mutation-map.svg", svgMarkup);
    }
  };

  if (collapseInput instanceof HTMLInputElement) {
    collapseInput.addEventListener("change", () => {
      state.collapsed = collapseInput.checked;
      if (state.collapsed && !state.activeGene) state.activeGene = tracks[0]?.gene || "";
      render();
    });
  }
  render();
}

function renderHivResistanceWorkspace({
  workspaceId,
  rawTableId,
  table,
  payload,
  mutationPanels,
}) {
  const workspace = document.getElementById(workspaceId);
  const rawTableNode = document.getElementById(rawTableId);
  if (!(workspace instanceof HTMLElement) || !(rawTableNode instanceof HTMLElement)) return;
  ensureHivResistanceStyles();
  rawTableNode.dataset.exportTitle = "HIVDB_药物耐药解释表";
  renderInteractiveContigTable(
    rawTableNode,
    Array.isArray(table?.columns) ? table.columns : [],
    Array.isArray(table?.rows) ? table.rows : [],
    rawTableId,
  );
  const data = buildHivResistanceData(table, payload, mutationPanels);
  if (!data.classSummaries.length) {
    workspace.innerHTML = `<div class="hiv-resistance-empty">当前样本缺少可视化所需的耐药结构化数据，已保留原始药物总表供继续核查。</div>`;
    return;
  }
  const classLabels = {
    NRTI: "核苷类逆转录酶抑制剂",
    NNRTI: "非核苷类逆转录酶抑制剂",
    PI: "蛋白酶抑制剂",
    INSTI: "整合酶链转移抑制剂",
    CAI: "衣壳抑制剂",
  };
  const geneDescriptions = {
    PR: "蛋白酶突变谱及 PI 相关解释线索",
    RT: "逆转录酶突变谱及 NRTI / NNRTI 相关解释线索",
    IN: "整合酶突变谱及 INSTI 相关解释线索",
    CA: "衣壳突变谱及 CAI 相关解释线索",
  };
  const algorithmName = String(data.algorithm?.name || data.algorithm?.family || "HIVDB").trim() || "HIVDB";
  const algorithmVersion = String(data.algorithm?.version || "").trim();
  const initialClass = data.classSummaries
    .slice()
    .sort((a, b) => getHivResistanceLevelRank(b.topLevel) - getHivResistanceLevelRank(a.topLevel))[0]?.className
    || data.classSummaries[0]?.className
    || "";
  const state = { activeClass: initialClass };

  workspace.innerHTML = `
    <div class="hiv-resistance-shell">
      <div class="hiv-resistance-overview">
        ${data.classSummaries.map((item) => `
          <article class="hiv-resistance-overview-card">
            <span>${escapeHtml(item.className)}</span>
            <strong>${escapeHtml(item.topLevel)}</strong>
            <p>主导药物：${escapeHtml(item.topDrug)}${item.resistantCount ? ` · ${escapeHtml(String(item.resistantCount))} 个药物达到低度及以上耐药` : " · 当前类群未见明确耐药抬升"}</p>
          </article>
        `).join("")}
      </div>
      <div class="hiv-resistance-tabbar" role="tablist" aria-label="HIV drug resistance class switcher">
        ${data.classSummaries.map((item) => `
          <button type="button" data-hiv-resistance-class="${escapeHtml(item.className)}" class="${item.className === state.activeClass ? "active" : ""}">
            ${escapeHtml(item.className)} · ${escapeHtml(item.topLevel)}
          </button>
        `).join("")}
      </div>
      <div id="${escapeHtml(`${workspaceId}-panel`)}"></div>
      ${data.genes.length || data.sequenceAlerts.length ? `
        <div class="hiv-resistance-gene-grid">
          ${data.genes.map((entry) => `
            <article class="hiv-resistance-gene-card">
              <div>
                <h4>${escapeHtml(entry.gene)} 耐药解释</h4>
                <p>${escapeHtml(geneDescriptions[entry.gene] || `${entry.gene} 突变与解释线索`)}</p>
              </div>
              <div class="hiv-resistance-gene-meta">
                <strong>${escapeHtml(String(entry.mutations.length || 0))} 个观测突变</strong>
                <div class="hiv-resistance-mutations">
                  ${(entry.mutations.length ? entry.mutations : ["None"]).slice(0, 24).map((mutation) => `
                    <span class="hiv-resistance-mutation-pill">${escapeHtml(mutation)}</span>
                  `).join("")}
                </div>
              </div>
              ${entry.comments.length ? `
                <details>
                  <summary>查看 ${escapeHtml(String(entry.comments.length))} 条解释说明</summary>
                  <ul class="hiv-resistance-comment-list">
                    ${entry.comments.map((comment) => `
                      <li>
                        <strong>${escapeHtml(String(comment?.condition || comment?.comment_id || entry.gene))}</strong><br>
                        ${escapeHtml(String(comment?.text || ""))}
                      </li>
                    `).join("")}
                  </ul>
                </details>
              ` : `<p>当前基因没有额外的规则注释，但观测突变已保留用于后续人工复核。</p>`}
            </article>
          `).join("")}
          ${data.sequenceAlerts.length ? `
            <article class="hiv-resistance-alert-card">
              <h5>序列告警</h5>
              <p>这些告警通常意味着移码、终止密码子或测序覆盖异常，需要在做治疗判断前回看原始比对。</p>
              <ul class="hiv-resistance-comment-list">
                ${data.sequenceAlerts.map((alert) => `<li>${escapeHtml(alert)}</li>`).join("")}
              </ul>
            </article>
          ` : ""}
        </div>
      ` : ""}
    </div>
  `;

  const panelNode = document.getElementById(`${workspaceId}-panel`);
  if (!(panelNode instanceof HTMLElement)) return;

  const renderClassPanel = () => {
    const className = state.activeClass;
    const items = data.grouped.get(className) || [];
    const summary = data.classSummaries.find((item) => item.className === className);
    const levelTop = summary?.topLevel || "Susceptible";
    const relevantScores = items.filter((item) => item.triggered_rules.length || item.result_comments.length);
    panelNode.innerHTML = `
      <article class="hiv-resistance-panel">
        <div class="hiv-resistance-panel-head">
          <div>
            <h4>${escapeHtml(className)} 耐药解释</h4>
            <p>${escapeHtml(classLabels[className] || "药物耐药解释")} · ${escapeHtml(algorithmName)}${algorithmVersion ? ` ${escapeHtml(algorithmVersion)}` : ""}</p>
          </div>
          <div class="hiv-resistance-panel-meta">
            <span class="hiv-resistance-level is-${normalizeHivResistanceLevelName(levelTop)}">${escapeHtml(levelTop)}</span>
            <strong>${escapeHtml(String(items.length))} 个药物</strong>
          </div>
        </div>
        <div class="hiv-resistance-panel-grid">
          <div class="hiv-resistance-drug-list">
            ${items.map((item) => `
              <article class="hiv-resistance-drug-row">
                <div class="hiv-resistance-drug-grid">
                  <div class="hiv-resistance-drug-name">
                    <strong>${escapeHtml(item.drug)}</strong>
                    <span>${escapeHtml(item.fullname || "药物全名缺失")}</span>
                  </div>
                  <div class="hiv-resistance-drug-score">
                    <strong>${escapeHtml(item.score || "0")}</strong>
                    <span>分数 · ${escapeHtml(item.sir || "--")}</span>
                  </div>
                  <div>
                    <span class="hiv-resistance-level is-${normalizeHivResistanceLevelName(item.level_name)}">${escapeHtml(item.level_name)}</span>
                  </div>
                </div>
                ${(item.triggered_rules.length || item.result_comments.length) ? `
                  <div class="hiv-resistance-drug-details">
                    ${item.triggered_rules.length ? `
                      <details>
                        <summary>查看突变打分规则</summary>
                        <ul class="hiv-resistance-detail-list">
                          ${item.triggered_rules.map((rule) => `<li>${escapeHtml(rule)}</li>`).join("")}
                        </ul>
                      </details>
                    ` : ""}
                    ${item.result_comments.length ? `
                      <details>
                        <summary>查看结果注释</summary>
                        <ul class="hiv-resistance-detail-list">
                          ${item.result_comments.map((comment) => `<li>${escapeHtml(comment)}</li>`).join("")}
                        </ul>
                      </details>
                    ` : ""}
                  </div>
                ` : ""}
              </article>
            `).join("")}
          </div>
          <aside class="hiv-resistance-side">
            <article class="hiv-resistance-score-card">
              <h5>突变打分摘要</h5>
              ${relevantScores.length ? `
                <div class="hiv-resistance-score-list">
                  ${relevantScores.map((item) => `
                    <div class="hiv-resistance-score-item">
                      <strong>${escapeHtml(item.drug)} · ${escapeHtml(item.level_name)}</strong>
                      <p>${item.triggered_rules.length ? escapeHtml(item.triggered_rules.slice(0, 2).join(" / ")) : "当前药物没有额外的触发规则说明。"}${item.result_comments.length ? `<br>${escapeHtml(item.result_comments.slice(0, 1).join(" "))}` : ""}</p>
                    </div>
                  `).join("")}
                </div>
              ` : `<p>这一类药物当前没有额外的打分展开项，主体判断主要集中在最终分数和耐药等级。</p>`}
            </article>
          </aside>
        </div>
      </article>
    `;
  };

  workspace.querySelectorAll("[data-hiv-resistance-class]").forEach((button) => {
    button.addEventListener("click", () => {
      const nextClass = String(button.getAttribute("data-hiv-resistance-class") || "").trim();
      if (!nextClass || nextClass === state.activeClass) return;
      state.activeClass = nextClass;
      workspace.querySelectorAll("[data-hiv-resistance-class]").forEach((node) => {
        node.classList.toggle("active", String(node.getAttribute("data-hiv-resistance-class") || "").trim() === state.activeClass);
      });
      renderClassPanel();
    });
  });
  renderClassPanel();
}

function buildInfluenzaIgvSrcdoc(task, igvConfig, initialLocus = "") {
  const taskId = String(task?.id || "").trim();
  const referenceUrl = buildReportAssetUrl(taskId, igvConfig.reference_asset);
  const referenceIndexUrl = buildReportAssetUrl(taskId, igvConfig.reference_index_asset);
  const bamUrl = buildReportAssetUrl(taskId, igvConfig.bam_asset);
  const bamIndexUrl = buildReportAssetUrl(taskId, igvConfig.bam_index_asset);
  const gffUrl = igvConfig.gff_asset ? buildReportAssetUrl(taskId, igvConfig.gff_asset) : "";
  const alignmentVisibilityWindow = Number(igvConfig?.alignment_visibility_window || 5000);
  const alignmentSamplingWindowSize = Number(igvConfig?.alignment_sampling_window_size || 50);
  const alignmentSamplingDepth = Number(igvConfig?.alignment_sampling_depth || 50);
  const alignmentDownsampleReads = String(igvConfig?.alignment_downsample_reads ?? "true") !== "false";
  const alignmentHeight = Number(igvConfig?.alignment_height || 220);
  return `<!DOCTYPE html>
<html lang="zh-CN">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IGV 比对结果</title>
    <style>
      html, body { margin: 0; padding: 0; height: 100%; background: #ffffff; font-family: sans-serif; }
      #igv-div { height: 100vh; width: 100%; }
    </style>
    <script src="/public/igv/igv.min.js"></script>
  </head>
  <body>
    <div id="igv-div"></div>
    <script>
      const initialLocus = ${JSON.stringify(String(initialLocus || ""))};
      const options = {
        showNavigation: true,
        showRuler: true,
        loadDefaultGenomes: false,
        reference: {
          fastaURL: ${JSON.stringify(referenceUrl)},
          indexURL: ${JSON.stringify(referenceIndexUrl)}
        },
        tracks: [
          {
            name: ${JSON.stringify(String(igvConfig.viewer_label || "参考比对视图"))},
            type: "alignment",
            format: "bam",
            url: ${JSON.stringify(bamUrl)},
            indexURL: ${JSON.stringify(bamIndexUrl)},
            height: ${JSON.stringify(alignmentHeight)},
            displayMode: "SQUISHED",
            visibilityWindow: ${JSON.stringify(alignmentVisibilityWindow)},
            samplingWindowSize: ${JSON.stringify(alignmentSamplingWindowSize)},
            samplingDepth: ${JSON.stringify(alignmentSamplingDepth)},
            downsampleReads: ${JSON.stringify(alignmentDownsampleReads)},
            showSoftClips: false,
            showInsertionText: false,
            showMismatches: false
          }${gffUrl ? `,
          {
            name: "VADR 注释",
            type: "annotation",
            format: "gff3",
            url: ${JSON.stringify(gffUrl)},
            displayMode: "EXPANDED",
            visibilityWindow: 300000000
          }` : ""}
        ]
      };
      const igvDiv = document.getElementById("igv-div");
      let browser = null;
      function jump(target) {
        const locus = String(target || "").trim();
        if (!locus) return false;
        const chrom = locus.split(":")[0]?.trim();
        if (chrom && browser?.genome && typeof browser.genome.getChromosome === "function" && !browser.genome.getChromosome(chrom)) {
          try {
            window.parent.postMessage({
              type: "influenza-igv-debug",
              stage: "locus_not_in_reference",
              locus,
              detail: "当前 IGV 参考序列不包含该 chromosome",
            }, "*");
          } catch (error) {}
          return false;
        }
        try {
          window.parent.postMessage({
            type: "influenza-igv-debug",
            stage: "jump_called",
            locus,
            detail: "iframe 已收到跳转请求",
          }, "*");
        } catch (error) {}
        if (browser && typeof browser.search === "function") {
          browser.search(locus);
          try {
            window.parent.postMessage({
              type: "influenza-igv-debug",
              stage: "browser_search",
              locus,
              detail: "iframe 已调用 browser.search()",
            }, "*");
          } catch (error) {}
          return true;
        }
        const inputNode = document.getElementsByClassName("igv-search-input")[0];
        const searchBtn = document.getElementsByClassName("igv-search-icon-container")[0];
        if (inputNode) {
          inputNode.value = locus;
          inputNode.dispatchEvent(new Event("input", { bubbles: true }));
          inputNode.dispatchEvent(new Event("change", { bubbles: true }));
        }
        if (searchBtn) {
          searchBtn.click();
          try {
            window.parent.postMessage({
              type: "influenza-igv-debug",
              stage: "search_button_click",
              locus,
              detail: "iframe 已触发搜索按钮点击",
            }, "*");
          } catch (error) {}
          return true;
        }
        return false;
      }
      window.jump = jump;
      window.addEventListener("message", (event) => {
        if (event?.data?.value) {
          try {
            window.parent.postMessage({
              type: "influenza-igv-debug",
              stage: "message_received",
              locus: String(event.data.value || "").trim(),
              detail: "iframe 已收到父页面 postMessage",
            }, "*");
          } catch (error) {}
          jump(event.data.value);
        }
      });
      igv.createBrowser(igvDiv, options).then((createdBrowser) => {
        browser = createdBrowser;
        window.browser = createdBrowser;
        try {
          window.parent.postMessage({
            type: "influenza-igv-debug",
            stage: "browser_ready",
            locus: initialLocus || "",
            detail: "IGV browser 已初始化完成",
          }, "*");
        } catch (error) {}
        if (initialLocus) jump(initialLocus);
      });
    </script>
  </body>
</html>`;
}

function initializeDeferredIgvEmbed({
  frameId,
  panelId,
  buttonId,
  task,
  igvView,
  initialLocus = "",
  mutationTableNode = null,
}) {
  const igvFrame = document.getElementById(frameId);
  const lazyPanel = document.getElementById(panelId);
  const loadButton = document.getElementById(buttonId);
  if (!(igvFrame instanceof HTMLIFrameElement) || !lazyPanel || !(loadButton instanceof HTMLButtonElement)) return;
  let requested = false;
  let loaded = false;

  const revealLoadedFrame = () => {
    lazyPanel.hidden = true;
    igvFrame.hidden = false;
    loadButton.disabled = false;
    loadButton.textContent = "重新加载 IGV";
  };

  const requestLoad = (targetLocus = "") => {
    const locus = String(targetLocus || initialLocus || "").trim();
    if (!requested) {
      requested = true;
      igvFrame.hidden = true;
      loadButton.disabled = true;
      loadButton.textContent = "IGV 加载中...";
      igvFrame.srcdoc = buildInfluenzaIgvSrcdoc(task || {}, igvView || {}, locus);
      return true;
    }
    if (loaded && locus && igvFrame.contentWindow) {
      try {
        igvFrame.contentWindow.postMessage({ value: locus }, window.location.origin);
      } catch (error) {}
    }
    return false;
  };

  igvFrame.addEventListener("load", () => {
    loaded = true;
    revealLoadedFrame();
  });

  loadButton.addEventListener("click", () => {
    if (loaded) {
      requested = false;
      loaded = false;
      lazyPanel.hidden = false;
      requestLoad(initialLocus);
      return;
    }
    requestLoad(initialLocus);
  });

  if (mutationTableNode instanceof HTMLElement) {
    initializeInfluenzaIgvLink(mutationTableNode, igvFrame, requestLoad);
  }
}

function showInfluenzaIgvToast(message, kind = "info") {
  let node = document.getElementById("influenza-igv-debug-toast");
  if (!node) {
    node = document.createElement("div");
    node.id = "influenza-igv-debug-toast";
    node.className = "igv-debug-toast";
    node.setAttribute("aria-live", "polite");
    document.body.appendChild(node);
  }
  node.dataset.kind = kind;
  node.textContent = message;
  node.classList.add("is-visible");
  window.clearTimeout(showInfluenzaIgvToast.timer);
  showInfluenzaIgvToast.timer = window.setTimeout(() => {
    node.classList.remove("is-visible");
  }, 3200);
}

function updateInfluenzaIgvDebug(patch = {}) {
  const clickedNode = document.getElementById("influenza-igv-debug-clicked");
  const sentNode = document.getElementById("influenza-igv-debug-sent");
  const receivedNode = document.getElementById("influenza-igv-debug-received");
  const statusNode = document.getElementById("influenza-igv-debug-status");
  const timeNode = document.getElementById("influenza-igv-debug-updated");
  if (!clickedNode || !sentNode || !receivedNode || !statusNode) return;
  if (Object.prototype.hasOwnProperty.call(patch, "clicked")) {
    clickedNode.textContent = String(patch.clicked || "--");
  }
  if (Object.prototype.hasOwnProperty.call(patch, "sent")) {
    sentNode.textContent = String(patch.sent || "--");
  }
  if (Object.prototype.hasOwnProperty.call(patch, "received")) {
    receivedNode.textContent = String(patch.received || "--");
  }
  if (Object.prototype.hasOwnProperty.call(patch, "status")) {
    statusNode.textContent = String(patch.status || "--");
  }
  if (timeNode) {
    timeNode.textContent = new Date().toLocaleTimeString("zh-CN", { hour12: false });
  }
}

function extractInfluenzaLocusFromRow(row) {
  if (!(row instanceof HTMLTableRowElement)) return "";
  const locus = String(row.dataset.igvLocus || "").trim();
  return locus;
}

window.__influenzaIgvSelectLocusImpl = null;
window.__influenzaIgvSelectLocus = function(locusValue) {
  const locus = String(locusValue || "").trim();
  if (!locus) return false;
  updateInfluenzaIgvDebug({
    clicked: locus,
    status: `已捕获表格点击：${locus}`,
  });
  if (typeof window.__influenzaIgvSelectLocusImpl === "function") {
    return window.__influenzaIgvSelectLocusImpl(locus);
  }
  showInfluenzaIgvToast(`已记录位点，但 IGV 还未完成绑定：${locus}`, "warn");
  return true;
};

function initializeInfluenzaIgvLink(mutationTableNode, igvFrame, ensureFrameReady = null) {
  if (!(mutationTableNode instanceof HTMLElement) || !(igvFrame instanceof HTMLIFrameElement)) return;
  if (mutationTableNode.dataset.igvLinkBound === "true") return;
  mutationTableNode.dataset.igvLinkBound = "true";
  let pendingTarget = "";
  updateInfluenzaIgvDebug({
    clicked: "--",
    sent: "--",
    received: "--",
    status: "等待点击突变位点",
  });
  const jumpTo = (target, attempt = 0) => {
    const locus = String(target || "").trim();
    if (!locus) return false;
    updateInfluenzaIgvDebug({
      sent: locus,
      status: `正在向 IGV 发送位点：${locus}${attempt > 0 ? `（重试 ${attempt}）` : ""}`,
    });
    showInfluenzaIgvToast(`正在跳转 IGV：${locus}`);
    const win = igvFrame.contentWindow;
    const doc = win?.document;
    if (win && typeof win.jump === "function") {
      try {
        const jumped = win.jump(locus);
        updateInfluenzaIgvDebug({
          received: jumped ? locus : "--",
          status: jumped ? `已调用 iframe 内 jump()：${locus}` : `IGV 当前参考不包含位点：${locus}`,
        });
        showInfluenzaIgvToast(jumped ? `已直接调用 IGV jump()：${locus}` : `IGV 当前参考不包含位点：${locus}`, jumped ? "success" : "warn");
        return true;
      } catch (error) {
        updateInfluenzaIgvDebug({
          status: `IGV jump() 调用失败，改用备用方式：${error?.message || error}`,
        });
        showInfluenzaIgvToast(`IGV jump() 调用失败，改用备用方式：${error?.message || error}`, "warn");
      }
    }
    if (win?.browser && typeof win.browser.search === "function") {
      try {
        const result = win.browser.search(locus);
        if (result && typeof result.then === "function") {
          result.then(() => {
            updateInfluenzaIgvDebug({
              received: locus,
              status: `IGV browser.search() 已执行：${locus}`,
            });
            showInfluenzaIgvToast(`已调用 IGV browser.search()：${locus}`, "success");
          }).catch((error) => {
            updateInfluenzaIgvDebug({
              status: `IGV browser.search() 失败：${error?.message || error}`,
            });
            showInfluenzaIgvToast(`IGV browser.search() 失败：${error?.message || error}`, "warn");
          });
        } else {
          updateInfluenzaIgvDebug({
            received: locus,
            status: `IGV browser.search() 已执行：${locus}`,
          });
          showInfluenzaIgvToast(`已调用 IGV browser.search()：${locus}`, "success");
        }
        return true;
      } catch (error) {
        updateInfluenzaIgvDebug({
          status: `IGV browser.search() 调用失败，改用备用方式：${error?.message || error}`,
        });
        showInfluenzaIgvToast(`IGV browser.search() 调用失败，改用备用方式：${error?.message || error}`, "warn");
      }
    }
    const inputNode = doc?.getElementsByClassName("igv-search-input")?.[0];
    const searchBtn = doc?.getElementsByClassName("igv-search-icon-container")?.[0];
    if (inputNode instanceof HTMLInputElement) {
      inputNode.value = locus;
      inputNode.dispatchEvent(new Event("input", { bubbles: true }));
      inputNode.dispatchEvent(new Event("change", { bubbles: true }));
      inputNode.dispatchEvent(new Event("keydown", { bubbles: true, key: "Enter", code: "Enter" }));
    }
    if (searchBtn instanceof HTMLElement) {
      searchBtn.click();
      searchBtn.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, view: window }));
      updateInfluenzaIgvDebug({
        received: locus,
        status: `已触发 IGV 搜索按钮：${locus}`,
      });
      showInfluenzaIgvToast(`已触发 IGV 搜索：${locus}`, "success");
      return true;
    }
    if (win) {
      win.postMessage({ value: locus }, window.location.origin);
      updateInfluenzaIgvDebug({
        status: `已向 iframe postMessage 位点：${locus}`,
      });
      showInfluenzaIgvToast(`已向 IGV 发送消息：${locus}`, "success");
      if (attempt < 6) {
        window.setTimeout(() => {
          jumpTo(locus, attempt + 1);
        }, 250 + attempt * 120);
      }
      return true;
    }
    if (attempt < 6) {
      window.setTimeout(() => jumpTo(locus, attempt + 1), 250 + attempt * 120);
      return true;
    }
    updateInfluenzaIgvDebug({
      status: "IGV 还没准备好，稍后会自动再试",
    });
    showInfluenzaIgvToast("IGV 还没准备好，稍后会自动再试", "warn");
    return false;
  };

  igvFrame.addEventListener("load", () => {
    if (pendingTarget) {
      const target = pendingTarget;
      pendingTarget = "";
      updateInfluenzaIgvDebug({
        status: `IGV iframe 已加载，准备定位：${target}`,
      });
      showInfluenzaIgvToast(`IGV 已加载，准备定位：${target}`);
      jumpTo(target);
    } else {
      updateInfluenzaIgvDebug({
        status: "IGV iframe 已加载，等待用户点击位点",
      });
      showInfluenzaIgvToast("IGV 已加载", "success");
    }
  });

  const handleSelect = (locusValue) => {
    const locus = String(locusValue || "").trim();
    if (!locus) return false;
    pendingTarget = locus;
    updateInfluenzaIgvDebug({
      clicked: locus,
      status: `已捕获表格点击：${locus}`,
    });
    showInfluenzaIgvToast(`已点选位点：${locus}`);
    if (typeof ensureFrameReady === "function") {
      ensureFrameReady(locus);
    }
    if (igvFrame.srcdoc || igvFrame.src) {
      jumpTo(locus);
    } else {
      updateInfluenzaIgvDebug({
        sent: locus,
        status: `正在按需加载 IGV：${locus}`,
      });
      showInfluenzaIgvToast(`正在加载 IGV：${locus}`);
    }
    return true;
  };
  window.__influenzaIgvSelectLocusImpl = handleSelect;

  window.addEventListener("message", (event) => {
    if (event.source !== igvFrame.contentWindow) return;
    const payload = event?.data;
    if (!payload || payload.type !== "influenza-igv-debug") return;
    updateInfluenzaIgvDebug({
      received: payload.locus || payload.stage || "--",
      status: payload.detail || payload.stage || "IGV 已返回调试信息",
    });
  });

  mutationTableNode.classList.add("influenza-igv-linked-table");
  mutationTableNode.addEventListener("click", (event) => {
    const target = event.target instanceof Element ? event.target : null;
    const row = target?.closest("tr[data-igv-locus]");
    if (!(row instanceof HTMLTableRowElement)) return;
    const locus = String(row.dataset.igvLocus || "").trim();
    if (!locus) return;
    handleSelect(locus);
  }, true);
}

async function renderCgviewSection(task, section, sections = {}) {
  renderMiniStatGrid("cgview-summary", [
    { label: "图谱数量", value: String(section?.summary?.map_count ?? 0) },
    { label: "主图状态", value: Array.isArray(section?.maps) && section.maps.some((item) => item.role === "main") ? "已就绪" : "未提供" },
  ]);
  const container = document.getElementById("cgview-viewer-card");
  if (!container) return;
  container.classList.remove("empty-box");
  const maps = Array.isArray(section?.maps) ? section.maps : [];
  if (!maps.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>当前没有可用于 CGView 的 GenBank 文件</strong>
        <p class="empty-copy">请确认单菌流程已产出 <code>{Sam}_prokka/main.gbk</code> 或对应质粒 <code>.gbk</code> 文件，然后这里会自动生成交互环形图。</p>
      </div>
    `;
    return;
  }
  const runtime = window.CGView || window.CGV;
  const parserRuntime = window.CGParse;
  if (!runtime?.Viewer || !parserRuntime?.CGViewBuilder) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>CGView 组件未加载成功</strong>
        <p class="empty-copy">请确认已将 <code>public/cgview/d3.min.js</code>、<code>public/cgview/cgparse.min.js</code> 和 <code>public/cgview/cgview.min.js</code> 放入项目，再刷新结果页。</p>
      </div>
    `;
    return;
  }

  const initialKey = maps.find((item) => item.role === "main")?.key || maps[0].key;
  container.innerHTML = `
    <div class="cgview-shell">
      <div class="cgview-toolbar">
        <div class="cgview-picker" id="cgview-picker">
          <span class="cgview-picker-label">当前图谱</span>
          <button type="button" id="cgview-picker-toggle" class="cgview-picker-toggle" aria-haspopup="listbox" aria-expanded="false">
            <span id="cgview-picker-value">${escapeHtml(maps.find((item) => item.key === initialKey)?.label || "")}</span>
          </button>
          <div id="cgview-picker-panel" class="cgview-picker-panel" hidden>
            <label class="cgview-picker-search">
              <input id="cgview-picker-search" type="search" placeholder="搜索 contig / 质粒">
            </label>
            <div id="cgview-picker-options" class="cgview-picker-options" role="listbox" aria-label="CGView 图谱列表"></div>
          </div>
        </div>
        <span id="cgview-map-note" class="card-tag">${escapeHtml(maps.find((item) => item.key === initialKey)?.asset_name || "")}</span>
      </div>
      <div id="cgview-stage" class="cgview-stage"></div>
      <div class="cgview-overlay-panel" id="cgview-overlay-panel">
        <div class="cgview-overlay-header">
          <span class="cgview-overlay-label">外层标记</span>
          <div class="cgview-overlay-actions">
            <button type="button" class="cgview-overlay-action is-active" data-overlay-action="labels">全部基因名</button>
            <button type="button" class="cgview-overlay-action" data-overlay-action="all">全选</button>
            <button type="button" class="cgview-overlay-action" data-overlay-action="none">清空</button>
          </div>
        </div>
        <div class="cgview-overlay-toggles" id="cgview-overlay-toggles">
          <button type="button" class="cgview-overlay-toggle is-active" data-overlay-source="portal-resistance">耐药基因</button>
          <button type="button" class="cgview-overlay-toggle is-active" data-overlay-source="portal-virulence">毒力基因</button>
          <button type="button" class="cgview-overlay-toggle is-active" data-overlay-source="portal-mge">移动元件</button>
        </div>
      </div>
      <div class="cgview-controls" id="cgview-controls">
        <button type="button" class="cgview-control-button" data-cgview-action="zoom-in">放大</button>
        <button type="button" class="cgview-control-button" data-cgview-action="zoom-out">缩小</button>
        <button type="button" class="cgview-control-button" data-cgview-action="reset">重置</button>
        <button type="button" class="cgview-control-button" data-cgview-action="download-png">导出 PNG</button>
        <button type="button" class="cgview-control-button" data-cgview-action="download-svg">导出 SVG</button>
      </div>
    </div>
  `;

  const stage = document.getElementById("cgview-stage");
  const note = document.getElementById("cgview-map-note");
  const picker = document.getElementById("cgview-picker");
  const pickerToggle = document.getElementById("cgview-picker-toggle");
  const pickerValue = document.getElementById("cgview-picker-value");
  const pickerPanel = document.getElementById("cgview-picker-panel");
  const pickerSearch = document.getElementById("cgview-picker-search");
  const pickerOptions = document.getElementById("cgview-picker-options");
  const controls = document.getElementById("cgview-controls");
  const overlayPanel = document.getElementById("cgview-overlay-panel");
  const overlayToggles = document.getElementById("cgview-overlay-toggles");
  if (!stage || !picker || !pickerToggle || !pickerValue || !pickerPanel || !pickerSearch || !pickerOptions || !controls || !overlayPanel || !overlayToggles) return;

  let activeKey = initialKey;
  let activeViewer = null;
  const overlaySources = ["portal-resistance", "portal-virulence", "portal-mge"];
  const activeOverlaySources = new Set(overlaySources);
  let showAllGeneLabels = true;

  const setPickerOpen = (open) => {
    picker.classList.toggle("is-open", open);
    pickerPanel.hidden = !open;
    pickerToggle.setAttribute("aria-expanded", open ? "true" : "false");
    if (open) {
      pickerSearch.focus();
      pickerSearch.select();
    }
  };

  const renderPickerOptions = (keyword = "") => {
    const normalizedKeyword = String(keyword || "").trim().toLowerCase();
    const visibleMaps = maps.filter((item) => {
      if (!normalizedKeyword) return true;
      const haystack = `${item.label} ${item.asset_name} ${item.key}`.toLowerCase();
      return haystack.includes(normalizedKeyword);
    });
    pickerOptions.innerHTML = visibleMaps.length ? visibleMaps.map((item) => `
      <button
        type="button"
        class="cgview-picker-option${item.key === activeKey ? " is-active" : ""}"
        data-map-key="${escapeHtml(item.key)}"
        role="option"
        aria-selected="${item.key === activeKey ? "true" : "false"}"
      >
        <strong>${escapeHtml(item.label)}</strong>
        <span>${escapeHtml(item.asset_name || "")}</span>
      </button>
    `).join("") : `<p class="cgview-picker-empty">没有匹配的图谱</p>`;
  };

  const syncOverlayButtons = (availableSources = overlaySources) => {
    const available = new Set(availableSources);
    overlayPanel.querySelectorAll('[data-overlay-action="labels"]').forEach((button) => {
      if (!(button instanceof HTMLElement)) return;
      button.classList.toggle("is-active", showAllGeneLabels);
    });
    overlayToggles.querySelectorAll("[data-overlay-source]").forEach((button) => {
      if (!(button instanceof HTMLElement)) return;
      const source = button.dataset.overlaySource || "";
      const enabled = available.has(source);
      button.classList.toggle("is-active", enabled && activeOverlaySources.has(source));
      button.toggleAttribute("disabled", !enabled);
    });
    overlayPanel.classList.toggle("is-empty", !available.size);
  };

  const cgviewBuilderConfig = {
    settings: {
      format: "circular",
      backgroundColor: "rgb(255,255,255)",
      geneticCode: 11,
    },
    legend: {
      visible: true,
    },
  };

  const drawMap = async (mapKey) => {
    const currentMap = maps.find((item) => item.key === mapKey) || maps[0];
    if (!currentMap) return;
    activeKey = currentMap.key;
    pickerValue.textContent = currentMap.label || "";
    renderPickerOptions(pickerSearch.value);
    if (note) note.textContent = currentMap.asset_name || "";
    stage.innerHTML = `<div class="empty-copy">正在载入 ${escapeHtml(currentMap.label)}...</div>`;
    try {
      let payload = null;
      let builderError = null;
      try {
        const gbkText = await fetchReportAsset(task?.id || "", currentMap.asset_name, "text");
        const builder = new parserRuntime.CGViewBuilder(gbkText, {
          config: cgviewBuilderConfig,
          excludeFeatures: ["source", "gene", "exon"],
          excludeQualifiers: ["translation"],
        });
        payload = normalizeCgviewPayload(
          appendCgviewOverlayFeatures(builder.toJSON(), currentMap, sections),
          currentMap.label,
        );
        if (!hasRenderableCgviewPayload(payload)) {
          throw new Error("CGViewBuilder 返回的图谱缺少可渲染的图例或轨道配置");
        }
      } catch (error) {
        builderError = error;
      }

      if (!payload) {
        payload = normalizeCgviewPayload(
          appendCgviewOverlayFeatures(
            await fetchReportAsset(task?.id || "", currentMap.asset_name, "json", { renderAs: "cgview-json" }),
            currentMap,
            sections,
          ),
          currentMap.label,
        );
        if (note) note.textContent = `${currentMap.asset_name} · fallback`;
      }

      const mountId = `cgview-canvas-${Date.now()}`;
      stage.innerHTML = `<div id="${mountId}" class="cgview-canvas"></div>`;
      const size = Math.max(640, Math.min(960, Math.round(stage.clientWidth || 820)));
      const viewer = new runtime.Viewer(`#${mountId}`, {
        width: size,
        height: size,
      });
      viewer.io.loadJSON(payload);
      syncOverlayButtons(getCgviewOverlayTrackSources(viewer));
      if (viewer.annotation?.update) {
        viewer.annotation.update({
          visible: true,
          onlyDrawFavorites: !showAllGeneLabels,
        });
      }
      applyCgviewOverlayVisibility(viewer, activeOverlaySources);
      activeViewer = viewer;
      if (builderError) {
        console.warn("CGViewBuilder 构建失败，已回退到服务端转换：", builderError);
      }
    } catch (error) {
      stage.innerHTML = `
        <div class="empty-table-state">
          <strong>CGView 图谱加载失败</strong>
          <p class="empty-copy">${escapeHtml(error?.message || "无法根据当前 GenBank 文件绘制环形图。")}</p>
        </div>
      `;
    }
  };

  pickerToggle.addEventListener("click", () => {
    setPickerOpen(!picker.classList.contains("is-open"));
  });
  pickerSearch.addEventListener("input", () => {
    renderPickerOptions(pickerSearch.value);
  });
  pickerOptions.addEventListener("click", (event) => {
    const option = event.target.closest("[data-map-key]");
    if (!(option instanceof HTMLElement)) return;
    setPickerOpen(false);
    drawMap(option.dataset.mapKey || "");
  });
  document.addEventListener("click", (event) => {
    if (!picker.contains(event.target)) {
      setPickerOpen(false);
    }
  });
  controls.addEventListener("click", (event) => {
    const button = event.target.closest("[data-cgview-action]");
    if (!(button instanceof HTMLElement) || !activeViewer) return;
    const currentMap = maps.find((item) => item.key === activeKey) || maps[0];
    const safeName = String(currentMap?.key || "cgview").replace(/[^a-zA-Z0-9._-]+/g, "-");
    switch (button.dataset.cgviewAction) {
      case "zoom-in":
        activeViewer.zoomIn(1.5, { duration: 220 });
        break;
      case "zoom-out":
        activeViewer.zoomOut(1.5, { duration: 220 });
        break;
      case "reset":
        activeViewer.reset(240);
        break;
      case "download-png":
        activeViewer.io.downloadImage(1800, 1800, `${safeName}.png`);
        break;
      case "download-svg":
        activeViewer.io.downloadSVG(`${safeName}.svg`);
        break;
      default:
        break;
    }
  });
  overlayPanel.addEventListener("click", (event) => {
    const actionButton = event.target.closest("[data-overlay-action]");
    if (actionButton instanceof HTMLElement) {
      if (actionButton.dataset.overlayAction === "labels") {
        showAllGeneLabels = !showAllGeneLabels;
        syncOverlayButtons(activeViewer ? getCgviewOverlayTrackSources(activeViewer) : overlaySources);
        if (activeViewer?.annotation?.update) {
          activeViewer.annotation.update({
            visible: true,
            onlyDrawFavorites: !showAllGeneLabels,
          });
          activeViewer.drawFull();
        }
      } else if (actionButton.dataset.overlayAction === "all") {
        overlaySources.forEach((source) => activeOverlaySources.add(source));
        syncOverlayButtons(activeViewer ? getCgviewOverlayTrackSources(activeViewer) : overlaySources);
        if (activeViewer) applyCgviewOverlayVisibility(activeViewer, activeOverlaySources);
      } else if (actionButton.dataset.overlayAction === "none") {
        activeOverlaySources.clear();
        syncOverlayButtons(activeViewer ? getCgviewOverlayTrackSources(activeViewer) : overlaySources);
        if (activeViewer) applyCgviewOverlayVisibility(activeViewer, activeOverlaySources);
      }
      return;
    }
    const toggle = event.target.closest("[data-overlay-source]");
    if (!(toggle instanceof HTMLElement) || toggle.hasAttribute("disabled")) return;
    const source = toggle.dataset.overlaySource || "";
    if (!source) return;
    if (activeOverlaySources.has(source)) {
      activeOverlaySources.delete(source);
    } else {
      activeOverlaySources.add(source);
    }
    syncOverlayButtons(activeViewer ? getCgviewOverlayTrackSources(activeViewer) : overlaySources);
    if (activeViewer) applyCgviewOverlayVisibility(activeViewer, activeOverlaySources);
  });

  renderPickerOptions();
  syncOverlayButtons();
  await drawMap(initialKey);
}

function renderPathoSourceDistanceBins(containerId, section) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const bars = Array.isArray(section?.bars) ? section.bars : [];
  if (!bars.length) {
    container.innerHTML = `<p class="empty-copy">未读取到距离分布数据。</p>`;
    return;
  }
  container.innerHTML = `
    <div class="patho-distribution-wrap">
      <div class="patho-distribution-total">总比较对数：<strong>${escapeHtml(String(section?.total_pairs ?? 0))}</strong></div>
      <div class="patho-distribution-bars">
        ${bars.map((item) => `
          <article class="patho-bar-item">
            <div class="patho-bar-head">
              <span>${escapeHtml(item.label)}</span>
              <strong>${escapeHtml(String(item.value))}</strong>
            </div>
            <div class="patho-bar-track" aria-hidden="true">
              <div class="patho-bar-fill" style="width:${Math.max(0, Math.min(100, Number(item.ratio) || 0))}%"></div>
            </div>
            <div class="patho-bar-meta">${escapeHtml(String(item.share))}%</div>
          </article>
        `).join("")}
      </div>
    </div>
  `;
}

function buildHeatValueColor(value, minValue, maxValue, hue) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "rgba(62, 84, 111, 0.06)";
  if (maxValue <= minValue) return `hsla(${hue}, 52%, 54%, 0.34)`;
  const ratio = Math.max(0, Math.min(1, (numeric - minValue) / (maxValue - minValue)));
  const lightness = 96 - ratio * 48;
  const alpha = 0.16 + ratio * 0.7;
  return `hsla(${hue}, 58%, ${lightness}%, ${alpha})`;
}

function buildHeatTextColor(value, minValue, maxValue) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "var(--report-muted)";
  if (maxValue <= minValue) return "var(--report-ink)";
  const ratio = Math.max(0, Math.min(1, (numeric - minValue) / (maxValue - minValue)));
  return ratio >= 0.52 ? "#f8f5ef" : "var(--report-ink)";
}

function parseHeatmapPreview(section) {
  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  if (columns.length < 2 || !rows.length) {
    return { samples: [], cells: [], min: 0, max: 0 };
  }
  const samples = columns.slice(1);
  const cells = rows.map((row) => {
    const sample = String(row?.[0] ?? "").trim();
    const values = samples.map((target, index) => {
      const raw = row?.[index + 1];
      const numeric = Number(raw);
      return {
        source: sample,
        target,
        display: raw == null || raw === "" ? "-" : String(raw),
        value: Number.isFinite(numeric) ? numeric : null,
      };
    });
    return { sample, values };
  });
  const numericValues = cells.flatMap((row) => row.values.map((item) => item.value)).filter((value) => Number.isFinite(value));
  return {
    samples,
    cells,
    min: numericValues.length ? Math.min(...numericValues) : 0,
    max: numericValues.length ? Math.max(...numericValues) : 0,
  };
}

function normalizePathoSampleName(value) {
  const raw = String(value || "").trim();
  if (!raw) return "";
  return raw
    .replace(/\.(raw\.)?f(ast)?a(stq)?$/i, "")
    .replace(/\.f(ast)?q(\.gz)?$/i, "")
    .replace(/\.f(ast)?a$/i, "")
    .replace(/^\s+|\s+$/g, "");
}

function extractPathoSampleNames(value) {
  if (Array.isArray(value)) {
    return Array.from(new Set(value.map((item) => normalizePathoSampleName(item)).filter(Boolean)));
  }
  const text = String(value || "").trim();
  if (!text) return [];
  return Array.from(
    new Set(
      text
        .split(/[\n,;/|]+/)
        .map((item) => normalizePathoSampleName(item))
        .filter(Boolean),
    ),
  );
}

function encodePathoSampleGroup(samples) {
  return JSON.stringify(extractPathoSampleNames(samples));
}

function decodePathoSampleGroup(rawValue) {
  if (!rawValue) return [];
  try {
    const parsed = JSON.parse(rawValue);
    return extractPathoSampleNames(parsed);
  } catch (_error) {
    return extractPathoSampleNames(rawValue);
  }
}

let activePathoSampleNames = [];
let pathoSampleBindingsReady = false;

function clearPathoSampleHighlights() {
  document.querySelectorAll(".is-patho-sample-focus").forEach((node) => node.classList.remove("is-patho-sample-focus"));
}

function pushSampleHighlightToFrame(frame, messageType, samples) {
  if (!frame) return;
  try {
    const frameWindow = frame.contentWindow;
    if (!frameWindow) return;
    if (typeof frameWindow.postMessage === "function") {
      frameWindow.postMessage({ type: messageType, samples }, "*");
    }
    if (typeof frameWindow.setPortalHighlightedSamples === "function") {
      frameWindow.setPortalHighlightedSamples(samples);
    }
  } catch (_error) {
    // ignore cross-frame highlight errors
  }
}

function syncPathoTreeHighlights(samples) {
  document.querySelectorAll(".patho-grapetree-frame").forEach((frame) => {
    const title = String(frame.getAttribute("title") || "");
    if (/GrapeTree/i.test(title)) {
      pushSampleHighlightToFrame(frame, "portal-grapetree-highlight", samples);
      return;
    }
    pushSampleHighlightToFrame(frame, "portal-itol-highlight", samples);
  });
}

function applyPathoSampleHighlights(samples) {
  const normalized = extractPathoSampleNames(samples);
  activePathoSampleNames = normalized;
  clearPathoSampleHighlights();
  if (!normalized.length) {
    syncPathoTreeHighlights([]);
    return;
  }
  const sampleSet = new Set(normalized);
  document.querySelectorAll("[data-patho-sample-name]").forEach((node) => {
    const sample = normalizePathoSampleName(node.getAttribute("data-patho-sample-name"));
    if (sampleSet.has(sample)) {
      node.classList.add("is-patho-sample-focus");
    }
  });
  document.querySelectorAll("[data-patho-sample-group]").forEach((node) => {
    const group = decodePathoSampleGroup(node.getAttribute("data-patho-sample-group"));
    if (group.some((sample) => sampleSet.has(sample))) {
      node.classList.add("is-patho-sample-focus");
    }
  });
  document.querySelectorAll("[data-patho-sample-a][data-patho-sample-b]").forEach((node) => {
    const sampleA = normalizePathoSampleName(node.getAttribute("data-patho-sample-a"));
    const sampleB = normalizePathoSampleName(node.getAttribute("data-patho-sample-b"));
    if (sampleSet.has(sampleA) || sampleSet.has(sampleB)) {
      node.classList.add("is-patho-sample-focus");
    }
  });
  syncPathoTreeHighlights(normalized);
}

function bindPathoSampleInteractions() {
  if (pathoSampleBindingsReady) return;
  pathoSampleBindingsReady = true;
  const extractFromNode = (node) => {
    if (!node) return [];
    if (node.hasAttribute("data-patho-sample-group")) {
      return decodePathoSampleGroup(node.getAttribute("data-patho-sample-group"));
    }
    if (node.hasAttribute("data-patho-sample-name")) {
      return extractPathoSampleNames(node.getAttribute("data-patho-sample-name"));
    }
    if (node.hasAttribute("data-patho-sample-a") || node.hasAttribute("data-patho-sample-b")) {
      return extractPathoSampleNames([
        node.getAttribute("data-patho-sample-a"),
        node.getAttribute("data-patho-sample-b"),
      ]);
    }
    return [];
  };
  const activateFromEvent = (event) => {
    const target = event.target instanceof Element
      ? event.target.closest("[data-patho-sample-group], [data-patho-sample-name], [data-patho-sample-a][data-patho-sample-b]")
      : null;
    if (!target) return;
    const samples = extractFromNode(target);
    if (!samples.length) return;
    applyPathoSampleHighlights(samples);
  };
  document.addEventListener("mouseenter", activateFromEvent, true);
  document.addEventListener("focusin", activateFromEvent, true);
  document.addEventListener("click", activateFromEvent, true);
}

function renderHeatmapMatrix(containerId, title, section, options = {}) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const { samples, cells, min, max } = parseHeatmapPreview(section);
  if (!samples.length || !cells.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>${escapeHtml(title)} 暂无数据</strong>
        <p class="empty-copy">当前没有可用于热图展示的矩阵结果。</p>
      </div>
    `;
    return;
  }
  const hue = options.hue || 208;
  const unit = options.unit || "";
  const defaultCell = cells[0]?.values?.[0] || null;
  const defaultDetail = defaultCell
    ? `${defaultCell.source} vs ${defaultCell.target}: ${defaultCell.display}${unit}`
    : "悬停任意格点查看样本对关系。";
  container.innerHTML = `
    <div class="patho-heatmap-card">
      <div class="patho-heatmap-head">
        <div>
          <strong>${escapeHtml(title)}</strong>
          <p>悬停查看样本对的具体距离，点击后可锁定当前说明。</p>
        </div>
        <div class="patho-heatmap-scale" aria-hidden="true">
          <span>${escapeHtml(String(min))}${unit}</span>
          <div class="patho-heatmap-scale-bar patho-heatmap-scale-bar-${hue}"></div>
          <span>${escapeHtml(String(max))}${unit}</span>
        </div>
      </div>
      <div class="patho-heatmap-detail" data-heatmap-detail>${escapeHtml(defaultDetail)}</div>
      <div class="patho-heatmap-scroll">
        <div class="patho-heatmap-grid" style="grid-template-columns: minmax(120px, 160px) repeat(${samples.length}, minmax(42px, 1fr));">
          <div class="patho-heatmap-corner">样本</div>
          ${samples.map((sample) => `<div class="patho-heatmap-col-header" data-patho-sample-name="${escapeHtml(normalizePathoSampleName(sample))}" title="${escapeHtml(sample)}">${escapeHtml(sample)}</div>`).join("")}
          ${cells.map((row) => `
            <div class="patho-heatmap-row-header" data-patho-sample-name="${escapeHtml(normalizePathoSampleName(row.sample))}" title="${escapeHtml(row.sample)}">${escapeHtml(row.sample)}</div>
            ${row.values.map((cell) => {
              const isDiagonal = cell.source === cell.target;
              const bg = isDiagonal
                ? "rgba(62, 84, 111, 0.08)"
                : buildHeatValueColor(cell.value, min, max, hue);
              const color = isDiagonal
                ? "var(--report-ink-soft)"
                : buildHeatTextColor(cell.value, min, max);
              const detail = `${cell.source} vs ${cell.target}: ${cell.display}${unit}`;
              return `
                <button
                  class="patho-heatmap-cell${isDiagonal ? " is-diagonal" : ""}"
                  type="button"
                  style="--cell-bg:${bg}; --cell-color:${color};"
                  data-heatmap-detail-text="${escapeHtml(detail)}"
                  data-patho-sample-a="${escapeHtml(normalizePathoSampleName(cell.source))}"
                  data-patho-sample-b="${escapeHtml(normalizePathoSampleName(cell.target))}"
                  title="${escapeHtml(detail)}"
                >${escapeHtml(cell.display)}</button>
              `;
            }).join("")}
          `).join("")}
        </div>
      </div>
    </div>
  `;
  const detailNode = container.querySelector("[data-heatmap-detail]");
  let lockedButton = null;
  container.querySelectorAll(".patho-heatmap-cell").forEach((button) => {
    button.addEventListener("mouseenter", () => {
      if (lockedButton) return;
      if (detailNode) detailNode.textContent = button.dataset.heatmapDetailText || defaultDetail;
    });
    button.addEventListener("focus", () => {
      if (lockedButton) return;
      if (detailNode) detailNode.textContent = button.dataset.heatmapDetailText || defaultDetail;
    });
    button.addEventListener("click", () => {
      if (lockedButton === button) {
        button.classList.remove("is-locked");
        lockedButton = null;
        if (detailNode) detailNode.textContent = defaultDetail;
        return;
      }
      if (lockedButton) lockedButton.classList.remove("is-locked");
      lockedButton = button;
      button.classList.add("is-locked");
      if (detailNode) detailNode.textContent = button.dataset.heatmapDetailText || defaultDetail;
    });
  });
  container.addEventListener("mouseleave", () => {
    if (!lockedButton && detailNode) detailNode.textContent = defaultDetail;
  });
}

function renderCompactPairList(containerId, title, section) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const columns = Array.isArray(section?.columns) ? section.columns : [];
  const rows = Array.isArray(section?.rows) ? section.rows : [];
  if (!columns.length || !rows.length) {
    container.innerHTML = `
      <div class="empty-table-state">
        <strong>${escapeHtml(title)} 暂无数据</strong>
        <p class="empty-copy">当前没有可展示的高相似样本对。</p>
      </div>
    `;
    return;
  }
  const topRows = rows.slice(0, 8);
  container.innerHTML = `
    <div class="patho-pair-list">
      ${topRows.map((row, index) => `
        <article class="patho-pair-item" data-patho-sample-group="${escapeHtml(encodePathoSampleGroup([row[0], row[1]]))}">
          <div class="patho-pair-rank">${index + 1}</div>
          <div class="patho-pair-main">
            <strong>${escapeHtml(String(row[0] || "-"))} <span>vs</span> ${escapeHtml(String(row[1] || "-"))}</strong>
            <p>ANI ${escapeHtml(String(row[2] || "-"))}% · Ref ${escapeHtml(String(row[3] || "-"))}% · Query ${escapeHtml(String(row[4] || "-"))}%</p>
          </div>
        </article>
      `).join("")}
    </div>
  `;
}

function inferPathoInterpretation(section) {
  const clusterCount = Number(section?.cluster?.summary?.cluster_count || 0);
  const largestCluster = Number(section?.cluster?.summary?.largest_cluster || 0);
  const singletonClusters = Number(section?.cluster?.summary?.singleton_clusters || 0);
  const snpSummary = section?.snp_matrix?.summary || {};
  const aniSummary = section?.ani?.summary || {};
  const mutationSummary = section?.mutations?.summary || {};
  const topAniRow = Array.isArray(section?.ani?.top_pairs?.rows) && section.ani.top_pairs.rows.length
    ? section.ani.top_pairs.rows[0]
    : null;

  const minSnp = Number(snpSummary?.min_distance ?? NaN);
  const medianSnp = Number(snpSummary?.median_distance ?? NaN);
  const maxAni = Number(aniSummary?.max_ani ?? NaN);
  const pairCount = Number(snpSummary?.pair_count || 0);
  const dominantMutationSample = String(mutationSummary?.max_mutation_sample || "").trim() || "--";
  const dominantMutationCount = Number(mutationSummary?.max_mutation_count || 0);

  let clusterLabel = "以散在样本为主";
  let clusterTone = "neutral";
  let clusterSummary = "当前没有明显的大簇结构，样本关系更接近散点分布。";
  if (largestCluster >= 5) {
    clusterLabel = "存在明确聚集簇";
    clusterTone = "high";
    clusterSummary = `当前共形成 ${clusterCount} 个簇，最大簇包含 ${largestCluster} 个样本，应优先围绕该簇做流调和时空核查。`;
  } else if (largestCluster >= 3) {
    clusterLabel = "存在中等规模近缘簇";
    clusterTone = "mid";
    clusterSummary = `当前形成 ${clusterCount} 个簇，最大簇为 ${largestCluster} 个样本，已具备进一步追踪价值。`;
  } else if (largestCluster === 2) {
    clusterLabel = "存在小范围近邻对";
    clusterTone = "mid";
    clusterSummary = `当前以双样本近邻对为主，暂未形成更大规模聚类。`;
  }

  let evidenceLabel = "传播关联证据偏弱";
  let evidenceTone = "neutral";
  let evidenceSummary = "建议以树拓扑和样本背景信息为主，谨慎解释传播关联。";
  if ((Number.isFinite(minSnp) && minSnp <= 10) || (Number.isFinite(maxAni) && maxAni >= 99.9)) {
    evidenceLabel = "传播关联证据较强";
    evidenceTone = "high";
    evidenceSummary = `最小 SNP 距离为 ${Number.isFinite(minSnp) ? minSnp : "--"}，最高 ANI 为 ${Number.isFinite(maxAni) ? maxAni.toFixed(2) : "--"}%，簇内近缘关系较强。`;
  } else if ((Number.isFinite(minSnp) && minSnp <= 25) || (Number.isFinite(maxAni) && maxAni >= 99.5)) {
    evidenceLabel = "传播关联证据中等";
    evidenceTone = "mid";
    evidenceSummary = `最小 SNP 距离为 ${Number.isFinite(minSnp) ? minSnp : "--"}，最高 ANI 为 ${Number.isFinite(maxAni) ? maxAni.toFixed(2) : "--"}%，建议结合流行病学背景判读。`;
  }

  const topPairLabel = topAniRow
    ? `${String(topAniRow[0] || "-")} / ${String(topAniRow[1] || "-")}`
    : "--";
  const topPairDetail = topAniRow
    ? `ANI ${String(topAniRow[2] || "-")}%`
    : "当前没有可用于优先追踪的高相似样本对。";
  const mutationDetail = dominantMutationSample !== "--" && dominantMutationCount > 0
    ? `${dominantMutationSample}（${dominantMutationCount} 个变异位点）`
    : "当前未形成明显的高突变样本。";

  const lead = clusterTone === "high"
    ? "先围绕主簇确认是否存在时空重叠，再用 SNP 与 ANI 双证据判断簇内传播可信度。"
    : clusterTone === "mid"
      ? "先确认近缘样本对和小规模簇的背景联系，再决定是否需要扩大追踪范围。"
      : "先排除散发样本中的偶然近邻关系，再决定是否进入深入传播调查。";

  return {
    lead,
    cards: [
      {
        title: "是否成簇",
        tone: clusterTone,
        badge: clusterLabel,
        summary: clusterSummary,
        metrics: [
          { label: "成簇数量", value: clusterCount || "--" },
          { label: "最大簇规模", value: largestCluster || "--" },
          { label: "单例簇", value: singletonClusters || "--" },
        ],
        link: "#section-patho-cluster",
        linkLabel: "查看成簇信息",
        highlightTargets: [
          "#section-patho-cluster",
          "#patho-cluster-summary",
          "#patho-cluster-table",
          "#section-patho-grapetree",
          "#patho-grapetree-card",
          "#section-patho-core-tree",
          "#patho-core-tree-card",
        ],
      },
      {
        title: "距离是否支持传播关联",
        tone: evidenceTone,
        badge: evidenceLabel,
        summary: evidenceSummary,
        metrics: [
          { label: "最小 SNP", value: Number.isFinite(minSnp) ? minSnp : "--" },
          { label: "中位 SNP", value: Number.isFinite(medianSnp) ? medianSnp : "--" },
          { label: "最高 ANI", value: Number.isFinite(maxAni) ? `${maxAni.toFixed(2)}%` : "--" },
        ],
        link: "#section-patho-distance",
        linkLabel: "查看距离证据",
        highlightTargets: [
          "#section-patho-distance",
          "#patho-distance-bins",
          "#patho-snp-summary",
          "#patho-snp-matrix-table",
          "#section-patho-ani",
          "#patho-ani-summary",
          "#patho-ani-top-pairs",
          "#patho-ani-preview-table",
        ],
      },
      {
        title: "下一步优先追踪",
        tone: "accent",
        badge: pairCount > 0 ? "优先看近邻样本对与高变异样本" : "等待更多比较样本",
        summary: `优先检查高相似样本对 ${topPairLabel}，并补看高突变样本 ${mutationDetail}。`,
        metrics: [
          { label: "高相似样本对", value: topPairLabel },
          { label: "相似性线索", value: topPairDetail },
          { label: "高突变样本", value: mutationDetail },
        ],
        link: "#section-patho-ani",
        linkLabel: "进入优先样本",
        highlightTargets: [
          "#section-patho-ani",
          "#patho-ani-top-pairs",
          "#patho-ani-preview-table",
          "#section-patho-mutation",
          "#patho-mutation-summary",
          "#patho-mutation-table",
          "#section-patho-grapetree",
          "#patho-grapetree-card",
        ],
      },
    ],
  };
}

function renderPathoInterpretationBand(containerId, section) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const interpretation = inferPathoInterpretation(section || {});
  container.innerHTML = `
    <article class="result-card patho-interpretation-card">
      <div class="patho-interpretation-topline">
        <span class="section-chip">Interpretation</span>
        <p>${escapeHtml(interpretation.lead)}</p>
      </div>
      <div class="patho-interpretation-grid">
        ${interpretation.cards.map((card) => `
          <article class="patho-interpretation-panel tone-${escapeHtml(card.tone)}" data-patho-link-card data-patho-primary-target="${escapeHtml(card.link)}" data-patho-highlight-targets="${escapeHtml(JSON.stringify(card.highlightTargets || []))}">
            <div class="patho-interpretation-head">
              <span>${escapeHtml(card.title)}</span>
              <strong>${escapeHtml(card.badge)}</strong>
            </div>
            <p class="patho-interpretation-summary">${escapeHtml(card.summary)}</p>
            <dl class="patho-interpretation-metrics">
              ${card.metrics.map((item) => `
                <div>
                  <dt>${escapeHtml(item.label)}</dt>
                  <dd>${escapeHtml(String(item.value))}</dd>
                </div>
              `).join("")}
            </dl>
            <a class="patho-interpretation-link" href="${escapeHtml(card.link)}">${escapeHtml(card.linkLabel)}</a>
          </article>
        `).join("")}
      </div>
    </article>
  `;
}

let activePathoHighlightTimer = null;

function clearPathoLinkedHighlights() {
  document.querySelectorAll(".is-patho-linked-focus").forEach((node) => node.classList.remove("is-patho-linked-focus"));
  document.querySelectorAll(".is-patho-link-active").forEach((node) => node.classList.remove("is-patho-link-active"));
}

function activatePathoLinkedHighlights(selectors = [], card = null) {
  clearPathoLinkedHighlights();
  if (card) {
    card.classList.add("is-patho-link-active");
  }
  selectors.forEach((selector) => {
    document.querySelectorAll(selector).forEach((node) => node.classList.add("is-patho-linked-focus"));
    if (selector.startsWith("#section-")) {
      const navLink = document.querySelector(`.report-nav-link[href="${selector}"]`);
      if (navLink) navLink.classList.add("is-patho-linked-focus");
    }
  });
  if (activePathoHighlightTimer) {
    window.clearTimeout(activePathoHighlightTimer);
  }
  activePathoHighlightTimer = window.setTimeout(() => {
    clearPathoLinkedHighlights();
    activePathoHighlightTimer = null;
  }, 2600);
}

function bindPathoInterpretationBand() {
  document.querySelectorAll("[data-patho-link-card]").forEach((card) => {
    if (card.dataset.pathoLinkBound === "true") return;
    card.dataset.pathoLinkBound = "true";
    const rawTargets = card.getAttribute("data-patho-highlight-targets") || "[]";
    let targets = [];
    try {
      targets = JSON.parse(rawTargets);
    } catch (_error) {
      targets = [];
    }
    const primaryTarget = card.getAttribute("data-patho-primary-target") || "";
    const focusTargets = Array.isArray(targets) ? targets.filter(Boolean) : [];
    const activate = () => activatePathoLinkedHighlights(focusTargets, card);
    card.addEventListener("mouseenter", activate);
    card.addEventListener("focusin", activate);
    card.addEventListener("click", (event) => {
      const targetNode = primaryTarget ? document.querySelector(primaryTarget) : null;
      activate();
      if (targetNode) {
        targetNode.scrollIntoView({ behavior: "smooth", block: "start" });
      }
      const clickedLink = event.target instanceof HTMLElement ? event.target.closest("a[href^='#']") : null;
      if (clickedLink) {
        event.preventDefault();
      }
    });
  });
}

function renderCommunityBetaPcoa(containerId, section) {
  return renderCommunityBetaOrdination(
    containerId,
    Array.isArray(section?.pcoa_points) ? section.pcoa_points : [],
    {
      title: "PCoA 交互散点图",
      emptyLabel: "PCoA",
      xLabel: "PCo1",
      yLabel: "PCo2",
      insight: `当前展示 ${(Array.isArray(section?.pcoa_points) ? section.pcoa_points.length : 0)} 个样本在 Bray-Curtis 距离下的前两主坐标轴位置；悬停点位可查看样本名和分组。`,
    },
  );
}

function renderCommunityBetaNmds(containerId, section) {
  return renderCommunityBetaOrdination(
    containerId,
    Array.isArray(section?.nmds_points) ? section.nmds_points : [],
    {
      title: "NMDS 交互散点图",
      emptyLabel: "NMDS",
      xLabel: "NMDS1",
      yLabel: "NMDS2",
      insight: `当前展示 ${(Array.isArray(section?.nmds_points) ? section.nmds_points.length : 0)} 个样本在 Bray-Curtis 距离下的 NMDS 前两轴位置；悬停点位可查看样本名和分组。`,
      stress: section?.summary?.nmds_stress,
    },
  );
}

function tokenizeNewick(text) {
  return String(text || "")
    .split(/\s*(;|\(|\)|,|:)\s*/)
    .filter((token) => token !== "");
}

function parseNewick(text) {
  const tokens = tokenizeNewick(text);
  let index = 0;
  function parseNode() {
    const token = tokens[index];
    if (token === "(") {
      index += 1;
      const children = [];
      while (index < tokens.length && tokens[index] !== ")") {
        if (tokens[index] === ",") {
          index += 1;
          continue;
        }
        children.push(parseNode());
      }
      if (tokens[index] === ")") {
        index += 1;
      }
      let name = "";
      let length = 0;
      if (index < tokens.length && !["(", ")", ",", ":", ";"].includes(tokens[index])) {
        name = tokens[index];
        index += 1;
      }
      if (tokens[index] === ":") {
        index += 1;
        length = Number.parseFloat(tokens[index] || "0");
        index += 1;
      }
      return { name, length: Number.isFinite(length) ? length : 0, children };
    }
    const name = token && ![",", ")", ";"].includes(token) ? token : "";
    if (name) {
      index += 1;
    }
    let length = 0;
    if (tokens[index] === ":") {
      index += 1;
      length = Number.parseFloat(tokens[index] || "0");
      index += 1;
    }
    return { name, length: Number.isFinite(length) ? length : 0, children: [] };
  }
  try {
    return parseNode();
  } catch (_error) {
    return null;
  }
}

function collectTreeLeaves(node, leaves = []) {
  if (!node) return leaves;
  if (!Array.isArray(node.children) || !node.children.length) {
    leaves.push(node);
    return leaves;
  }
  node.children.forEach((child) => collectTreeLeaves(child, leaves));
  return leaves;
}

function annotateTreeLayout(node, distance = 0, depth = 0, state = { leafIndex: 0 }) {
  node._distance = distance;
  node._depth = depth;
  if (!Array.isArray(node.children) || !node.children.length) {
    node._y = state.leafIndex;
    state.leafIndex += 1;
    return;
  }
  node.children.forEach((child) => annotateTreeLayout(child, distance + (Number(child.length) || 0), depth + 1, state));
  node._y = node.children.reduce((sum, child) => sum + child._y, 0) / node.children.length;
}

function buildTreeLines(node, scaleX, scaleY, lines, labels) {
  const x = scaleX(node._distance || 0);
  const y = scaleY(node._y || 0);
  const children = Array.isArray(node.children) ? node.children : [];
  if (!children.length) {
    labels.push({
      x: x + 8,
      y: y + 4,
      text: node.name || "未命名样本",
      internal: false,
      sample: normalizePathoSampleName(node.name || ""),
    });
    lines.push(`<circle cx="${x}" cy="${y}" r="2.5" class="patho-tree-node"></circle>`);
    return;
  }
  const childYs = [];
  children.forEach((child) => {
    const childX = scaleX(child._distance || 0);
    const childY = scaleY(child._y || 0);
    childYs.push(childY);
    lines.push(`<line x1="${x}" y1="${childY}" x2="${childX}" y2="${childY}" class="patho-tree-branch"></line>`);
    buildTreeLines(child, scaleX, scaleY, lines, labels);
  });
  if (childYs.length) {
    lines.push(`<line x1="${x}" y1="${Math.min(...childYs)}" x2="${x}" y2="${Math.max(...childYs)}" class="patho-tree-branch"></line>`);
  }
  if (node.name) {
    labels.push({ x: x + 6, y: y - 6, text: node.name, internal: true });
  }
}

function isTreeSampleFocus(candidateSample, focusSamples) {
  const candidate = normalizePathoSampleName(candidateSample).toLowerCase();
  const normalizedFocus = extractPathoSampleNames(focusSamples).map((item) => item.toLowerCase()).filter(Boolean);
  if (!candidate || !normalizedFocus.length) return false;
  return normalizedFocus.some((focus) => {
    if (!focus) return false;
    if (candidate === focus) return true;
    const delimitedMatch = ["_", "-", ".", "/", " "].some((delimiter) => (
      candidate.startsWith(`${focus}${delimiter}`)
      || candidate.endsWith(`${delimiter}${focus}`)
      || candidate.includes(`${delimiter}${focus}${delimiter}`)
    ));
    if (delimitedMatch) return true;
    return focus.length >= 8 && candidate.includes(focus);
  });
}

function renderNewickTreeCard(containerId, treeSection, options = {}) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (String(treeSection?.status || "") !== "ready" || !String(treeSection?.newick || "").trim()) {
    container.innerHTML = `
      <article class="result-card">
        <div class="card-head">
          <div class="card-title-stack">
            <span class="section-chip">Tree</span>
            <h3>${escapeHtml(treeSection?.label || "树文件")}</h3>
          </div>
          <span class="card-tag">${escapeHtml(treeSection?.file_name || "-")}</span>
        </div>
        <div class="empty-box">
          <p>当前文件为空或未生成，无法渲染树图。</p>
        </div>
      </article>
    `;
    return;
  }
  const tree = parseNewick(treeSection.newick);
  if (!tree) {
    container.innerHTML = `
      <article class="result-card">
        <div class="card-head">
          <div class="card-title-stack">
            <span class="section-chip">Tree</span>
            <h3>${escapeHtml(treeSection?.label || "树文件")}</h3>
          </div>
          <span class="card-tag">${escapeHtml(treeSection?.file_name || "-")}</span>
        </div>
        <div class="empty-box">
          <p>Newick 解析失败，当前仅保留原始文本。</p>
          <pre class="patho-tree-raw">${escapeHtml(String(treeSection.newick || ""))}</pre>
        </div>
      </article>
    `;
    return;
  }
  annotateTreeLayout(tree);
  const leaves = collectTreeLeaves(tree, []);
  const rowHeight = options.rowHeight || 18;
  const innerHeight = Math.max(180, leaves.length * rowHeight + 28);
  const labelsWidth = options.labelsWidth || 180;
  const drawingWidth = options.drawingWidth || 1180;
  const plotWidth = drawingWidth - labelsWidth - 28;
  const distances = leaves.map((leaf) => Number(leaf._distance) || 0);
  const depths = leaves.map((leaf) => Number(leaf._depth) || 0);
  const maxDistance = Math.max(0, ...distances);
  const maxDepth = Math.max(1, ...depths);
  const scaleX = (value) => 12 + ((maxDistance > 0 ? value / maxDistance : value / maxDepth) * plotWidth);
  const scaleY = (value) => 16 + value * rowHeight;
  const lines = [];
  const labels = [];
  buildTreeLines(tree, scaleX, scaleY, lines, labels);
  const focusSampleNames = extractPathoSampleNames(
    options.focusSampleNames && (Array.isArray(options.focusSampleNames) ? options.focusSampleNames : [options.focusSampleNames]).length
      ? options.focusSampleNames
      : [
          currentReportData?.task?.sample_display_name || "",
          currentReportData?.task?.sample_name || "",
          currentReportData?.task?.name || "",
        ],
  );
  container.innerHTML = `
    <article class="result-card patho-tree-card">
      <div class="card-head">
        <div class="card-title-stack">
          <span class="section-chip">Tree</span>
          <h3>${escapeHtml(treeSection.label || "树图")}</h3>
        </div>
        <span class="card-tag">${escapeHtml(treeSection.file_name || "-")}</span>
      </div>
      <div class="patho-tree-meta">
        <span>叶节点：<strong>${escapeHtml(String(treeSection.leaf_count || leaves.length))}</strong></span>
        <span>字符数：<strong>${escapeHtml(String(treeSection.char_count || String(treeSection.newick || "").length))}</strong></span>
      </div>
      <div class="patho-tree-frame">
        <svg class="patho-tree-svg" viewBox="0 0 ${drawingWidth} ${innerHeight}" preserveAspectRatio="xMinYMin meet" role="img" aria-label="${escapeHtml(treeSection.label || 'Tree')}">
          ${lines.join("")}
          ${labels.map((label) => `
            <text x="${label.x}" y="${label.y}" class="${label.internal ? "patho-tree-label patho-tree-label-internal" : `patho-tree-label${isTreeSampleFocus(label.sample || label.text || "", focusSampleNames) ? " is-patho-sample-focus" : ""}`}"${label.internal ? "" : ` data-patho-sample-name="${escapeHtml(label.sample || "")}"`}>${escapeHtml(label.text)}</text>
          `).join("")}
        </svg>
      </div>
      <details class="patho-tree-details">
        <summary>查看原始 Newick</summary>
        <pre class="patho-tree-raw">${escapeHtml(String(treeSection.newick || ""))}</pre>
      </details>
    </article>
  `;
}

function buildGrapeTreeMetadataTsv(mlstSection) {
  const columns = Array.isArray(mlstSection?.columns) ? mlstSection.columns.map((value) => String(value || "").trim()) : [];
  const rows = Array.isArray(mlstSection?.rows) ? mlstSection.rows : [];
  if (columns.length < 3 || !rows.length) return "";
  const metadataHeaders = ["ID", ...columns];
  const lines = [metadataHeaders.join("\t")];
  rows.forEach((row) => {
    const cells = Array.isArray(row) ? row.map((value) => String(value ?? "").trim()) : [];
    const sampleName = cells[0] || "";
    if (!sampleName) return;
    lines.push([sampleName, ...cells].join("\t"));
  });
  return lines.join("\n");
}

function buildITOLMetadataTsv(mlstSection) {
  return buildGrapeTreeMetadataTsv(mlstSection);
}

function enhancePathoClusterTable(section) {
  const rows = Array.isArray(section?.cluster?.table?.rows) ? section.cluster.table.rows : [];
  const tableRows = document.querySelectorAll("#patho-cluster-table tbody tr");
  tableRows.forEach((rowNode, index) => {
    const row = rows[index] || [];
    const sampleGroup = encodePathoSampleGroup(row[2]);
    if (!sampleGroup || sampleGroup === "[]") return;
    rowNode.setAttribute("data-patho-sample-group", sampleGroup);
    const sampleCell = rowNode.children[2];
    if (sampleCell instanceof HTMLElement) {
      sampleCell.setAttribute("data-patho-sample-group", sampleGroup);
    }
  });
}

function enhancePathoMutationTable(section) {
  const rows = Array.isArray(section?.mutations?.table?.rows) ? section.mutations.table.rows : [];
  const tableRows = document.querySelectorAll("#patho-mutation-table tbody tr");
  tableRows.forEach((rowNode, index) => {
    const row = rows[index] || [];
    const sampleGroup = encodePathoSampleGroup(row[5]);
    if (!sampleGroup || sampleGroup === "[]") return;
    rowNode.setAttribute("data-patho-sample-group", sampleGroup);
    const sampleCell = rowNode.children[5];
    if (sampleCell instanceof HTMLElement) {
      sampleCell.setAttribute("data-patho-sample-group", sampleGroup);
    }
  });
}

function renderOfficialGrapeTreeCard(containerId, treeSection, options = {}) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (String(treeSection?.status || "") !== "ready" || !String(treeSection?.newick || "").trim()) {
    renderNewickTreeCard(containerId, treeSection, { rowHeight: 18, drawingWidth: 1280, labelsWidth: 220 });
    return;
  }
  const iframeId = `grapetree-frame-${Math.random().toString(36).slice(2, 10)}`;
  const taskId = String(treeSection?.task_id || treeSection?.report_id || treeSection?.file_name || "default");
  const fileName = String(treeSection.file_name || "grapetree.nwk");
  const storageKey = `task:${taskId}:grapetree:${fileName}`;
  const sourceUrl = `/public/GrapeTree-master/MSTree_holder.html?portalStorageKey=${encodeURIComponent(storageKey)}`;
  container.innerHTML = `
    <article class="result-card patho-tree-card patho-grapetree-card">
      <div class="card-head">
        <div class="card-title-stack">
          <span class="section-chip">Official Viewer</span>
          <h3>${escapeHtml(treeSection.label || "GrapeTree 最小生成树")}</h3>
        </div>
        <span class="card-tag">${escapeHtml(treeSection.file_name || "-")}</span>
      </div>
      <div class="patho-tree-meta">
        <span>叶节点：<strong>${escapeHtml(String(treeSection.leaf_count || "--"))}</strong></span>
        <span>载入方式：<strong>public/GrapeTree-master</strong></span>
        <a class="patho-grapetree-link" href="${sourceUrl}" target="_blank" rel="noopener noreferrer">单独打开 Grapetree</a>
      </div>
      <div class="patho-grapetree-frame-wrap">
        <iframe id="${iframeId}" class="patho-grapetree-frame" src="${sourceUrl}" title="GrapeTree Official Viewer"></iframe>
      </div>
      <details class="patho-tree-details">
        <summary>查看原始 Newick</summary>
        <pre class="patho-tree-raw">${escapeHtml(String(treeSection.newick || ""))}</pre>
      </details>
    </article>
  `;
  const iframe = document.getElementById(iframeId);
  if (!iframe) return;
  const treeText = String(treeSection.newick || "");
  const metadataText = String(options?.metadataText || "");
  const metadataCategory = String(options?.metadataCategory || "");
  const tryInject = () => {
    try {
      const frameWindow = iframe.contentWindow;
      if (!frameWindow || typeof frameWindow.distributeFile !== "function") {
        return false;
      }
      if (typeof frameWindow.setPortalPersistenceKey === "function") {
        frameWindow.setPortalPersistenceKey(storageKey);
      }
      if (typeof frameWindow.setPortalOriginalTree === "function") {
        frameWindow.setPortalOriginalTree(treeText, fileName);
      }
      if (metadataText && typeof frameWindow.setPortalInjectedMetadata === "function") {
        frameWindow.setPortalInjectedMetadata(metadataText, metadataCategory || "ST");
      }
      if (typeof frameWindow.restorePortalTreeState === "function" && frameWindow.restorePortalTreeState()) {
        return true;
      }
      frameWindow.distributeFile(treeText, fileName);
      return true;
    } catch (_error) {
      return false;
    }
  };
  iframe.addEventListener("load", () => {
    if (tryInject()) return;
    let attempt = 0;
    const timer = window.setInterval(() => {
      attempt += 1;
      if (tryInject() || attempt >= 40) {
        window.clearInterval(timer);
      }
    }, 250);
  }, { once: true });
}

function renderITOLTreeCard(containerId, treeSection, options = {}) {
  const container = document.getElementById(containerId);
  if (!container) return;
  if (String(treeSection?.status || "") !== "ready" || !String(treeSection?.newick || "").trim()) {
    renderNewickTreeCard(containerId, treeSection, { rowHeight: 18, drawingWidth: 1280, labelsWidth: 220 });
    return;
  }
  const iframeId = `itol-frame-${Math.random().toString(36).slice(2, 10)}`;
  const taskId = String(treeSection?.task_id || treeSection?.report_id || treeSection?.file_name || "default");
  const fileName = String(treeSection.file_name || "tree.nwk");
  const storageKey = `task:${taskId}:itol:${fileName}`;
  const viewerVersion = "20260331-5";
  const treeData = treeSection?.itol_tree_json || null;
  try {
    if (treeData) {
      window.localStorage.setItem(`${storageKey}:treeData`, JSON.stringify(treeData));
    } else {
      window.localStorage.removeItem(`${storageKey}:treeData`);
    }
  } catch (_error) {
    // ignore storage errors and let iframe fallback to empty tree payload
  }
  const sourceUrl = `/public/itol/index.html?portalStorageKey=${encodeURIComponent(storageKey)}&v=${encodeURIComponent(viewerVersion)}`;
  container.innerHTML = `
    <article class="result-card patho-tree-card patho-grapetree-card">
      <div class="card-head">
        <div class="card-title-stack">
          <span class="section-chip">iTOL Viewer</span>
          <h3>${escapeHtml(treeSection.label || "核心系统发育树")}</h3>
        </div>
        <span class="card-tag">${escapeHtml(fileName)}</span>
      </div>
      <div class="patho-tree-meta">
        <span>叶节点：<strong>${escapeHtml(String(treeSection.leaf_count || "--"))}</strong></span>
        <span>载入方式：<strong>public/itol</strong></span>
        <a class="patho-grapetree-link" href="${sourceUrl}" target="_blank" rel="noopener noreferrer">单独打开 iTOL</a>
      </div>
      <div class="patho-grapetree-frame-wrap">
        <iframe id="${iframeId}" class="patho-grapetree-frame" src="${sourceUrl}" title="iTOL Viewer"></iframe>
      </div>
      <details class="patho-tree-details">
        <summary>查看原始 Newick</summary>
        <pre class="patho-tree-raw">${escapeHtml(String(treeSection.newick || ""))}</pre>
      </details>
    </article>
  `;
  const iframe = document.getElementById(iframeId);
  if (!iframe) return;
  const tryInject = () => {
    try {
      const frameWindow = iframe.contentWindow;
      if (frameWindow && typeof frameWindow.postMessage === "function") {
        frameWindow.postMessage({
          type: "portal-itol-tree",
          treeData,
        }, "*");
      }
      if (!frameWindow || typeof frameWindow.setViewerPayload !== "function") {
        return false;
      }
      frameWindow.setViewerPayload({
        storageKey,
        treeData,
      });
      return true;
    } catch (_error) {
      return false;
    }
  };
  iframe.addEventListener("load", () => {
    tryInject();
    let attempt = 0;
    const timer = window.setInterval(() => {
      attempt += 1;
      tryInject();
      if (attempt >= 16) {
        window.clearInterval(timer);
      }
    }, 400);
  }, { once: true });
}

function setPathoSourceReportChrome(task) {
  const shell = document.querySelector(".report-shell");
  if (shell) {
    shell.dataset.reportKind = "pathosource";
  }
  document.querySelectorAll(".report-scene-switcher").forEach((node) => node.classList.add("hidden"));
  const clinicalNav = document.getElementById("nav-group-clinical-wrapper");
  if (clinicalNav) clinicalNav.classList.add("hidden");
  const title = document.querySelector(".report-title-block h1");
  if (title) {
    title.textContent = "分子溯源分析结果";
  }
  const subtitle = document.querySelector(".report-subtitle");
  if (subtitle) {
    subtitle.textContent = "本报告围绕成簇关系、SNP 距离、ANI、MLST 及多树拓扑组织结果，面向实验室溯源判读与聚类调查。";
  }
}

function renderPathoSourceReport(data) {
  const task = data?.task || {};
  const section = data?.sections?.pathosource || {};
  setPathoSourceReportChrome(task);
  const nav = document.querySelector(".report-nav");
  if (nav) {
    nav.innerHTML = buildPathoSourceNav();
  }
  const content = document.querySelector(".report-content");
  if (content) {
    content.innerHTML = buildPathoSourceLayout();
  }
  fillTaskMeta(task);
  renderPathoInterpretationBand("patho-interpretation-band", section || {});
  document.getElementById("overview-metrics").innerHTML = buildMetricCards(data.overview_metrics || []);
  renderMiniStatGrid("patho-cluster-summary", [
    { label: "成簇数量", value: String(section?.cluster?.summary?.cluster_count ?? "--") },
    { label: "最大簇规模", value: String(section?.cluster?.summary?.largest_cluster ?? "--") },
    { label: "单例簇数量", value: String(section?.cluster?.summary?.singleton_clusters ?? "--") },
  ]);
  buildTableCard("patho-cluster-table", "成簇信息整理", section?.cluster?.table?.columns || [], section?.cluster?.table?.rows || []);
  renderPathoSourceDistanceBins("patho-distance-bins", section?.distance_bins || {});
  renderMiniStatGrid("patho-snp-summary", [
    { label: "样本数量", value: String(section?.snp_matrix?.summary?.sample_count ?? "--") },
    { label: "比较对数", value: String(section?.snp_matrix?.summary?.pair_count ?? "--") },
    { label: "最小 SNP", value: String(section?.snp_matrix?.summary?.min_distance ?? "--") },
    { label: "最大 SNP", value: String(section?.snp_matrix?.summary?.max_distance ?? "--") },
    { label: "中位 SNP", value: String(section?.snp_matrix?.summary?.median_distance ?? "--") },
    { label: "零距离对", value: String(section?.snp_matrix?.summary?.zero_distance_pairs ?? "--") },
  ]);
  renderHeatmapMatrix("patho-snp-matrix-table", "SNP 距离矩阵热图", section?.snp_matrix?.preview || {}, { hue: 208, unit: " SNP" });
  renderMiniStatGrid("patho-ani-summary", [
    { label: "样本数量", value: String(section?.ani?.summary?.sample_count ?? "--") },
    { label: "比较对数", value: String(section?.ani?.summary?.pair_count ?? "--") },
    { label: "最低 ANI", value: section?.ani?.summary?.min_ani == null ? "--" : `${section.ani.summary.min_ani}%` },
    { label: "最高 ANI", value: section?.ani?.summary?.max_ani == null ? "--" : `${section.ani.summary.max_ani}%` },
    { label: "中位 ANI", value: section?.ani?.summary?.median_ani == null ? "--" : `${section.ani.summary.median_ani}%` },
  ]);
  renderCompactPairList("patho-ani-top-pairs", "ANI 最高样本对", section?.ani?.top_pairs || {});
  renderHeatmapMatrix("patho-ani-preview-table", "ANI 矩阵热图", section?.ani?.preview || {}, { hue: 32, unit: "%" });
  renderOfficialGrapeTreeCard(
    "patho-grapetree-card",
    {
      ...(section?.trees?.grapetree || {}),
      task_id: task?.id || "",
    },
    {},
  );
  renderMiniStatGrid("patho-mlst-summary", [
    { label: "样本数量", value: String(section?.mlst?.summary?.sample_count ?? "--") },
    { label: "MLST 方案", value: String(section?.mlst?.summary?.scheme_count ?? "--") },
    { label: "ST 数量", value: String(section?.mlst?.summary?.st_count ?? "--") },
    { label: "优势 ST", value: String(section?.mlst?.summary?.dominant_st ?? "--") },
  ]);
  renderNewickTreeCard("patho-mlst-tree-card", section?.trees?.mlst || {}, { rowHeight: 20, drawingWidth: 1120, labelsWidth: 180 });
  const pathoMlstContainer = document.getElementById("patho-mlst-table");
  if (pathoMlstContainer) {
    pathoMlstContainer.dataset.exportTitle = "PathoSource_MLST统计结果";
    renderInteractiveContigTable(
      pathoMlstContainer,
      section?.mlst?.columns || [],
      section?.mlst?.rows || [],
      "patho-mlst-table",
    );
  }
  renderMiniStatGrid("patho-mutation-summary", [
    { label: "总位点数", value: String(section?.mutations?.summary?.site_count ?? "--") },
    { label: "发生变异位点", value: String(section?.mutations?.summary?.mutated_site_count ?? "--") },
    { label: "出现变异样本", value: String(section?.mutations?.summary?.mutated_sample_count ?? "--") },
    { label: "最高突变样本", value: String(section?.mutations?.summary?.max_mutation_sample ?? "--") },
    { label: "最高突变数", value: String(section?.mutations?.summary?.max_mutation_count ?? "--") },
  ]);
  const pathoMutationContainer = document.getElementById("patho-mutation-table");
  if (pathoMutationContainer) {
    pathoMutationContainer.dataset.exportTitle = "PathoSource_突变位点统计";
    renderInteractiveContigTable(
      pathoMutationContainer,
      section?.mutations?.table?.columns || [],
      section?.mutations?.table?.rows || [],
      "patho-mutate-table",
    );
  }
  renderITOLTreeCard(
    "patho-core-tree-card",
    {
      ...(section?.trees?.core || {}),
      task_id: task?.id || "",
    },
    {
      title: "核心 SNP 系统发育树",
    },
  );
  enhancePathoClusterTable(section || {});
  enhancePathoMutationTable(section || {});
  bindPathoInterpretationBand();
  bindPathoSampleInteractions();
  initializeReportNav();
  bindChartExportButtons();
}

function renderCommunityReport(data) {
  const task = data?.task || {};
  const section = data?.sections?.community || {};
  const outputRows = Array.isArray(section?.outputs?.rows) ? section.outputs.rows : [];
  const normalizedOutputs = outputRows.map((row) => ({
    label: String(row?.[0] || "").trim(),
    status: String(row?.[1] || "").trim().toLowerCase(),
    statusLabel: String(row?.[1] || "").trim() || "未知",
    path: String(row?.[2] || "").trim(),
  }));
  setCommunityReportChrome(task);
  const nav = document.querySelector(".report-nav");
  if (nav) {
    nav.innerHTML = buildCommunityNav();
  }
  const content = document.querySelector(".report-content");
  if (content) {
    content.innerHTML = buildCommunityLayout();
  }
  fillTaskMeta(task);
  document.getElementById("overview-metrics").innerHTML = buildMetricCards(data.overview_metrics || []);
  renderMiniStatGrid("community-summary-grid", [
    { label: "样本 ID 列", value: String(section?.summary?.sample_id_column ?? "--") },
    { label: "分组列", value: String(section?.summary?.group_column ?? "--") },
    { label: "流程模式", value: String(section?.summary?.workflow_mode ?? "--") },
    { label: "demux 样本数", value: String(section?.summary?.demux_sample_count ?? "--") },
    { label: "统计层级", value: String(section?.summary?.taxonomy_level ?? "--") },
    { label: "标准化方式", value: String(section?.summary?.normalization ?? "--") },
  ]);
  initializeReportNav();
  bindChartExportButtons();
  scheduleCommunitySectionRender("section-community-qc", () => {
  renderMiniStatGrid("community-demux-grid", [
    { label: "demux 样本数", value: String(section?.summary?.demux_sample_count ?? "--") },
    { label: "总 forward reads", value: String(section?.demux_preview?.summary?.forward_total ?? "--") },
    { label: "总 reverse reads", value: String(section?.demux_preview?.summary?.reverse_total ?? "--") },
    { label: "建议 trunc-len-f", value: String(section?.summary?.demux_trunc_len_f ?? "--") },
    { label: "建议 trunc-len-r", value: String(section?.summary?.demux_trunc_len_r ?? "--") },
    { label: "建议 depth", value: String(section?.summary?.demux_sampling_depth ?? "--") },
  ]);
  renderMiniStatGrid("community-denoise-grid", [
    { label: "去噪样本数", value: String(section?.denoise_preview?.summary?.sample_count ?? "--") },
    { label: "平均过滤保留率", value: section?.denoise_preview?.summary?.avg_pass_filter_pct != null ? `${section.denoise_preview.summary.avg_pass_filter_pct}%` : "--" },
    { label: "平均非嵌合保留率", value: section?.denoise_preview?.summary?.avg_non_chimeric_pct != null ? `${section.denoise_preview.summary.avg_non_chimeric_pct}%` : "--" },
    { label: "merged 总数", value: String(section?.denoise_preview?.summary?.merged_total ?? "--") },
  ]);
  renderCommunityQcRarefaction("community-qc-rarefaction-chart", section?.alpha || {});
  buildTableCard("community-demux-table", "demux.qzv 样本读数预览", section?.demux_preview?.columns || [], section?.demux_preview?.rows || []);
  buildTableCard("community-denoise-table", "denoising-stats-dada2.qzv 去噪结果预览", section?.denoise_preview?.columns || [], section?.denoise_preview?.rows || []);
  renderCommunityQcAssetCards(
    "community-qc-assets",
    task.id,
    normalizedOutputs.filter((item) => ["demux.qzv", "denoising-stats-dada2.qzv"].includes(item.label) && item.status === "ready"),
    section,
  );
  });
  scheduleCommunitySectionRender("section-community-composition", () => {
  renderMiniStatGrid("community-taxonomy-grid", [
    { label: "分类层级数", value: String(section?.taxa_abundance?.summary?.level_count ?? "--") },
    { label: "汇总条目数", value: String(section?.taxa_abundance?.summary?.row_count ?? "--") },
    { label: "统计层级", value: String(section?.summary?.taxonomy_level ?? "--") },
  ]);
  renderCommunityAssetLinks(
    "community-composition-assets",
    task.id,
    normalizedOutputs.filter((item) => ["taxonomy.qzv", "taxa-barplot.qzv", "taxonomy.tsv"].includes(item.label) && item.status === "ready"),
    "当前没有可展示的物种组成输出",
    "后续生成 taxonomy 与 taxa barplot 结果后，这里会显示对应文件入口。",
  );
  bindCommunityRankTabs(section);
  });
  scheduleCommunitySectionRender("section-community-alpha", () => {
  renderMiniStatGrid("community-alpha-grid", [
    { label: "稀释深度", value: String(section?.alpha?.summary?.selected_depth ?? section?.summary?.demux_sampling_depth ?? "--") },
    { label: "Shannon p 值", value: section?.alpha?.summary?.shannon_pvalue != null ? String(section.alpha.summary.shannon_pvalue) : "--" },
    { label: "特征数 p 值", value: section?.alpha?.summary?.observed_features_pvalue != null ? String(section.alpha.summary.observed_features_pvalue) : "--" },
    { label: "分组数", value: String(section?.alpha?.summary?.group_count ?? "--") },
  ]);
  renderCommunityAlphaWorkspace(section?.alpha || {});
  buildTableCard("community-alpha-table", "Alpha 样本预览", section?.alpha?.sample_columns || [], section?.alpha?.sample_rows || []);
  bindChartExportButtons();
  });
  scheduleCommunitySectionRender("section-community-beta", () => {
  renderMiniStatGrid("community-beta-grid", [
    { label: "距离指标", value: String(section?.beta?.summary?.measure ?? "--") },
    { label: "PERMANOVA R2", value: section?.beta?.summary?.permanova_r2 != null ? String(section.beta.summary.permanova_r2) : "--" },
    { label: "PERMANOVA p 值", value: section?.beta?.summary?.permanova_p != null ? String(section.beta.summary.permanova_p) : "--" },
    { label: "ANOSIM R", value: section?.beta?.summary?.anosim_r != null ? String(section.beta.summary.anosim_r) : "--" },
    { label: "ANOSIM p 值", value: section?.beta?.summary?.anosim_p != null ? String(section.beta.summary.anosim_p) : "--" },
    { label: "NMDS stress", value: section?.beta?.summary?.nmds_stress != null ? String(section.beta.summary.nmds_stress) : "--" },
    { label: "Betadisper p 值", value: section?.beta?.summary?.betadisper_p != null ? String(section.beta.summary.betadisper_p) : "--" },
  ]);
  buildTableCard(
    "community-beta-stats-table",
    "Beta 显著性统计",
    ["统计项", "主要指标", "说明"],
    [
      [
        "PERMANOVA",
        `R2=${section?.beta?.summary?.permanova_r2 ?? "--"}；p=${section?.beta?.summary?.permanova_p ?? "--"}`,
        "用于检验不同分组之间的整体群落组成差异。",
      ],
      [
        "ANOSIM",
        `R=${section?.beta?.summary?.anosim_r ?? "--"}；p=${section?.beta?.summary?.anosim_p ?? "--"}`,
        "用于评价组间排序分离程度，R 越大说明组间分离越明显。",
      ],
      [
        "Betadisper",
        `F=${section?.beta?.summary?.betadisper_f ?? "--"}；p=${section?.beta?.summary?.betadisper_p ?? "--"}`,
        "用于检验各分组内部离散度是否一致，便于辅助解释 PERMANOVA。",
      ],
    ],
  );
  buildTableCard(
    "community-beta-distance-table",
    "组内距离与分组统计",
    section?.beta?.group_distances?.columns?.length ? section.beta.group_distances.columns : (section?.beta?.group_counts?.columns || []),
    section?.beta?.group_distances?.rows?.length ? section.beta.group_distances.rows : (section?.beta?.group_counts?.rows || []),
  );
  renderCommunityBetaPcoa("community-beta-pcoa-figure", section?.beta || {});
  renderCommunityBetaNmds("community-beta-nmds-figure", section?.beta || {});
  renderCommunityBetaHeatmap("community-beta-distance-figure", section?.beta || {});
  renderCommunityBetaClusterComposition("community-beta-composition-figure", section?.beta || {}, section?.taxa_abundance || {});
  bindChartExportButtons();
  });
  scheduleCommunitySectionRender("section-community-biomarker", () => {
  renderMiniStatGrid("community-biomarker-summary-grid", [
    { label: "LEfSe 条目", value: String(section?.differential?.summary?.lefse_feature_count ?? "--") },
    { label: "LEfSe 显著项", value: String(section?.differential?.summary?.lefse_significant_count ?? "--") },
    { label: "最高 LDA", value: section?.differential?.summary?.lefse_lda_max != null ? String(section.differential.summary.lefse_lda_max) : "--" },
    { label: "RF 条目", value: String(section?.differential?.summary?.rf_feature_count ?? "--") },
    { label: "最高重要性", value: section?.differential?.summary?.rf_top_importance != null ? String(section.differential.summary.rf_top_importance) : "--" },
    { label: "当前预览", value: String(section?.differential?.summary?.preview_mode ?? "--") },
  ]);
  renderCommunityBiomarkerBars("community-biomarker-lefse-chart", section?.differential?.lefse?.rows || [], {
    kind: "lefse",
    label: "LEfSe LDA 排序",
    xLabel: "LDA score",
  });
  renderCommunityBiomarkerBars("community-biomarker-rf-chart", section?.differential?.rf?.rows || [], {
    kind: "rf",
    label: "Random Forest 特征重要性",
    xLabel: "Feature importance",
  });
  buildTableCard("community-differential-table", "LEfSe 结果预览", section?.differential?.lefse?.columns || section?.differential?.columns || [], section?.differential?.lefse?.rows || []);
  buildTableCard("community-rf-table", "RF 结果预览", section?.differential?.rf?.columns || [], section?.differential?.rf?.rows || []);
  bindChartExportButtons();
  });
  scheduleCommunitySectionRender("section-community-network", () => {
  renderMiniStatGrid("community-network-summary-grid", [
    { label: "节点数", value: String(section?.network?.summary?.node_count ?? "--") },
    { label: "边数", value: String(section?.network?.summary?.edge_count ?? "--") },
    { label: "模块数", value: String(section?.network?.summary?.module_count ?? "--") },
    { label: "平均 Degree", value: section?.network?.summary?.avg_degree != null ? String(section.network.summary.avg_degree) : "--" },
    { label: "网络密度", value: section?.network?.summary?.density != null ? String(section.network.summary.density) : "--" },
    { label: "模块度", value: section?.network?.summary?.modularity != null ? String(section.network.summary.modularity) : "--" },
  ]);
  renderCommunityNetworkGraph("community-network-graph", section?.network || {});
  renderCommunityNetworkHubBars("community-network-hubs", section?.network?.nodes || []);
  renderCommunityNetworkRoleScatter("community-network-role-scatter", section?.network?.nodes || []);
  renderCommunityNetworkLayeredRoles("community-network-role-layered", section?.network?.nodes || []);
  buildTableCard("community-network-module-table", "模块统计", section?.network?.module_preview?.columns || [], section?.network?.module_preview?.rows || []);
  buildTableCard("community-network-node-table", "关键节点预览", section?.network?.node_preview?.columns || [], section?.network?.node_preview?.rows || []);
  buildTableCard("community-network-edge-table", "关键边预览", section?.network?.edge_preview?.columns || [], section?.network?.edge_preview?.rows || []);
  bindChartExportButtons();
  });
  scheduleCommunitySectionRender("section-community-metadata", () => {
  buildTableCard("community-metadata-table", "元数据概况", section?.metadata?.columns || [], section?.metadata?.rows || []);
  });
  scheduleCommunitySectionRender("section-community-notes", () => {
  buildTableCard("community-modules-table", "模块规划", section?.modules?.columns || [], section?.modules?.rows || []);
  buildTableCard("community-outputs-table", "输出产物", section?.outputs?.columns || [], section?.outputs?.rows || []);
  renderCommunityNotes("community-notes-card", section?.notes || []);
  });
  const qcNote = document.getElementById("community-qc-note");
  if (qcNote) {
    const passPct = section?.denoise_preview?.summary?.avg_pass_filter_pct;
    const nonChimericPct = section?.denoise_preview?.summary?.avg_non_chimeric_pct;
    qcNote.innerHTML = `
      <p>当前质控区已直接读取 <code>demux.qzv</code> 和 <code>denoising-stats-dada2.qzv</code> 内部表格，而不只是提供文件入口。</p>
      <p>本次数据平均过滤保留率约为 <strong>${escapeHtml(passPct != null ? `${passPct}%` : "--")}</strong>，平均非嵌合保留率约为 <strong>${escapeHtml(nonChimericPct != null ? `${nonChimericPct}%` : "--")}</strong>。</p>
    `;
  }
}

function scheduleCommunitySectionRender(sectionId, renderFn) {
  if (typeof renderFn !== "function") return;
  const section = document.getElementById(sectionId);
  let executed = false;
  const run = () => {
    if (executed) return;
    executed = true;
    window.requestAnimationFrame(() => {
      window.setTimeout(() => {
        try {
          renderFn();
        } catch (error) {
          console.error(error);
        }
      }, 0);
    });
  };
  if (!section) {
    run();
    return;
  }
  if (!("IntersectionObserver" in window)) {
    window.setTimeout(run, 0);
    return;
  }
  const observer = new IntersectionObserver((entries) => {
    if (!entries.some((entry) => entry.isIntersecting)) return;
    observer.disconnect();
    run();
  }, {
    rootMargin: "320px 0px",
  });
  observer.observe(section);
}

async function loadReport() {
  const shell = document.querySelector(".report-shell");
  if (!shell) return;
  let data = window.__EMBEDDED_REPORT_DATA__ || currentReportData;
  if (!data) {
    const endpoint = new URL(shell.dataset.reportEndpoint, window.location.origin);
    const currentSample = new URLSearchParams(window.location.search).get("sample");
    if (currentSample) {
      endpoint.searchParams.set("sample", currentSample);
    }
    const response = await fetch(endpoint.toString(), { credentials: "same-origin" });
    data = await response.json();
    if (!response.ok) throw new Error(data.error || "结果数据加载失败");
  }
  currentReportData = data;
  let reportKind = String(data?.task?.report_kind || "").trim() || "default";
  if (isSarsCov2NextcladeReport(data)) reportKind = "sars-cov-2";
  else if (isHmpvNextcladeReport(data)) reportKind = "hmpv";
  else if (isDenvNextcladeReport(data)) reportKind = "denv";
  else if (isZikavNextcladeReport(data)) reportKind = "zikav";
  else if (isChikvNextcladeReport(data)) reportKind = "chikv";
  else if (isEbolaNextcladeReport(data)) reportKind = "ebola";
  else if (isHpivTypingReport(data)) reportKind = "hpiv";
  else if (isHadvTypingReport(data)) reportKind = "hadv";
  else if (isNorovirusTypingReport(data)) reportKind = "norovirus";
  else if (isEnterovirusTypingReport(data)) reportKind = "enterovirus";
  else if (isHepatovirusTypingReport(data)) reportKind = "hepatovirus";
  else if (isHivTypingReport(data)) reportKind = "hiv";
  else if (isBandavirusTypingReport(data)) reportKind = "bandavirus";
  else if (isOrthohantavirusTypingReport(data)) reportKind = "orthohantavirus";
  else if (isOrthoebolavirusTypingReport(data)) reportKind = "orthoebolavirus";
  else if (isAstroviridaeTypingReport(data)) reportKind = "astroviridae";
  else if (isRhinovirusTypingReport(data)) reportKind = "rhinovirus";
  else if (isSeasonalHcovTypingReport(data)) reportKind = "seasonal_hcov";
  else if (isRotavirusTypingReport(data)) reportKind = "rotavirus";
  else if (isRsvNextcladeReport(data)) reportKind = "rsv";
  else if (isMonkeypoxNextcladeReport(data)) reportKind = "monkeypox";
  else if (isInfluenzaTypingReport(data)) reportKind = "influenza";
  shell.dataset.reportKind = reportKind;
  if (isCommunityReport(data)) {
    renderCommunityReport(data);
    shell.dataset.reportReady = "true";
    return;
  }
  if (isPathoSourceReport(data)) {
    renderPathoSourceReport(data);
    shell.dataset.reportReady = "true";
    return;
  }
  const isMetaMethod = isMetaReport(data);
  const isVirusReport = isVirusFocusedReport(data);
  [
    document.getElementById("assembly-coverage-card"),
    document.getElementById("contig-depth-relationship-card"),
    document.getElementById("contig-length-depth-scatter-card"),
  ].forEach((node) => {
    if (!node) return;
    node.classList.toggle("hidden", isMetaMethod);
  });
  if (isMetaMethod) {
    [
      document.getElementById("section-checkm"),
      document.getElementById("section-gene-annotation"),
    ].forEach((node) => node?.remove());
    [
      document.querySelector('.report-nav-link[href="#section-checkm"]'),
      document.querySelector('.report-nav-link[href="#section-gene-annotation"]'),
    ].forEach((node) => node?.remove());
  }
  fillTaskMeta(data.task || {});
  if (applyMultiSampleLandingLayout(data)) {
    const titleNode = document.getElementById("report-sample-title");
    const copyNode = document.getElementById("report-sample-copy");
    if (titleNode) titleNode.textContent = data?.task?.name || data?.task?.id || "多样本结果概览";
    if (copyNode) copyNode.textContent = "当前处于批次概览模式；点击任一样本名可进入完整单样本报告。";
    renderMultiSampleOverview(data);
    shell.dataset.reportReady = "true";
    return;
  }
  applySarsCov2ReportChrome(data);
  applyRsvReportChrome(data);
  applyMonkeypoxReportChrome(data);
  applyInfluenzaReportChrome(data);
  bindReportSceneSwitcher(data);
  document.getElementById("overview-metrics").innerHTML = buildMetricCards(data.overview_metrics || []);
  renderRawQc(data.sections || {});
  renderTaxonomyRiskSummary(data.sections?.species_identification?.risk_summary || {});
  renderTaxonomyInterpretation(data.sections?.species_identification?.interpretation || {});
  renderSpeciesIdentification(data.sections?.species_identification || {});
  renderTaxonomyRarefaction(data.sections?.species_identification?.rarefaction || {});
  renderTaxonomyAbundance(data.sections?.species_identification?.abundance || {});
  renderBinningSection(data.task || {}, data.sections?.binning_results || {});
  buildTableCard("assembly-summary-table", "组装后信息统计", data.sections?.assembly?.summary?.columns || [], data.sections?.assembly?.summary?.rows || []);
  renderAssemblyCoverage(data.sections?.assembly?.coverage || {});
  renderContigDepthRelationship(data.sections?.assembly?.contig_depth_relationship || {});
  renderContigLengthDepthScatter(data.sections?.assembly?.contig_depth_relationship?.length_depth_scatter || {});
  if (!isMetaMethod && !isVirusReport) {
    buildTableCard("contig-annotation-table", "各个 Contig 注释结果", data.sections?.assembly?.contig_annotation?.columns || [], data.sections?.assembly?.contig_annotation?.rows || []);
    await renderCgviewSection(data.task || {}, data.sections?.assembly?.cgview || {}, data.sections || {});
  }
  if (!isMetaMethod && !isVirusReport) {
    buildTableCard("checkm-table", "CheckM 统计结果", data.sections?.assembly?.checkm?.columns || [], data.sections?.assembly?.checkm?.rows || []);
    buildTableCard("gene-annotation-summary-table", "基因注释统计", data.sections?.assembly?.gene_annotation_summary?.columns || [], data.sections?.assembly?.gene_annotation_summary?.rows || []);
  } else {
    const contigAnnotationSection = document.getElementById("section-contig-annotation");
    const checkmSection = document.getElementById("section-checkm");
    const geneAnnotationSection = document.getElementById("section-gene-annotation");
    const cgviewSection = document.getElementById("section-cgview");
    if (isVirusReport) {
      contigAnnotationSection?.remove();
      checkmSection?.remove();
      geneAnnotationSection?.remove();
      cgviewSection?.remove();
    }
    const checkmTable = document.getElementById("checkm-table");
    const geneAnnotationTable = document.getElementById("gene-annotation-summary-table");
    const contigAnnotationTable = document.getElementById("contig-annotation-table");
    const cgviewCard = document.getElementById("cgview-viewer-card");
    if (contigAnnotationTable) contigAnnotationTable.innerHTML = "";
    if (checkmTable) checkmTable.innerHTML = "";
    if (geneAnnotationTable) geneAnnotationTable.innerHTML = "";
    if (cgviewCard) cgviewCard.innerHTML = "";
  }
  renderGeneLengthDistribution(data.sections?.assembly?.gene_length_distribution || {});
  renderResistanceVirulenceOverview(data.sections?.resistance_virulence?.overview || {});
  renderNeisseriaAmrSection(data.sections?.tb_amr || data.sections?.mlst?.neisseria_amr || {});
  buildTableCard("rv-summary-table", "耐药毒力结果汇总", data.sections?.resistance_virulence?.summary?.columns || [], data.sections?.resistance_virulence?.summary?.rows || []);
  buildTableCard("virulence-table", "毒力元件", data.sections?.resistance_virulence?.virulence_elements?.columns || [], data.sections?.resistance_virulence?.virulence_elements?.rows || []);
  renderCategoryGeneRelationship("virulence-relationship-chart", data.sections?.resistance_virulence?.virulence_relationship || {});
  buildTableCard("resistance-table", "耐药元件", data.sections?.resistance_virulence?.resistance_elements?.columns || [], data.sections?.resistance_virulence?.resistance_elements?.rows || []);
  renderCategoryGeneRelationship("resistance-relationship-chart", data.sections?.resistance_virulence?.resistance_relationship || {});
  renderMlstSection(data.sections?.mlst || {});
  renderSerotypeSection(data.sections?.serotype || {});
  buildTableCard("priority-serotype-table", "关注毒力血清型", data.sections?.priority_serotype?.columns || [], data.sections?.priority_serotype?.rows || []);
  renderMgeOverview(data.sections?.mge_monitoring || {});
  buildTableCard("mge-resistance-table", "耐药相关移动元件监测", data.sections?.mge_monitoring?.resistance?.columns || [], data.sections?.mge_monitoring?.resistance?.rows || []);
  buildTableCard("mge-virulence-table", "毒力相关移动元件监测", data.sections?.mge_monitoring?.virulence?.columns || [], data.sections?.mge_monitoring?.virulence?.rows || []);
  initializeInteractiveCharts();
  bindChartExportButtons();
  initializeNextcladeGeneSummary();
  if (isSarsCov2NextcladeReport(data) || isHmpvNextcladeReport(data) || isDenvNextcladeReport(data) || isZikavNextcladeReport(data) || isChikvNextcladeReport(data) || isEbolaNextcladeReport(data) || isHpivTypingReport(data) || isHadvTypingReport(data) || isNorovirusTypingReport(data) || isEnterovirusTypingReport(data) || isHepatovirusTypingReport(data) || isHivTypingReport(data) || isBandavirusTypingReport(data) || isOrthohantavirusTypingReport(data) || isOrthoebolavirusTypingReport(data) || isAstroviridaeTypingReport(data) || isRhinovirusTypingReport(data) || isSeasonalHcovTypingReport(data) || isRotavirusTypingReport(data) || isRsvNextcladeReport(data) || isMonkeypoxNextcladeReport(data) || isInfluenzaTypingReport(data)) initializeReportNav();
  shell.dataset.reportReady = "true";
}

function initializeReportNav() {
  const groups = Array.from(document.querySelectorAll('[data-nav-group]'));
  const topLinks = Array.from(document.querySelectorAll('.report-nav > .report-nav-group > a.report-nav-link[href^="#"]'));
  const subLinks = Array.from(document.querySelectorAll('.report-subnav a[href^="#"]'));
  const toggles = Array.from(document.querySelectorAll('[data-nav-toggle]'));
  const navRoot = document.querySelector('.report-nav');
  const currentSectionLabel = document.getElementById('report-current-section');
  const previousNavState = window.__reportNavState;
  if (previousNavState) {
    window.removeEventListener('scroll', previousNavState.onScroll);
    window.removeEventListener('resize', previousNavState.onResize);
    window.removeEventListener('hashchange', previousNavState.onHashChange);
  }
  const collectNavTargetNodes = () => {
    const targetIds = new Set();
    [...topLinks, ...subLinks].forEach((link) => {
      const href = String(link.getAttribute('href') || '').trim();
      if (!href.startsWith('#') || href.length <= 1) return;
      targetIds.add(href.slice(1));
    });
    toggles.forEach((toggle) => {
      const targetId = String(toggle.dataset.navSection || '').trim();
      if (targetId) targetIds.add(targetId);
    });
    return Array.from(targetIds)
      .map((id) => document.getElementById(id))
      .filter((node) => node instanceof HTMLElement);
  };
  let syncScheduled = false;

  const openGroup = (group) => {
    const toggle = group.querySelector('[data-nav-toggle]');
    const subnav = group.querySelector('.report-subnav');
    if (!toggle || !subnav) return;
    group.classList.add('is-open');
    toggle.setAttribute('aria-expanded', 'true');
    subnav.hidden = false;
  };

  const closeGroup = (group) => {
    const toggle = group.querySelector('[data-nav-toggle]');
    const subnav = group.querySelector('.report-subnav');
    if (!toggle || !subnav) return;
    group.classList.remove('is-open');
    toggle.setAttribute('aria-expanded', 'false');
    subnav.hidden = true;
  };

  const closeOtherGroups = (activeGroup = null) => {
    groups.forEach((group) => {
      if (activeGroup instanceof HTMLElement && group === activeGroup) return;
      closeGroup(group);
    });
  };

  const clearActive = () => {
    [...topLinks, ...subLinks, ...toggles].forEach((node) => node.classList.remove('is-active'));
  };

  const revealNavNode = (node) => {
    if (!(node instanceof HTMLElement)) return;
    node.scrollIntoView({ block: 'nearest', inline: 'nearest' });
    if (navRoot instanceof HTMLElement) {
      navRoot.scrollIntoView({ block: 'nearest', inline: 'nearest' });
    }
  };

  const getNavNodeLabel = (node) => {
    if (!(node instanceof HTMLElement)) return "";
    return String(node.textContent || "")
      .replace(/\s+/g, " ")
      .trim();
  };

  const updateCurrentSectionLabel = (node, sectionId = "") => {
    if (!(currentSectionLabel instanceof HTMLElement)) return;
    const label = getNavNodeLabel(node) || (sectionId ? sectionId.replace(/^section-/, "") : "报告首页");
    currentSectionLabel.textContent = `当前位置：${label}`;
  };

  const syncActiveState = () => {
    syncScheduled = false;
    const threshold = 168;
    const sectionNodes = collectNavTargetNodes();
    const visibleSections = sectionNodes.filter((section) => {
      if (!(section instanceof HTMLElement)) return false;
      if (section.closest('.hidden')) return false;
      const style = window.getComputedStyle(section);
      return style.display !== 'none' && style.visibility !== 'hidden';
    });
    let current = visibleSections[0]?.id || '';
    let bestDistance = Number.POSITIVE_INFINITY;
    visibleSections.forEach((section) => {
      const rect = section.getBoundingClientRect();
      const topDelta = rect.top - threshold;
      if (topDelta <= 0) {
        const passedDistance = Math.abs(topDelta);
        if (passedDistance <= bestDistance) {
          bestDistance = passedDistance;
          current = section.id;
        }
        return;
      }
      if (!current || topDelta < bestDistance) {
        bestDistance = topDelta;
        current = section.id;
      }
    });
    clearActive();
    let matchedLink = document.querySelector(`.report-subnav a[href="#${current}"]`);
    if (matchedLink) {
      const group = matchedLink.closest('[data-nav-group]');
      if (group) openGroup(group);
      closeOtherGroups(group instanceof HTMLElement ? group : null);
      matchedLink.classList.add('is-active');
      const toggle = group?.querySelector('[data-nav-toggle]');
      toggle?.classList.add('is-active');
      updateCurrentSectionLabel(matchedLink, current);
      revealNavNode(toggle || matchedLink);
      return;
    }
    const topLink = document.querySelector(`.report-nav > .report-nav-group > a.report-nav-link[href="#${current}"]`);
    if (topLink) {
      closeOtherGroups();
      topLink.classList.add('is-active');
      updateCurrentSectionLabel(topLink, current);
      revealNavNode(topLink);
      return;
    }
    const toggle = toggles.find((node) => node.dataset.navSection === current);
    if (toggle) {
      const group = toggle.closest('[data-nav-group]');
      if (group) openGroup(group);
      closeOtherGroups(group instanceof HTMLElement ? group : null);
      toggle.classList.add('is-active');
      updateCurrentSectionLabel(toggle, current);
      revealNavNode(toggle);
    }
  };

  const scheduleSyncActiveState = () => {
    if (syncScheduled) return;
    syncScheduled = true;
    window.requestAnimationFrame(syncActiveState);
  };

  groups.forEach((group) => {
    const toggle = group.querySelector('[data-nav-toggle]');
    const subnav = group.querySelector('.report-subnav');
    if (!toggle || !subnav) return;
    closeGroup(group);
    if (toggle.dataset.navBound === '1') return;
    toggle.dataset.navBound = '1';
    toggle.addEventListener('click', () => {
      const isAlreadyOpen = group.classList.contains('is-open');
      if (isAlreadyOpen) closeGroup(group);
      else openGroup(group);
      clearActive();
      toggle.classList.add('is-active');
      if (!isAlreadyOpen) {
        revealNavNode(toggle);
        return;
      }
      const targetId = toggle.dataset.navSection;
      if (targetId) {
        const section = document.getElementById(targetId);
        if (section) {
          section.scrollIntoView({ behavior: 'smooth', block: 'start' });
          history.replaceState(null, '', `#${targetId}`);
          window.setTimeout(scheduleSyncActiveState, 80);
        }
      }
    });
  });

  [...topLinks, ...subLinks].forEach((link) => {
    if (link.dataset.navBound === '1') return;
    link.dataset.navBound = '1';
    link.addEventListener('click', () => {
      window.setTimeout(scheduleSyncActiveState, 40);
    });
  });

  window.addEventListener('scroll', scheduleSyncActiveState, { passive: true });
  window.addEventListener('resize', scheduleSyncActiveState, { passive: true });
  window.addEventListener('hashchange', scheduleSyncActiveState);
  window.__reportNavState = {
    onScroll: scheduleSyncActiveState,
    onResize: scheduleSyncActiveState,
    onHashChange: scheduleSyncActiveState,
  };
  scheduleSyncActiveState();
}

function initializeTopbarAutoHide() {
  const topbar = document.querySelector(".report-topbar");
  if (!(topbar instanceof HTMLElement)) return;

  let lastScrollY = window.scrollY || 0;
  let ticking = false;
  let pointerNearTop = false;
  const stableTopThreshold = 180;

  const getRevealZoneHeight = () => Math.max(56, Math.min(120, topbar.offsetHeight + 18));

  const syncTopbar = () => {
    const currentScrollY = window.scrollY || 0;
    if (currentScrollY <= stableTopThreshold) {
      topbar.classList.remove("is-collapsed");
      lastScrollY = currentScrollY;
      ticking = false;
      return;
    }
    const delta = currentScrollY - lastScrollY;
    const shouldCollapse = currentScrollY > stableTopThreshold && delta > 10;
    const shouldExpand = currentScrollY < 64 || (pointerNearTop && currentScrollY > 48);
    const shouldStayCollapsed = currentScrollY > stableTopThreshold && !pointerNearTop;

    if (shouldCollapse) {
      topbar.classList.add("is-collapsed");
    } else if (shouldExpand) {
      topbar.classList.remove("is-collapsed");
    } else if (shouldStayCollapsed) {
      topbar.classList.add("is-collapsed");
    }

    lastScrollY = currentScrollY;
    ticking = false;
  };

  window.addEventListener("scroll", () => {
    if (ticking) return;
    ticking = true;
    window.requestAnimationFrame(syncTopbar);
  }, { passive: true });

  window.addEventListener("mousemove", (event) => {
    const nextPointerNearTop = event.clientY <= getRevealZoneHeight();
    if (nextPointerNearTop === pointerNearTop) return;
    pointerNearTop = nextPointerNearTop;
    window.requestAnimationFrame(syncTopbar);
  }, { passive: true });

  window.addEventListener("mouseleave", () => {
    if (!pointerNearTop) return;
    pointerNearTop = false;
    window.requestAnimationFrame(syncTopbar);
  });
}


document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll(".report-tab-button").forEach((button) => {
    button.addEventListener("click", () => {
      const target = button.dataset.reportTab;
      document.querySelectorAll(".report-tab-button").forEach((item) => item.classList.toggle("active", item === button));
      document.querySelectorAll(".report-tab-panel").forEach((panel) => {
        panel.classList.toggle("active", panel.dataset.reportPanel === target);
      });
      initializeInteractiveCharts();
    });
  });
  document.querySelectorAll("[data-report-export-format]").forEach((button) => {
    button.addEventListener("click", async () => {
      try {
        await exportReportPage(button.dataset.reportExportFormat || "html");
      } catch (error) {
        console.error(error);
      }
    });
  });
  initializeReportSearch();
  initializeReportNav();
  initializeTopbarAutoHide();
  loadReport().catch((error) => {
    console.error(error);
  });
});
