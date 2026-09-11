from pathlib import Path


ROOT = Path(__file__).parents[1]
APP_JS = (ROOT / "bac_analysis_portal" / "static" / "app.js").read_text(encoding="utf-8")
STYLES = (ROOT / "bac_analysis_portal" / "static" / "styles.css").read_text(encoding="utf-8")


def test_prediction_view_shows_summary_distributions_and_detail_table() -> None:
    assert "预测标签分布" in APP_JS
    assert "关注级别分布" in APP_JS
    assert "本批次研判结论" in APP_JS
    assert "置信度只表示模型对结论的确定程度" in APP_JS
    assert "后续研判" in APP_JS
    assert "变化证据" in APP_JS
    assert "显著偏移点" in APP_JS
    assert "偏移点明细" in APP_JS
    assert "历史常规区间" in APP_JS
    assert "只能判断单个样本类别" in APP_JS
    assert "自动发现偏移信号" in APP_JS
    assert "扫描对象" in APP_JS
    assert "偏移信号雷达" in APP_JS
    assert "signal_id" in APP_JS
    assert "预测结果明细" in APP_JS
    assert "主要依据" in APP_JS
    assert "modeling-prediction-page-size" in APP_JS
    assert "modeling-prediction-table" in APP_JS


def test_prediction_visualization_has_structured_styles() -> None:
    assert ".modeling-prediction-bars" in STYLES
    assert ".modeling-prediction-risk" in STYLES
    assert ".modeling-prediction-toolbar" in STYLES
    assert ".modeling-prediction-conclusion" in STYLES
    assert ".modeling-prediction-drivers" in STYLES
    assert ".modeling-trend-comparison" in STYLES
    assert ".modeling-trend-chart" in STYLES
    assert ".modeling-trend-indicators" in STYLES
    assert ".modeling-trend-deviation-table" in STYLES
    assert ".modeling-prediction-model-capability" in STYLES
    assert ".modeling-scan-types" in STYLES
    assert ".modeling-signal-radar-table" in STYLES
