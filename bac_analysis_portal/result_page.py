from __future__ import annotations

from urllib.parse import quote


def build_result_back_target(task_id: str, *, return_to: str = "", sample_key: str = "", database_section: str = "") -> dict[str, str]:
    normalized_task_id = str(task_id or "").strip()
    normalized_return_to = str(return_to or "").strip().lower()
    normalized_sample_key = str(sample_key or "").strip()
    normalized_database_section = str(database_section or "").strip() or "database-samples-panel"
    if normalized_return_to in {"batch", "overview", "multi"}:
        return {
            "href": f"/tasks/{quote(normalized_task_id, safe='')}/result-page?return_to=queue&task={quote(normalized_task_id, safe='')}",
            "label": "返回批次概览",
        }
    if normalized_return_to == "database":
        href = f"/workstation?tab=database&database_section={quote(normalized_database_section, safe='')}"
        if normalized_sample_key:
            href += f"&sample_key={quote(normalized_sample_key, safe='')}"
        return {"href": href, "label": "返回样本列表"}
    return {
        "href": f"/workstation?tab=queue&task={quote(normalized_task_id, safe='')}",
        "label": "返回任务列表",
    }

def _inject_result_preview_style(html: str, task_id: str = "") -> str:
    preview_style = """
<style id="portal-result-preview-style">
html, body {
  margin: 0 !important;
  padding: 0 !important;
  width: 100% !important;
  min-height: 100% !important;
  overflow-x: auto !important;
  background: #ffffff !important;
}
body {
  font-size: 16px !important;
}
.container-fluid.main-container,
div.main-container {
  width: min(1680px, calc(100vw - 96px)) !important;
  max-width: min(1680px, calc(100vw - 96px)) !important;
  margin: 0 auto !important;
  padding: 18px 24px 32px !important;
  box-sizing: border-box !important;
}
.container-fluid.main-container > .row {
  display: flex !important;
  align-items: flex-start !important;
  gap: 24px !important;
  margin: 0 !important;
}
.container-fluid.main-container > .row::before,
.container-fluid.main-container > .row::after {
  display: none !important;
}
.container-fluid.main-container > .row > [class*="col-"] {
  float: none !important;
  padding-left: 0 !important;
  padding-right: 0 !important;
}
.container-fluid.main-container > .row > .col-xs-12.col-sm-4.col-md-3 {
  flex: 0 0 280px !important;
  width: 280px !important;
  max-width: 280px !important;
}
.container-fluid.main-container > .row > .toc-content {
  flex: 1 1 auto !important;
  width: auto !important;
  max-width: none !important;
  min-width: 0 !important;
  padding: 0 !important;
}
div.tocify, #section-TOC {
  position: sticky !important;
  top: 82px !important;
  width: 280px !important;
  max-width: 280px !important;
  max-height: calc(100vh - 100px) !important;
  margin: 0 !important;
  overflow: auto !important;
}
.html-widget, .plotly, .datatables, table, img, svg, canvas {
  max-width: 100% !important;
}
pre {
  white-space: pre-wrap !important;
  word-break: break-word !important;
}
@media (max-width: 1024px) {
  .container-fluid.main-container,
  div.main-container {
    width: calc(100vw - 28px) !important;
    max-width: calc(100vw - 28px) !important;
    padding: 12px !important;
  }
  .container-fluid.main-container > .row {
    display: block !important;
  }
  .container-fluid.main-container > .row > .col-xs-12.col-sm-4.col-md-3,
  .container-fluid.main-container > .row > .toc-content {
    width: 100% !important;
    max-width: none !important;
  }
  div.tocify, #section-TOC {
    position: relative !important;
    top: auto !important;
    width: 100% !important;
    max-width: none !important;
    max-height: none !important;
    margin-bottom: 16px !important;
  }
}
</style>
"""
    if 'id="portal-result-preview-style"' in html:
        return html
    if '</head>' in html:
        return html.replace('</head>', f'{preview_style}</head>', 1)
    return preview_style + html

def _inject_result_page_chrome(html: str, task: dict, *, back_href: str = "", back_label: str = "") -> str:
    title = str(task.get("name") or task.get("id") or "任务结果")
    status = str(task.get("status") or "未知")
    target = build_result_back_target(str(task.get("id") or "")) if not back_href else {"href": back_href, "label": back_label or "返回任务列表"}
    chrome_style = f"""
<style id="portal-result-page-chrome">
body {{
  padding-top: 64px !important;
  scroll-padding-top: 76px !important;
}}
#portal-result-page-bar {{
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  z-index: 9999;
  height: 64px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
  padding: 0 18px;
  border-bottom: 1px solid rgba(35, 54, 88, 0.12);
  background: rgba(246, 248, 252, 0.98);
  backdrop-filter: blur(10px);
  box-sizing: border-box;
}}
#portal-result-page-left {{
  display: flex;
  align-items: center;
  gap: 12px;
  min-width: 0;
}}
#portal-result-float-back {{
  position: fixed;
  top: 12px;
  left: 18px;
  z-index: 10002;
  border: 1px solid rgba(35, 54, 88, 0.14);
  background: #ffffff;
  color: #17233a;
  border-radius: 14px;
  padding: 10px 14px;
  font: inherit;
  cursor: pointer;
  box-shadow: 0 10px 26px rgba(23, 35, 58, 0.12);
}}
#portal-result-title {{
  font-size: 20px;
  font-weight: 700;
  color: #17233a;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}}
#portal-result-status {{
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-width: 96px;
  height: 40px;
  padding: 0 14px;
  border-radius: 999px;
  background: #355c94;
  color: #fff;
  font-weight: 700;
  font-size: 14px;
}}
@media (max-width: 720px) {{
  body {{
    padding-top: 84px !important;
  }}
  #portal-result-page-bar {{
    height: 84px;
    padding: 10px 12px;
    align-items: flex-start;
  }}
  #portal-result-page-left {{
    flex-direction: column;
    align-items: flex-start;
    gap: 8px;
    padding-top: 38px;
  }}
  #portal-result-float-back {{
    top: 10px;
    left: 12px;
    padding: 9px 12px;
  }}
}}
</style>
<script>
window.addEventListener('DOMContentLoaded', function () {{
  var backButton = document.getElementById('portal-result-float-back');
    if (backButton) {{
      backButton.addEventListener('click', function () {{
      window.location.href = {target["href"]!r}
    }});
  }}
}});
</script>
"""
    chrome_body = f"""
<button id="portal-result-float-back" type="button">{target["label"]}</button>
<div id="portal-result-page-bar">
  <div id="portal-result-page-left">
    <div id="portal-result-title">{title}</div>
  </div>
  <div id="portal-result-status">{status}</div>
</div>
"""
    if '</head>' in html and 'portal-result-page-chrome' not in html:
        html = html.replace('</head>', chrome_style + '\n</head>', 1)
    if '<body>' in html and 'portal-result-page-bar' not in html:
        html = html.replace('<body>', '<body>\n' + chrome_body, 1)
    return html

def _build_result_placeholder_page(task: dict, *, back_href: str = "", back_label: str = "") -> str:
    title = str(task.get("name") or task.get("id") or "任务结果")
    status = str(task.get("status") or "未知")
    task_id = str(task.get("id") or "")
    target = build_result_back_target(task_id) if not back_href else {"href": back_href, "label": back_label or "返回任务列表"}
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    body {{ margin: 0; font-family: sans-serif; background: #f3f6fb; color: #17233a; }}
    .bar {{ display:flex; justify-content:space-between; align-items:center; padding:12px 18px; border-bottom:1px solid rgba(35,54,88,.12); background:rgba(246,248,252,.98); }}
    .left {{ display:flex; align-items:center; gap:12px; }}
    button {{ border:1px solid rgba(35,54,88,.14); background:#fff; border-radius:12px; padding:10px 14px; cursor:pointer; }}
    .status {{ display:inline-flex; align-items:center; justify-content:center; min-width:96px; height:40px; padding:0 14px; border-radius:999px; background:#6a7488; color:#fff; font-weight:700; }}
    .body {{ padding:40px 24px; }}
  </style>
</head>
<body>
  <div class="bar">
    <div class="left"><button onclick="window.location.href={target["href"]!r}">{target["label"]}</button><strong>{title}</strong></div>
    <div class="status">{status}</div>
  </div>
  <div class="body">当前没有可展示的结果页面。</div>
</body>
</html>"""
