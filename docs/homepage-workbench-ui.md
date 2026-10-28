# 首页：紧凑型分析工作台

## 设计目标

- 首先选择分析模块，再设置输入与输出；高级参数按需调整。
- 保留五种分析模块、原有提交校验和模块权限，不以视觉状态代替校验结果。
- 用列表、分隔线和字级表达层次，减少嵌套卡片、渐变、投影和重复说明。

## 实现位置

- `bac_analysis_portal/templates/index.html`：分组导航、模块单选列表、四步进度和流程说明区域。
- `bac_analysis_portal/static/css/domains/submission-home.css`：以 `.workstation-ui` 为作用域覆盖旧版样式，不影响登录页。
- `bac_analysis_portal/static/workbench/domains/submission/home.js`：五种模块的静态流程说明；不决定分析参数或执行步骤。
- `bac_analysis_portal/static/app.js`：沿用原有模块与步骤状态，更新说明、当前步骤和首步摘要可见性。

## 交互验收

- 模块单选支持键盘；选中态、右侧说明和底部名称一致。
- 原生表单字段、步骤按钮、示例入口和权限控制继续沿用原来的 ID。
- 第一页隐藏尚未形成的任务摘要，后续步骤恢复显示。
- 手机端保留可横向滚动的全部导航入口，正文不产生横向溢出。
- 折叠导航提供可访问名称与悬停标题；当前步骤提供 `aria-current="step"`。
- 流程说明注明“实际步骤以所选运行方案为准”，不表示任务已完成或已校验。

## 验证命令

```sh
rtk proxy node --check bac_analysis_portal/static/app.js
rtk proxy node --check bac_analysis_portal/static/workbench/domains/submission/home.js
rtk proxy .venv_web/bin/python -m pytest -q tests/test_homepage_workbench_ui.py tests/test_release_task_submission.py tests/test_queue_closure_overview.py tests/test_closure_actions.py
rtk git diff --check
```

浏览器验收尺寸：1366×768、1440×900、1920×1080、2560×1440、3440×1440、390×844。
桌面检查五种模块和下一步按钮的首屏可见性；手机检查下一步可达及页面无横向溢出。
