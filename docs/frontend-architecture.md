# 前端模块架构

工作台使用原生 ES Modules，不引入前端构建步骤。`static/app.js` 是组合入口，负责创建
runtime、注入共享能力并注册域模块；业务事件与生命周期进入 `static/workbench/domains/`。

## 模块归属

- `workbench/core/api.js`：JSON 请求、错误解析和 CSRF 请求入口。
- `workbench/core/ui.js`：无业务归属的 UI 工具。
- `workbench/core/navigation.js`：Tab 生命周期注册与激活。
- `workbench/core/path-browser.js`：路径浏览器目标注册与选择分发。
- `workbench/domains/project/`：项目管理。
- `workbench/domains/admin/`：后台管理。
- `workbench/domains/audit/`：审计追踪。

域模块统一暴露 `init(runtime)`、`activate()` 和 `refresh()`。域模块不得互相导入，只能通过
注入的 `runtime.selectors` 和 `runtime.actions` 使用共享状态与动作。新增业务状态不得挂到
新的全局变量；新增域事件绑定不得写回 `app.js`。

项目域的 `management.js` 拥有项目聚合、甘特图、负责人、里程碑和归档行为；
`storage.js` 负责兼容既有 localStorage key。项目实现不得重新写回入口。

需要使用文件浏览器的域通过 `runtime.pathBrowser.registerTarget(selector, handler)` 注册目标；
通用路径浏览器不得包含后台、项目或其他业务字段的赋值特判。

## 模板与样式

已迁移域的模板放在 `templates/workbench/`，由 `index.html` 使用 Jinja partial 引入。
单域样式放在 `static/css/domains/`；共享控件和跨域选择器继续保留在 `styles.css`。

当一个域文件难以快速理解或超过约 1000 行时，按单一工作流继续拆子模块。禁止为了拆分
制造只转发一次调用的浅模块。

## 质量闸门

`frontend-size-budgets.json` 记录旧巨型文件的迁移后行数上限。更新预算只能伴随文件继续
缩小，不得用提高预算绕过检查。

统一验证命令：

```bash
python -m pytest tests -q
find bac_analysis_portal/static/workbench -name "*.js" -print0 |
  while IFS= read -r -d "" file; do node --input-type=module --check < "$file"; done
node --input-type=module --check < bac_analysis_portal/static/app.js
```
