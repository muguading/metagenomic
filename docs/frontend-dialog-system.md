# 前端弹窗系统

本项目弹窗统一从 `bac_analysis_portal/static/app.js` 的工作台弹窗工具派生。新增确认、输入、表单类弹窗时，不要直接拼一套新的 `browser-modal` 外壳。

## 统一入口

| 工具 | 用途 |
| --- | --- |
| `openWorkbenchConfirmModal(options)` | 确认类弹窗，返回 `Promise<boolean>` |
| `openWorkbenchFormModal(options)` | 表单类弹窗，返回 `Promise<payload \| null>` |
| `showModalElement(modal)` | 打开已有静态弹窗，统一处理 `hidden`、`aria-hidden` 与滚动锁 |
| `hideModalElement(modal)` | 关闭已有静态弹窗，统一释放滚动锁 |

## 使用规则

1. 不使用原生 `alert`、`confirm`、`prompt`。
2. 不直接操作弹窗的 `aria-hidden` 和 `modal-scroll-locked`。
3. 静态模板弹窗必须通过 `showModalElement` / `hideModalElement` 打开关闭。
4. 动态确认和输入弹窗必须通过 `openWorkbenchConfirmModal` / `openWorkbenchFormModal`。
5. 删除、停止、清空等高风险动作使用 `tone: "danger"`。
6. 会丢失修改、覆盖状态或需要人工确认的动作使用 `tone: "warning"`。
7. 表单校验放在 `validateSubmit`，返回值整理放在 `transformSubmit`。

## 审计命令

```bash
rtk grep -n "\\balert\\(|\\bconfirm\\(|\\bprompt\\(" bac_analysis_portal/static/app.js bac_analysis_portal/templates/index.html
rtk grep -n "setAttribute\\(\\\"aria-hidden\" bac_analysis_portal/static/app.js
rtk grep -n "Modal\\?\\.classList\\.(add|remove)\\(\\\"hidden\\\"\\)|Modal\\.classList\\.(add|remove)\\(\\\"hidden\\\"\\)" bac_analysis_portal/static/app.js
rtk grep -n "modal-scroll-locked" bac_analysis_portal/static/app.js
```

通过标准：除统一工具函数内部外，不应出现直接开关弹窗可见性、滚动锁或原生弹窗调用。
