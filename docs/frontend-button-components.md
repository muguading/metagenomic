# 前端按钮组件库

本项目按钮统一从 `bac_analysis_portal/static/styles.css` 的 button token 与 `.ui-button` 组件层派生。新增页面优先使用这里的类名；历史类名已做兼容映射，逐步迁移即可。

## 基础用法

```html
<button class="ui-button ui-button--primary" type="button">保存</button>
<button class="ui-button ui-button--secondary" type="button">取消</button>
<button class="ui-button ui-button--warning" type="button">需要关注</button>
<button class="ui-button ui-button--danger" type="button">删除</button>
```

## 语义变体

| 类名 | 用途 |
| --- | --- |
| `.ui-button--primary` | 页面主动作，每个操作区尽量只出现一个 |
| `.ui-button--secondary` | 常规动作、工具条动作、弹窗次动作 |
| `.ui-button--quiet` | 低强调动作，适合列表内弱操作 |
| `.ui-button--warning` | 需要关注、可能改变闭环状态但非删除的动作 |
| `.ui-button--danger` | 删除、终止、清空等不可逆或高风险动作 |
| `.ui-button--link` | 行内文字动作 |

## 尺寸与布局

| 类名 | 用途 |
| --- | --- |
| `.ui-button--sm` | 工具条、表格行内操作 |
| `.ui-button--lg` | 面板头部或弹窗主动作 |
| `.ui-button--pill` | 胶囊型短按钮 |
| `.ui-button--block` | 移动端或表单中铺满容器 |
| `.ui-button--icon` | 纯图标按钮 |

## 兼容类

以下历史类名已经接入按钮 token，可继续工作，但新代码应优先使用 `.ui-button`：

| 历史类名 | 等价语义 |
| --- | --- |
| `.primary-button` | `.ui-button ui-button--primary` |
| `.ghost-button` | `.ui-button ui-button--secondary` |
| `.ghost-button.danger` | `.ui-button ui-button--danger` |
| `.danger-button` | `.ui-button ui-button--danger` |

## 迁移规则

1. 新增按钮必须写 `type="button"`，提交按钮除外。
2. 不要在业务模块里重新定义按钮颜色、阴影和禁用态。
3. 模块内只允许覆盖尺寸、布局宽度或特殊排列。
4. 危险按钮必须使用 danger 语义，不要只靠红色文本表达风险。
5. 手机端不要隐藏关键动作；使用折叠、底部操作区或 `.ui-button--block` 适配。

## 审计清单

迁移旧按钮时按下面顺序判断：

1. 是否只是普通动作：改为 `.ui-button ui-button--secondary`。
2. 是否是当前区域唯一主动作：改为 `.ui-button ui-button--primary`。
3. 是否是删除、停止、清空：改为 `.ui-button ui-button--danger`。
4. 是否是行内轻量操作：改为 `.ui-button ui-button--quiet ui-button--sm`。
5. 是否是文字链接：改为 `.ui-button ui-button--link` 或 `.inline-link-button`。

可用下面的搜索辅助发现残留自绘按钮：

```bash
rtk grep -n "button.*background\\|button.*box-shadow\\|button.*border-color" bac_analysis_portal/static/styles.css
rtk grep -n "class=\\\"[^\"]*button" bac_analysis_portal/templates/index.html bac_analysis_portal/static/app.js
```

发现残留后，优先删除颜色、阴影、禁用态覆盖；只保留模块需要的宽度、高度、排列和特殊图标结构。
