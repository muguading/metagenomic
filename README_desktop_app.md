# 桌面版启动

当前项目已经支持通过 `pywebview + Flask` 方式作为桌面应用运行。

## 安装依赖

```bash
./.venv_web/bin/python -m pip install -r requirements-web.txt
```

## 启动桌面 App

```bash
./.venv_web/bin/python run_bac_analysis_desktop.py
```

启动后会：

- 先显示统一的连接与登录入口
- 可选择连接服务器，或本机启动
- 若连接服务器，会直接跳到远端 `/login` 并自动提交账号密码
- 若本机启动，会先起本地服务再自动登录

## 一键打包 macOS App

```bash
bash scripts/build_mac_desktop_app.sh
```

如果你希望打包完成后立即打开应用：

```bash
bash scripts/build_mac_desktop_app.sh --open
```

脚本会自动完成：
- 用现有 `病原微生物分析工作台.spec` 重新打包
- 生成 `dist/病原微生物分析工作台.app`
- 覆盖安装到 `/Applications/病原微生物分析工作台.app`

首次打包前请确认：

- `requirements-web.txt` 中依赖已安装完成
- 优先准备好 `./.venv_web` 虚拟环境，脚本会优先使用它
- macOS 当前用户对 `/Applications` 有写权限

## Windows 打包准备

Windows 版建议在 Windows 10/11 机器上构建，不要在 macOS 上直接交叉打包。

### 1. 安装依赖

在 Windows PowerShell 中执行：

```powershell
python -m pip install -r requirements-web.txt
```

如果你希望固定到项目虚拟环境：

```powershell
python -m venv .venv_web
.\.venv_web\Scripts\python.exe -m pip install -r requirements-web.txt
```

### 2. 生成 Windows 桌面程序

```powershell
powershell -ExecutionPolicy Bypass -File .\build_windows_desktop_app.ps1
```

脚本会：

- 使用 `desktop_app_windows.spec`
- 生成 `dist\PathogenWorkbench\PathogenWorkbench.exe`
- 如果系统已安装 Inno Setup，并且命令 `iscc` 可用，会继续生成安装包

如果你只想先生成可执行目录，不做安装包：

```powershell
powershell -ExecutionPolicy Bypass -File .\build_windows_desktop_app.ps1 -SkipInstaller
```

### 3. 生成 Windows 安装包

安装 [Inno Setup](https://jrsoftware.org/isinfo.php) 后，脚本会自动调用：

- 安装脚本模板：`scripts/windows_installer.iss`
- 默认输出目录：`dist_windows_installer`

### 4. Windows 图标说明

当前仓库内只有 `app_icon.png` 和 `app_icon.icns`。

- 如果没有 `bac_analysis_portal/static/app_icon.ico`，Windows 可执行文件和安装包仍可构建
- 但会使用默认图标，不够像正式交付件

建议在 Windows 机器上补一个：

`bac_analysis_portal/static/app_icon.ico`

这样 `PyInstaller` 和 `Inno Setup` 会自动使用它
