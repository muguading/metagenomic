param(
    [switch]$SkipInstaller
)

$ErrorActionPreference = "Stop"

$ProjectRoot = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$AppName = "PathogenWorkbench"
$SpecFile = Join-Path $ProjectRoot "desktop_app_windows.spec"
$RequirementsFile = Join-Path $ProjectRoot "requirements-web.txt"
$VenvPython = Join-Path $ProjectRoot ".venv_web\Scripts\python.exe"
$PyInstallerDist = Join-Path $ProjectRoot "dist\$AppName"
$InstallerScript = Join-Path $ProjectRoot "scripts\windows_installer.iss"
$InstallerOutput = Join-Path $ProjectRoot "dist_windows_installer"

function Resolve-PythonCommand {
  if (Test-Path $VenvPython) {
    return @{
      Executable = $VenvPython
      PrefixArgs = @()
      Display = $VenvPython
    }
  }

  $python = Get-Command python -ErrorAction SilentlyContinue
  if ($python) {
    return @{
      Executable = $python.Source
      PrefixArgs = @()
      Display = $python.Source
    }
  }

  $pyLauncher = Get-Command py -ErrorAction SilentlyContinue
  if ($pyLauncher) {
    return @{
      Executable = $pyLauncher.Source
      PrefixArgs = @("-3")
      Display = "$($pyLauncher.Source) -3"
    }
  }

  throw "Python was not found. Install Python 3 or create .venv_web first."
}

$PythonCommand = Resolve-PythonCommand

function Invoke-Python {
  param(
    [string[]]$PythonArgs
  )
  $allArgs = @()
  if ($PythonCommand.PrefixArgs) {
    $allArgs += $PythonCommand.PrefixArgs
  }
  if ($PythonArgs) {
    $allArgs += $PythonArgs
  }
  & $PythonCommand.Executable @allArgs
}

Write-Host "==> Project root: $ProjectRoot"
Write-Host "==> Python: $($PythonCommand.Display)"
Write-Host "==> Checking dependencies..."

try {
  Invoke-Python -PythonArgs @("-c", "import PyInstaller, flask, webview")
} catch {
  throw "Missing packaging dependencies. Run: python -m pip install -r $RequirementsFile"
}

if (!(Test-Path $SpecFile)) {
  throw "Windows spec file was not found: $SpecFile"
}

$IcoPath = Join-Path $ProjectRoot "bac_analysis_portal\static\app_icon.ico"
if (!(Test-Path $IcoPath)) {
  Write-Warning "app_icon.ico was not found. The Windows executable will use the default icon."
}

Write-Host "==> Building Windows desktop app..."
Invoke-Python -PythonArgs @("-m", "PyInstaller", $SpecFile, "--noconfirm")

if (!(Test-Path $PyInstallerDist)) {
  throw "Build finished but output folder was not found: $PyInstallerDist"
}

Write-Host "==> Windows app output: $PyInstallerDist"

if ($SkipInstaller) {
  Write-Host "==> Installer generation skipped"
  exit 0
}

$Iscc = Get-Command iscc -ErrorAction SilentlyContinue
if (-not $Iscc) {
  Write-Warning "Inno Setup compiler (iscc) was not found. Installer generation was skipped."
  Write-Warning "After installing Inno Setup, run: iscc `"$InstallerScript`""
  exit 0
}

if (!(Test-Path $InstallerOutput)) {
  New-Item -ItemType Directory -Path $InstallerOutput | Out-Null
}

Write-Host "==> Building Windows installer..."
& $Iscc.Source "/DAppVersion=0.1.0" "/DProjectRoot=$ProjectRoot" "/DAppDist=$PyInstallerDist" "/DOutputDir=$InstallerOutput" $InstallerScript

Write-Host "==> Installer output: $InstallerOutput"
