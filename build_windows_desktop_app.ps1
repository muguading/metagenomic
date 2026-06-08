param(
    [switch]$SkipInstaller
)

$ScriptPath = Join-Path $PSScriptRoot "scripts\build_windows_desktop_app.ps1"

if (!(Test-Path $ScriptPath)) {
    throw "Script was not found: $ScriptPath"
}

& powershell -ExecutionPolicy Bypass -File $ScriptPath @PSBoundParameters
