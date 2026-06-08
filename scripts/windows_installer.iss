#define MyAppName "病原微生物分析工作台"
#define MyAppExeName "PathogenWorkbench"
#ifndef AppVersion
  #define AppVersion "0.1.0"
#endif
#ifndef ProjectRoot
  #define ProjectRoot "."
#endif
#ifndef AppDist
  #define AppDist AddBackslash(ProjectRoot) + "dist\" + MyAppExeName
#endif
#ifndef OutputDir
  #define OutputDir AddBackslash(ProjectRoot) + "dist_windows_installer"
#endif
#define SetupIconCandidate AddBackslash(ProjectRoot) + "bac_analysis_portal\static\app_icon.ico"

[Setup]
AppId={{8C9D7E0E-8D44-4D75-9BC7-4ED1A2D1B6F0}
AppName={#MyAppName}
AppVersion={#AppVersion}
AppPublisher=Metagenomic
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
OutputDir={#OutputDir}
OutputBaseFilename={#MyAppName}_Setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
PrivilegesRequired=admin
UninstallDisplayIcon={app}\{#MyAppExeName}.exe
#if FileExists(SetupIconCandidate)
SetupIconFile={#SetupIconCandidate}
#endif

[Languages]
Name: "chinesesimp"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "创建桌面快捷方式"; GroupDescription: "附加任务:"

[Files]
Source: "{#AppDist}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}.exe"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}.exe"; Description: "启动 {#MyAppName}"; Flags: nowait postinstall skipifsilent
