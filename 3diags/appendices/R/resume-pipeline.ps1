# Resume interrupted full pipeline (retros → collect → hindcasts → ensemble → profiles → render → verify).
# Skips retro peels that already have Report.sso (r4ss::retro default).
param(
  [string]$AppendixDir = "appendices",
  [string]$Study = $(if ($env:DIAG_STUDY) { $env:DIAG_STUDY } else { "sma-2026" }),
  [string]$Ss3Exe = $(if ($env:SS3_EXE) { $env:SS3_EXE } else { "ss3.exe" })
)

$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..\..")
$env:SS3_EXE = $Ss3Exe
$env:DIAG_STUDY = $Study

$logDir = Join-Path $AppendixDir "data/pipeline"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$log = Join-Path $logDir "resume.log"

function Log-Step($msg) {
  $line = "=== $msg $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ==="
  Write-Host $line
  Add-Content -Path $log -Value $line
}

function Run-R($label, [string[]]$RArgs) {
  Log-Step $label
  $cmd = "Rscript $($RArgs -join ' ')"
  Add-Content -Path $log -Value $cmd
  $prevEap = $ErrorActionPreference
  $ErrorActionPreference = "Continue"
  & Rscript @RArgs *>&1 | Tee-Object -FilePath $log -Append
  $exit = $LASTEXITCODE
  $ErrorActionPreference = $prevEap
  if ($exit -ne 0) { throw "$label failed (exit $exit)" }
}

Log-Step "Resume pipeline start (study=$Study, SS3_EXE=$Ss3Exe)"

$cfgLines = & Rscript (Join-Path $AppendixDir "R/read-study-paths.R") --appendix-dir $AppendixDir --study $Study
$parts = $cfgLines -split '\|'
$ssBase = $parts[0]
$retroScenarios = $parts[1] -split ','
$nRetroYears = 11

$studyArgs = @("--study", $Study)
$retroDone = $true
foreach ($sc in $retroScenarios) {
  $done = (Get-ChildItem "$ssBase/$sc/retrospectives" -Directory -EA SilentlyContinue |
    Where-Object { Test-Path (Join-Path $_.FullName "Report.sso") }).Count
  if ($done -lt $nRetroYears) { $retroDone = $false }
}
$oneStep = Join-Path $AppendixDir "data/retrospectives/one-step.csv"

if (-not $retroDone) {
  Run-R "Retrospectives --run" @(
    (Join-Path $AppendixDir "R/run-retrospectives.R"),
    "--appendix-dir", $AppendixDir,
    @studyArgs,
    "--run"
  )
} else {
  Log-Step "Retrospectives --run skipped (peels complete)"
}

if (-not (Test-Path $oneStep)) {
  Run-R "Retrospectives --collect" @(
    (Join-Path $AppendixDir "R/run-retrospectives.R"),
    "--appendix-dir", $AppendixDir,
    @studyArgs,
    "--collect"
  )
} else {
  Log-Step "Retrospectives --collect skipped (one-step.csv exists)"
}

Run-R "Pipeline remainder (--skip-retros)" @(
  (Join-Path $AppendixDir "R/run-full-pipeline.R"),
  "--appendix-dir", $AppendixDir,
  @studyArgs,
  "--skip-retros"
)

Log-Step "Resume pipeline complete"
