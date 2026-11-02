# Runs main.py in Abaqus/CAE and shows the progress log of the calculation
# live in the console (the print output of Abaqus/CAE is not shown while it
# runs). Called by run_abaqus.bat (-Mode noGUI) and run_abaqus_gui.bat
# (-Mode script).
#   -InputDir   input folder with Input.py (relative to the project folder or
#               absolute); passed to main.py in MONOPILE_INPUT_DIR
# main.py writes the log to the file given in the environment variable
# MONOPILE_PROGRESS_LOG (see methods/progress.py).
param(
    [ValidateSet("noGUI", "script")]
    [string]$Mode = "noGUI",
    [string]$InputDir = "Inputs"
)

Set-Location $PSScriptRoot

# Input folder: absolute path, must contain Input.py
if ([System.IO.Path]::IsPathRooted($InputDir)) {
    $inputPath = $InputDir
} else {
    $inputPath = Join-Path $PSScriptRoot $InputDir
}
if (-not (Test-Path (Join-Path $inputPath "Input.py"))) {
    Write-Host "Input folder '$InputDir' not found or without Input.py." `
        -ForegroundColor Red
    exit 2
}
$inputPath = (Resolve-Path $inputPath).Path
Write-Host "Input folder: $inputPath"
$env:MONOPILE_INPUT_DIR = $inputPath

$log = Join-Path $env:TEMP ("monopile_progress_{0}.log" -f $PID)
Remove-Item $log -ErrorAction SilentlyContinue
$env:MONOPILE_PROGRESS_LOG = $log

$proc = Start-Process -FilePath "cmd.exe" `
    -ArgumentList "/c abaqus cae $Mode=main.py" -NoNewWindow -PassThru
$null = $proc.Handle    # needed to read the exit code later

$shown = 0
function Show-NewLines {
    if (-not (Test-Path $log)) { return }
    $lines = @(Get-Content $log -ErrorAction SilentlyContinue)
    for ($i = $script:shown; $i -lt $lines.Count; $i++) {
        $line = $lines[$i]
        if ($line -match "FAILED|ABORTED|NOT converged|not converged") {
            Write-Host $line -ForegroundColor Red
        } elseif ($line -match "FINISHED|All calculations finished") {
            Write-Host $line -ForegroundColor Green
        } else {
            Write-Host $line
        }
        if ($Mode -eq "script" -and $line -match "All calculations finished") {
            Write-Host "Abaqus/CAE stays open; close it to end this window."
        }
    }
    $script:shown = $lines.Count
}

while (-not $proc.HasExited) {
    Show-NewLines
    Start-Sleep -Milliseconds 1000
}
Show-NewLines
Remove-Item $log -ErrorAction SilentlyContinue
exit $proc.ExitCode
