param(
    [string]$Executable = "build\Release\pqss-phase2-repair-halfedge.exe",
    [Parameter(Mandatory = $true)][string]$Phase1Halfedge,
    [Parameter(Mandatory = $true)][string]$OutputRoot,
    [Parameter(Mandatory = $true)][string]$ModelId,
    [Parameter(Mandatory = $true)][double]$MaximumDirectedHausdorff,
    [int]$AbortWorkingSetMiB = 1850
)

$ErrorActionPreference = "Stop"
$executablePath = (Resolve-Path -LiteralPath $Executable).Path
$halfedgePath = (Resolve-Path -LiteralPath $Phase1Halfedge).Path
$outputRootPath = [IO.Path]::GetFullPath($OutputRoot)
$outputPath = Join-Path (Join-Path $outputRootPath "models") $ModelId

if (Test-Path -LiteralPath $outputPath) {
    throw "Output directory already exists: $outputPath"
}
if (-not (Test-Path -LiteralPath $outputRootPath)) {
    New-Item -ItemType Directory -Path $outputRootPath -Force | Out-Null
}
$logRoot = Split-Path -Path $outputPath -Parent
if (-not (Test-Path -LiteralPath $logRoot)) {
    New-Item -ItemType Directory -Path $logRoot -Force | Out-Null
}

$stdoutPath = "$outputPath.stdout.txt"
$stderrPath = "$outputPath.stderr.txt"
$peakPath = "$outputPath.peak_memory_mib.txt"

$arguments = @(
    "--phase1-halfedge", $halfedgePath,
    "--output-root", $outputRootPath,
    "--model-id", $ModelId,
    "--maximum-directed-hausdorff",
    $MaximumDirectedHausdorff.ToString("R", [Globalization.CultureInfo]::InvariantCulture)
)

$process = Start-Process -FilePath $executablePath -ArgumentList $arguments `
    -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath `
    -WindowStyle Hidden -PassThru

$peakBytes = 0L
$aborted = $false
while (-not $process.HasExited) {
    Start-Sleep -Milliseconds 250
    $process.Refresh()
    if ($process.HasExited) { break }
    $peakBytes = [Math]::Max($peakBytes, $process.WorkingSet64)
    if ($process.WorkingSet64 -ge $AbortWorkingSetMiB * 1MB) {
        Stop-Process -Id $process.Id -Force
        $aborted = $true
        break
    }
}
$process.WaitForExit()
$process.Refresh()
$peakBytes = [Math]::Max($peakBytes, $process.PeakWorkingSet64)
$peakMiB = [double]($peakBytes / 1MB)
$peakMiB.ToString("R", [Globalization.CultureInfo]::InvariantCulture) |
    Set-Content -LiteralPath $peakPath -Encoding ASCII

if ($aborted) {
    throw "phase2 repair aborted at working-set limit: $AbortWorkingSetMiB MiB; peak=$peakMiB MiB; stderr=$stderrPath"
}
if (-not (Test-Path -LiteralPath (Join-Path $outputPath "model.json"))) {
    throw "phase2 repair failed with exit code $($process.ExitCode); peak=$peakMiB MiB; stderr=$stderrPath"
}

Write-Output "phase2_repair_completed peak_memory_mib=$peakMiB output=$outputPath"
