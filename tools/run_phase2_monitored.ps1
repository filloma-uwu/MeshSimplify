param(
    [string]$Executable = "$PSScriptRoot\..\build_cpp_phase1\Release\pqss-hausdorff-simplify.exe",
    [Parameter(Mandatory = $true)][string]$Phase1Halfedge,
    [Parameter(Mandatory = $true)][string]$Source,
    [Parameter(Mandatory = $true)][string]$OutputDirectory,
    [Parameter(Mandatory = $true)][string]$ModelId,
    [Parameter(Mandatory = $true)][double]$MaximumDirectedHausdorff,
    [int]$AbortWorkingSetMiB = 1850,
    [switch]$NoCoplanarUnion,
    [switch]$NoPlanarConvexification,
    [switch]$NoConvexHull,
    [switch]$NoAdaptiveConvexCover,
    [switch]$NoBoundedBoxCover,
    [switch]$NoKdop,
    [switch]$NoBox
)

$ErrorActionPreference = "Stop"
$executablePath = (Resolve-Path $Executable).Path
$halfedgePath = (Resolve-Path $Phase1Halfedge).Path
$sourcePath = (Resolve-Path $Source).Path
$outputPath = [IO.Path]::GetFullPath($OutputDirectory)
if (Test-Path -LiteralPath $outputPath) {
    throw "Output directory already exists: $outputPath"
}
$logRoot = Split-Path $outputPath -Parent
if (-not (Test-Path -LiteralPath $logRoot)) {
    New-Item -ItemType Directory -Path $logRoot -Force | Out-Null
}
$stdoutPath = "$outputPath.stdout.txt"
$stderrPath = "$outputPath.stderr.txt"
$peakPath = "$outputPath.peak_memory_mib.txt"

$arguments = @(
    "--phase1-halfedge", $halfedgePath,
    "--source", $sourcePath,
    "--output-dir", $outputPath,
    "--model-id", $ModelId,
    "--maximum-directed-hausdorff",
    $MaximumDirectedHausdorff.ToString("R", [Globalization.CultureInfo]::InvariantCulture)
)
if ($NoCoplanarUnion) { $arguments += "--no-coplanar-union" }
if ($NoPlanarConvexification) { $arguments += "--no-planar-convexification" }
if ($NoConvexHull) { $arguments += "--no-convex-hull" }
if ($NoAdaptiveConvexCover) { $arguments += "--no-adaptive-convex-cover" }
if ($NoBoundedBoxCover) { $arguments += "--no-bounded-box-cover" }
if ($NoKdop) { $arguments += "--no-kdop" }
if ($NoBox) { $arguments += "--no-box" }

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
    Set-Content $peakPath -Encoding ASCII

if ($aborted) {
    throw "phase2 aborted at working-set limit: $AbortWorkingSetMiB MiB; peak=$peakMiB MiB"
}
if (-not (Test-Path -LiteralPath (Join-Path $outputPath "model.json"))) {
    throw "phase2 failed with exit code $($process.ExitCode); peak=$peakMiB MiB; stderr=$stderrPath"
}
Write-Output "phase2_completed peak_memory_mib=$peakMiB output=$outputPath"
