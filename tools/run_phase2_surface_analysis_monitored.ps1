param(
    [string]$Executable = "$PSScriptRoot\..\build_cpp_phase1\Release\pqss-halfedge-surface-analyze.exe",
    [Parameter(Mandatory = $true)][string]$Phase1Halfedge,
    [Parameter(Mandatory = $true)][string]$OutputDirectory,
    [Parameter(Mandatory = $true)][double]$MaximumDirectedHausdorff,
    [int]$AbortWorkingSetMiB = 500,
    [int]$MaximumRuntimeSeconds = 300,
    [int]$NoProgressTimeoutSeconds = 120
)

$ErrorActionPreference = "Stop"
$executablePath = (Resolve-Path $Executable).Path
$halfedgePath = (Resolve-Path $Phase1Halfedge).Path
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
    "--output-dir", $outputPath,
    "--maximum-directed-hausdorff",
    $MaximumDirectedHausdorff.ToString("R", [Globalization.CultureInfo]::InvariantCulture)
)
$process = Start-Process -FilePath $executablePath -ArgumentList $arguments `
  -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath `
  -WindowStyle Hidden -PassThru

$peakBytes = 0L
$aborted = $false
$abortReason = ""
$started = [DateTime]::UtcNow
$lastProgress = $started
$lastLogWrite = [DateTime]::MinValue
$nextHeartbeat = $started.AddSeconds(10)
while (-not $process.HasExited) {
    Start-Sleep -Milliseconds 100
    $process.Refresh()
    if ($process.HasExited) { break }
    $peakBytes = [Math]::Max($peakBytes, $process.WorkingSet64)
    $now = [DateTime]::UtcNow
    if (Test-Path -LiteralPath $stderrPath) {
        $writeTime = (Get-Item -LiteralPath $stderrPath).LastWriteTimeUtc
        if ($writeTime -gt $lastLogWrite) {
            $lastLogWrite = $writeTime
            $lastProgress = $now
        }
    }
    if ($now -ge $nextHeartbeat) {
        $elapsed = [int]($now - $started).TotalSeconds
        $workingSetMiB = [Math]::Round($process.WorkingSet64 / 1MB, 1)
        Write-Output "phase2_heartbeat elapsed_seconds=$elapsed working_set_mib=$workingSetMiB"
        $nextHeartbeat = $now.AddSeconds(10)
    }
    if ($process.WorkingSet64 -ge $AbortWorkingSetMiB * 1MB) {
        $aborted = $true
        $abortReason = "working set reached $AbortWorkingSetMiB MiB"
    }
    elseif (($now - $started).TotalSeconds -ge $MaximumRuntimeSeconds) {
        $aborted = $true
        $abortReason = "runtime reached $MaximumRuntimeSeconds seconds"
    }
    elseif (($now - $lastProgress).TotalSeconds -ge $NoProgressTimeoutSeconds) {
        $aborted = $true
        $abortReason = "no log progress for $NoProgressTimeoutSeconds seconds"
    }
    if ($aborted) {
        Stop-Process -Id $process.Id -Force
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
    throw "phase2 surface analysis aborted: $abortReason; peak=$peakMiB MiB"
}
if (-not (Test-Path -LiteralPath (Join-Path $outputPath "model.json"))) {
    throw "phase2 surface analysis failed; peak=$peakMiB MiB; stderr=$stderrPath"
}
Write-Output "phase2_surface_analysis_completed peak_memory_mib=$peakMiB output=$outputPath"
