param(
    [string]$Executable = "$PSScriptRoot\..\build\Release\pqss-primitive-mesh-analyze.exe",
    [string]$InputRoot = "$PSScriptRoot\..\test_data\real_scene\source_pool",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\surface_primitives_uniform",
    [int[]]$ModelIds = @(2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20),
    [int]$MaxParallel = 1,
    [double]$MaxProcessMemoryGB = 2.0,
    [double]$MinimumFreeMemoryGB = 8.0,
    [string]$Algorithm = "UnifiedSurfacePrimitiveOptimizer"
)

$ErrorActionPreference = "Stop"
if ($MaxParallel -lt 1) { throw "MaxParallel must be positive" }
if ($MaxProcessMemoryGB -le 0.0 -or $MaxProcessMemoryGB -gt 2.0 -or
    $MinimumFreeMemoryGB -le 0.0) {
    throw "process memory must be in (0, 2] GB and free-memory safety must be positive"
}

$executablePath = (Resolve-Path -LiteralPath $Executable).Path
$inputPath = (Resolve-Path -LiteralPath $InputRoot).Path
$outputPath = [System.IO.Path]::GetFullPath($OutputRoot)
New-Item -ItemType Directory -Path $outputPath -Force | Out-Null

$pending = [System.Collections.Generic.Queue[int]]::new()
foreach ($modelId in $ModelIds) {
    $metadata = Join-Path $outputPath "$modelId\model.json"
    if (-not (Test-Path -LiteralPath $metadata)) { $pending.Enqueue($modelId) }
}

$running = @{}
$nextSystemMemoryCheck = Get-Date
while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    while ($pending.Count -gt 0 -and $running.Count -lt $MaxParallel) {
        $modelId = $pending.Dequeue()
        $directory = Join-Path $outputPath $modelId
        New-Item -ItemType Directory -Path $directory -Force | Out-Null
        $arguments = @(
            "--input", (Join-Path $inputPath "$modelId.obj"),
            "--output-dir", $directory,
            "--unified-candidates",
            "--maximum-process-memory-gb", $MaxProcessMemoryGB
        )
        $process = Start-Process -FilePath $executablePath -ArgumentList $arguments `
            -RedirectStandardOutput (Join-Path $directory "analysis.stdout.txt") `
            -RedirectStandardError (Join-Path $directory "analysis.stderr.txt") `
            -WindowStyle Hidden -PassThru
        $running[$process.Id] = [pscustomobject]@{
            Process = $process
            ModelId = $modelId
            PeakWorkingSetBytes = [int64]0
            PeakPrivateBytes = [int64]0
        }
    }

    Start-Sleep -Milliseconds 200
    foreach ($processId in @($running.Keys)) {
        $job = $running[$processId]
        $job.Process.Refresh()
        $job.PeakWorkingSetBytes = [math]::Max(
            $job.PeakWorkingSetBytes, $job.Process.WorkingSet64)
        $job.PeakPrivateBytes = [math]::Max(
            $job.PeakPrivateBytes, $job.Process.PrivateMemorySize64)
        if ($job.PeakWorkingSetBytes -gt $MaxProcessMemoryGB * 1GB -or
            $job.PeakPrivateBytes -gt $MaxProcessMemoryGB * 1GB) {
            foreach ($active in $running.Values) {
                if (-not $active.Process.HasExited) { $active.Process.Kill() }
            }
            throw "Model $($job.ModelId) exceeded the $MaxProcessMemoryGB GB process-memory limit (working set or private bytes)"
        }
        if (-not $job.Process.HasExited) { continue }
        $job.Process.WaitForExit()
        $peakMemoryMB = [math]::Round($job.PeakWorkingSetBytes / 1MB, 1)
        $peakPrivateMB = [math]::Round($job.PeakPrivateBytes / 1MB, 1)
        Set-Content -LiteralPath (Join-Path $outputPath "$($job.ModelId)\peak_memory_mb.txt") `
            -Value $peakMemoryMB -Encoding ASCII
        Set-Content -LiteralPath (Join-Path $outputPath "$($job.ModelId)\peak_private_memory_mb.txt") `
            -Value $peakPrivateMB -Encoding ASCII
        $metadata = Join-Path $outputPath "$($job.ModelId)\model.json"
        if (-not (Test-Path -LiteralPath $metadata)) {
            $errorPath = Join-Path $outputPath "$($job.ModelId)\analysis.stderr.txt"
            $errorText = if (Test-Path -LiteralPath $errorPath) {
                Get-Content -LiteralPath $errorPath -Raw
            } else { "no stderr was produced" }
            throw "Model $($job.ModelId) failed: $errorText"
        }
        $running.Remove($processId)
        Write-Output "completed_model=$($job.ModelId) peak_working_set_mb=$peakMemoryMB peak_private_mb=$peakPrivateMB"
    }

    if ((Get-Date) -ge $nextSystemMemoryCheck) {
        $operatingSystem = Get-CimInstance Win32_OperatingSystem
        $freeMemoryGB = $operatingSystem.FreePhysicalMemory / 1MB
        if ($freeMemoryGB -lt $MinimumFreeMemoryGB) {
            foreach ($active in $running.Values) {
                if (-not $active.Process.HasExited) { $active.Process.Kill() }
            }
            throw "System free memory fell below $MinimumFreeMemoryGB GB"
        }
        $nextSystemMemoryCheck = (Get-Date).AddSeconds(2)
    }
}

$manifestModels = foreach ($modelId in $ModelIds) {
    [ordered]@{ id = $modelId; metadata = "$modelId/model.json" }
}
$manifest = [ordered]@{
    algorithm = $Algorithm
    complete = $true
    note = "One Release binary and identical unified-candidate options for every model."
    model_count = $ModelIds.Count
    models = @($manifestModels)
}
$manifest | ConvertTo-Json -Depth 5 | Set-Content `
    -LiteralPath (Join-Path $outputPath "viewer_manifest.json") -Encoding ASCII

Write-Output "manifest=$(Join-Path $outputPath 'viewer_manifest.json')"
