param(
    [string]$Executable = "$PSScriptRoot\..\build\Release\pqss-primitive-mesh-analyze.exe",
    [string]$InputRoot = "$PSScriptRoot\..\test_data\real_scene\source_pool",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\primitive_mesh_cpp_uniform_official_scope_001",
    [int]$MaxParallel = 1
)

$ErrorActionPreference = "Stop"
$modelIds = @(2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20)
$executablePath = (Resolve-Path $Executable).Path
$inputPath = (Resolve-Path $InputRoot).Path
$outputPath = [IO.Path]::GetFullPath($OutputRoot)
$modelsPath = Join-Path $outputPath "models"
New-Item -ItemType Directory -Path $modelsPath -Force | Out-Null

$pending = [Collections.Generic.Queue[int]]::new()
foreach ($modelId in $modelIds) {
    if (-not (Test-Path (Join-Path $modelsPath "$modelId\model.json"))) {
        $pending.Enqueue($modelId)
    }
}
$running = @{}
$failures = @()
while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    while ($pending.Count -gt 0 -and $running.Count -lt $MaxParallel) {
        $modelId = $pending.Dequeue()
        $directory = Join-Path $modelsPath "$modelId"
        New-Item -ItemType Directory -Path $directory -Force | Out-Null
        $process = Start-Process -FilePath $executablePath -ArgumentList @(
            "--input", (Join-Path $inputPath "$modelId.obj"),
            "--output-dir", $directory,
            "--primitive-types", "polygon,surface",
            "--round-surface-segments", "24",
            "--maximum-process-memory-gb", "2"
        ) -RedirectStandardOutput (Join-Path $directory "analysis.stdout.txt") `
          -RedirectStandardError (Join-Path $directory "analysis.stderr.txt") `
          -WindowStyle Hidden -PassThru
        $running[$process.Id] = [pscustomobject]@{
            Process = $process; ModelId = $modelId; Directory = $directory
            PeakWorkingSetBytes = 0L; PeakPrivateBytes = 0L
        }
    }
    Start-Sleep -Milliseconds 200
    foreach ($processId in @($running.Keys)) {
        $job = $running[$processId]
        $job.Process.Refresh()
        if (-not $job.Process.HasExited) {
            $job.PeakWorkingSetBytes = [math]::Max(
                $job.PeakWorkingSetBytes, $job.Process.WorkingSet64)
            $job.PeakPrivateBytes = [math]::Max(
                $job.PeakPrivateBytes, $job.Process.PrivateMemorySize64)
            continue
        }
        $job.Process.WaitForExit()
        $job.Process.Refresh()
        $job.PeakWorkingSetBytes = [math]::Max(
            $job.PeakWorkingSetBytes, $job.Process.PeakWorkingSet64)
        ($job.PeakWorkingSetBytes / 1MB).ToString(
            "R", [Globalization.CultureInfo]::InvariantCulture) | Set-Content `
            (Join-Path $job.Directory "peak_memory_mb.txt") -Encoding ASCII
        ($job.PeakPrivateBytes / 1MB).ToString(
            "R", [Globalization.CultureInfo]::InvariantCulture) | Set-Content `
            (Join-Path $job.Directory "peak_private_memory_mb.txt") -Encoding ASCII
        if (-not (Test-Path (Join-Path $job.Directory "model.json"))) {
            $errorText = Get-Content (Join-Path $job.Directory "analysis.stderr.txt") -Raw
            $failures += "Model $($job.ModelId) failed: $errorText"
            Write-Error -ErrorAction Continue $failures[-1]
        }
        else {
            Write-Output "completed_model=$($job.ModelId)"
        }
        $running.Remove($processId)
    }
}

if ($failures.Count -gt 0) {
    throw ($failures -join [Environment]::NewLine)
}

$models = foreach ($modelId in $modelIds) {
    [ordered]@{ id = $modelId; metadata = "models/$modelId/model.json" }
}
$manifest = [ordered]@{
    algorithm = "CppStreamingAdjacentSurfaceFixedPointV6"
    complete = $true
    model_count = $modelIds.Count
    models = @($models)
}
$manifest | ConvertTo-Json -Depth 5 | Set-Content `
    (Join-Path $outputPath "viewer_manifest.json") -Encoding ASCII
Write-Output "manifest=$(Join-Path $outputPath 'viewer_manifest.json')"
