param(
    [string]$Executable = "$PSScriptRoot\..\build\Release\pqss-primitive-mesh-analyze.exe",
    [string]$InputRoot = "$PSScriptRoot\..\test_data\real_scene\source_pool",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\primitive_mesh_cpp_official_scope_strengths_003",
    [int]$MaxParallelModels = 11
)

$ErrorActionPreference = "Stop"
$modelIds = @(2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20)
$generatedIds = @(4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20)
$executablePath = (Resolve-Path -LiteralPath $Executable).Path
$inputPath = (Resolve-Path -LiteralPath $InputRoot).Path
$outputPath = [System.IO.Path]::GetFullPath($OutputRoot)
$modelsPath = Join-Path $outputPath "models"
New-Item -ItemType Directory -Path $modelsPath -Force | Out-Null

if ($MaxParallelModels -lt 1) { throw "MaxParallelModels must be positive" }
$pending = [System.Collections.Generic.Queue[int]]::new()
foreach ($modelId in $generatedIds) {
    if (-not (Test-Path (Join-Path $modelsPath "$modelId\viewer_manifest.json"))) {
        $pending.Enqueue($modelId)
    }
}
$running = @{}
$startedAt = Get-Date
while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    while ($pending.Count -gt 0 -and $running.Count -lt $MaxParallelModels) {
        $modelId = $pending.Dequeue()
        $directory = Join-Path $modelsPath "$modelId"
        New-Item -ItemType Directory -Path $directory -Force | Out-Null
        $process = Start-Process -FilePath $executablePath -ArgumentList @(
            "--input", (Join-Path $inputPath "$modelId.obj"),
            "--output-dir", $directory,
            "--primitive-types", "polygon,frustum",
            "--cone-segments", "24",
            "--strength-sweep"
        ) -RedirectStandardOutput (Join-Path $directory "sweep.stdout.txt") `
          -RedirectStandardError (Join-Path $directory "sweep.stderr.txt") `
          -WindowStyle Hidden -PassThru
        $running[$process.Id] = [pscustomobject]@{
            Process = $process; ModelId = $modelId; Directory = $directory
        }
    }
    Start-Sleep -Milliseconds 250
    foreach ($processId in @($running.Keys)) {
        $job = $running[$processId]
        if (-not $job.Process.HasExited) { continue }
        $job.Process.WaitForExit()
        if (-not (Test-Path (Join-Path $job.Directory "viewer_manifest.json"))) {
            $errorText = Get-Content (Join-Path $job.Directory "sweep.stderr.txt") -Raw
            throw "Model $($job.ModelId) failed: $errorText"
        }
        $running.Remove($processId)
        $elapsed = [math]::Round(((Get-Date) - $startedAt).TotalSeconds, 1)
        Write-Output "completed_model=$($job.ModelId) elapsed_seconds=$elapsed"
    }
}

$externalRoots = @{
    2 = "/outputs/primitive_mesh_cpp_model2_sweep_batched_atomic_boxes_probe"
    3 = "/outputs/primitive_mesh_cpp_model3_sweep_batched_atomic_boxes_probe"
}
$variants = foreach ($modelId in $modelIds) {
    foreach ($index in 0..100) {
        if ($externalRoots.ContainsKey($modelId)) {
            $metadata = "$($externalRoots[$modelId])/s$($index.ToString('D3'))/model.json"
        } else {
            $metadata = "models/$modelId/s$($index.ToString('D3'))/model.json"
        }
        [ordered]@{ model_id = $modelId; strength = $index / 100.0; metadata = $metadata }
    }
}
$models = foreach ($modelId in $modelIds) {
    if ($externalRoots.ContainsKey($modelId)) {
        $metadata = "$($externalRoots[$modelId])/s040/model.json"
    } else {
        $metadata = "models/$modelId/s040/model.json"
    }
    [ordered]@{ id = $modelId; metadata = $metadata }
}
$manifest = [ordered]@{
    algorithm = "CppHierarchicalPrimitiveMeshAnalysis"
    complete = $true
    model_count = $modelIds.Count
    default_strength = 0.41
    models = @($models)
    strength_variants = @($variants)
}
$manifest | ConvertTo-Json -Depth 5 | Set-Content `
    -LiteralPath (Join-Path $outputPath "viewer_manifest.json") -Encoding ASCII
Write-Output "manifest=$(Join-Path $outputPath 'viewer_manifest.json')"
