param(
    [double]$Strength = 0.50,
    [string]$Executable = "$PSScriptRoot\..\build\Release\pqss-primitive-mesh-analyze.exe",
    [string]$InputRoot = "$PSScriptRoot\..\test_data\real_scene\source_pool",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\primitive_mesh_cpp_official_scope_s050",
    [int]$MaxParallel = 13
)

$ErrorActionPreference = "Stop"
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$modelIds = @(2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20)
$executablePath = (Resolve-Path -LiteralPath $Executable).Path
$inputPath = (Resolve-Path -LiteralPath $InputRoot).Path
$outputPath = [System.IO.Path]::GetFullPath($OutputRoot)
$strengthText = $Strength.ToString("0.00", $culture)

if ($Strength -lt 0.0 -or $Strength -gt 1.0) {
    throw "Strength must be in [0, 1]"
}
if ($MaxParallel -lt 1) {
    throw "MaxParallel must be positive"
}

New-Item -ItemType Directory -Path $outputPath -Force | Out-Null
$modelsPath = Join-Path $outputPath "models"
New-Item -ItemType Directory -Path $modelsPath -Force | Out-Null

$pending = [System.Collections.Generic.Queue[int]]::new()
foreach ($modelId in $modelIds) {
    $metadata = Join-Path $modelsPath "$modelId\model.json"
    if (-not (Test-Path -LiteralPath $metadata)) {
        $pending.Enqueue($modelId)
    }
}

$running = @{}
while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    while ($pending.Count -gt 0 -and $running.Count -lt $MaxParallel) {
        $modelId = $pending.Dequeue()
        $directory = Join-Path $modelsPath "$modelId"
        New-Item -ItemType Directory -Path $directory -Force | Out-Null
        $arguments = @(
            "--input", (Join-Path $inputPath "$modelId.obj"),
            "--output-dir", $directory,
            "--primitive-types", "polygon,surface",
            "--round-surface-segments", "24",
            "--analysis-strength", $strengthText
        )
        $process = Start-Process -FilePath $executablePath -ArgumentList $arguments `
            -RedirectStandardOutput (Join-Path $directory "analysis.stdout.txt") `
            -RedirectStandardError (Join-Path $directory "analysis.stderr.txt") `
            -WindowStyle Hidden -PassThru
        $running[$process.Id] = [pscustomobject]@{ Process = $process; ModelId = $modelId }
    }

    Start-Sleep -Milliseconds 200
    foreach ($processId in @($running.Keys)) {
        $job = $running[$processId]
        if (-not $job.Process.HasExited) { continue }
        $job.Process.WaitForExit()
        $metadata = Join-Path $modelsPath "$($job.ModelId)\model.json"
        if (-not (Test-Path -LiteralPath $metadata)) {
            $errorText = Get-Content (Join-Path $modelsPath "$($job.ModelId)\analysis.stderr.txt") -Raw
            throw "Model $($job.ModelId) failed: $errorText"
        }
        $running.Remove($processId)
        Write-Output "completed_model=$($job.ModelId)"
    }
}

$manifestModels = foreach ($modelId in $modelIds) {
    [ordered]@{ id = $modelId; metadata = "models/$modelId/model.json" }
}
$manifest = [ordered]@{
    algorithm = "CppHierarchicalPrimitiveMeshAnalysis"
    complete = $true
    model_count = $modelIds.Count
    analysis_strength = $Strength
    models = @($manifestModels)
}
$manifest | ConvertTo-Json -Depth 5 | Set-Content `
    -LiteralPath (Join-Path $outputPath "viewer_manifest.json") -Encoding ASCII

Write-Output "manifest=$(Join-Path $outputPath 'viewer_manifest.json')"
