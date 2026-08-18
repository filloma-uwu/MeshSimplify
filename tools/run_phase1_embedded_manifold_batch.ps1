param(
    [string]$Executable = "$PSScriptRoot\..\build_cpp_phase1\Release\pqss-topology-fill.exe",
    [string]$InputRoot = "$PSScriptRoot\..\test_data\real_scene\source_pool",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\phase1_official_uniform_embedded_manifold_v53",
    [int]$MaxParallel = 4
)

$ErrorActionPreference = "Stop"
$modelIds = @(19, 5, 16, 17, 14, 18, 3, 4, 2, 15, 13, 20, 12)
$executablePath = (Resolve-Path $Executable).Path
$inputPath = (Resolve-Path $InputRoot).Path
$outputPath = [IO.Path]::GetFullPath($OutputRoot)
if (Test-Path -LiteralPath $outputPath) {
    throw "Output root already exists: $outputPath"
}
$modelsPath = New-Item -ItemType Directory -Path (Join-Path $outputPath "models")
$pending = [Collections.Generic.Queue[int]]::new()
foreach ($modelId in $modelIds) { $pending.Enqueue($modelId) }
$running = @{}
$failures = @()

while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    while ($pending.Count -gt 0 -and $running.Count -lt $MaxParallel) {
        $modelId = $pending.Dequeue()
        $directory = Join-Path $modelsPath.FullName "$modelId"
        $process = Start-Process -FilePath $executablePath -ArgumentList @(
            "--input", (Join-Path $inputPath "$modelId.obj"),
            "--output-dir", $directory,
            "--maximum-grid-voxels", "12000000",
            "--padding", "4",
            "--maximum-steps", "32"
        ) -RedirectStandardOutput (Join-Path $outputPath "model_${modelId}.stdout.txt") `
          -RedirectStandardError (Join-Path $outputPath "model_${modelId}.stderr.txt") `
          -WindowStyle Hidden -PassThru
        $running[$process.Id] = [pscustomobject]@{
            Process = $process
            ModelId = $modelId
            PeakWorkingSetBytes = 0L
        }
        Write-Output "started_model=$modelId pid=$($process.Id)"
    }
    Start-Sleep -Milliseconds 250
    foreach ($processId in @($running.Keys)) {
        $job = $running[$processId]
        $job.Process.Refresh()
        if (-not $job.Process.HasExited) {
            $job.PeakWorkingSetBytes = [math]::Max(
                $job.PeakWorkingSetBytes, $job.Process.WorkingSet64)
            continue
        }
        $job.Process.WaitForExit()
        $job.Process.Refresh()
        $job.PeakWorkingSetBytes = [math]::Max(
            $job.PeakWorkingSetBytes, $job.Process.PeakWorkingSet64)
        $peakPath = Join-Path $outputPath "model_$($job.ModelId).peak_memory_mib.txt"
        ([double]($job.PeakWorkingSetBytes / 1MB)).ToString(
            "R", [Globalization.CultureInfo]::InvariantCulture) |
            Set-Content $peakPath -Encoding ASCII
        $metadata = Join-Path $modelsPath.FullName "$($job.ModelId)\model.json"
        if (-not (Test-Path -LiteralPath $metadata)) {
            $errorPath = Join-Path $outputPath "model_$($job.ModelId).stderr.txt"
            $errorTail = Get-Content $errorPath -Tail 1
            $failures += "model $($job.ModelId): $errorTail"
            Write-Output "failed_model=$($job.ModelId) error=$errorTail"
        }
        else {
            Write-Output "completed_model=$($job.ModelId) peak_memory_mib=$([math]::Round($job.PeakWorkingSetBytes / 1MB, 1))"
        }
        $running.Remove($processId)
    }
}

if ($failures.Count -gt 0) { throw ($failures -join [Environment]::NewLine) }
$binaryHash = (Get-FileHash $executablePath -Algorithm SHA256).Hash
$manifestModels = foreach ($modelId in ($modelIds | Sort-Object)) {
    [ordered]@{ id = $modelId; metadata = "models/$modelId/model.json" }
}
$manifest = [ordered]@{
    algorithm = "CppEmbeddedManifoldPhase1HalfedgeV53"
    binary_sha256 = $binaryHash
    complete = $true
    model_count = $manifestModels.Count
    models = @($manifestModels)
    options = [ordered]@{
        maximum_grid_voxels = 12000000
        padding = 4
        maximum_steps = 32
        target_betti = @(1, 0, 0)
        phase1_representation = "PQSSHED1 validated embedded closed halfedge"
        unique_geometric_vertices = $true
        duplicate_faces = 0
        geometric_edge_uses = 2
        cpu_threads_per_model = 1
        maximum_parallel_models = $MaxParallel
        process_memory_limit_mib = 2048
    }
}
$manifest | ConvertTo-Json -Depth 6 | Set-Content `
    (Join-Path $outputPath "viewer_manifest.json") -Encoding ASCII
Write-Output "manifest=$(Join-Path $outputPath 'viewer_manifest.json')"
