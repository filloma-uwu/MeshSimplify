param(
    [int]$ModelId = 2,
    [string]$Executable = "$PSScriptRoot\..\build\Release\pqss-primitive-mesh-analyze.exe",
    [string]$InputObj = "$PSScriptRoot\..\test_data\real_scene\source_pool\2.obj",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\primitive_mesh_cpp_model2_strengths_002",
    [int]$MaxParallel = 10
)

$ErrorActionPreference = "Stop"
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$executablePath = (Resolve-Path -LiteralPath $Executable).Path
$inputPath = (Resolve-Path -LiteralPath $InputObj).Path
$outputPath = [System.IO.Path]::GetFullPath($OutputRoot)

if ($MaxParallel -lt 1) {
    throw "MaxParallel must be positive"
}

New-Item -ItemType Directory -Path $outputPath -Force | Out-Null
$pending = [System.Collections.Generic.Queue[int]]::new()
0..100 | ForEach-Object {
    $modelJson = Join-Path $outputPath ("s{0:D3}\model.json" -f $_)
    if (-not (Test-Path -LiteralPath $modelJson)) {
        $pending.Enqueue($_)
    }
}
$running = @{}
$completed = 101 - $pending.Count
$startedAt = Get-Date

while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    while ($pending.Count -gt 0 -and $running.Count -lt $MaxParallel) {
        $index = $pending.Dequeue()
        $strength = ($index / 100.0).ToString("0.00", $culture)
        $name = "s{0:D3}" -f $index
        $directory = Join-Path $outputPath $name
        New-Item -ItemType Directory -Path $directory -Force | Out-Null
        $stdout = Join-Path $directory "analysis.stdout.txt"
        $stderr = Join-Path $directory "analysis.stderr.txt"
        $arguments = @(
            "--input", $inputPath,
            "--output-dir", $directory,
            "--primitive-types", "polygon,frustum",
            "--cone-segments", "24",
            "--analysis-strength", $strength
        )
        $process = Start-Process -FilePath $executablePath -ArgumentList $arguments `
            -RedirectStandardOutput $stdout -RedirectStandardError $stderr `
            -WindowStyle Hidden -PassThru
        $running[$process.Id] = [pscustomobject]@{
            Process = $process
            Index = $index
            Name = $name
            Stderr = $stderr
        }
    }

    Start-Sleep -Milliseconds 250
    foreach ($id in @($running.Keys)) {
        $job = $running[$id]
        if (-not $job.Process.HasExited) {
            continue
        }
        $job.Process.WaitForExit()
        $modelJson = Join-Path $outputPath "$($job.Name)\model.json"
        if (-not (Test-Path -LiteralPath $modelJson)) {
            $message = Get-Content -LiteralPath $job.Stderr -Raw
            throw "Strength index $($job.Index) did not produce model.json: $message"
        }
        $running.Remove($id)
        $completed++
        $elapsed = [math]::Round(((Get-Date) - $startedAt).TotalSeconds, 1)
        Write-Output "completed=$completed/101 strength=$([math]::Round($job.Index / 100.0, 2)) elapsed_seconds=$elapsed"
    }
}

foreach ($index in 0..100) {
    $modelJson = Join-Path $outputPath ("s{0:D3}\model.json" -f $index)
    if (-not (Test-Path -LiteralPath $modelJson)) {
        throw "Sweep is incomplete: missing $modelJson"
    }
}

$variants = foreach ($index in 0..100) {
    [ordered]@{
        model_id = $ModelId
        strength = $index / 100.0
        metadata = "s{0:D3}/model.json" -f $index
    }
}
$manifest = [ordered]@{
    algorithm = "CppHierarchicalPrimitiveMeshAnalysis"
    complete = $true
    model_count = 1
    models = @([ordered]@{ id = $ModelId; metadata = "s050/model.json" })
    strength_variants = @($variants)
}
$manifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $outputPath "viewer_manifest.json") -Encoding ASCII
Write-Output "manifest=$(Join-Path $outputPath 'viewer_manifest.json')"
