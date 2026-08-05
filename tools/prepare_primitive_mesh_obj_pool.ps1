param(
    [string]$AnalysisRoot = "$PSScriptRoot\..\outputs\primitive_mesh_cpp_official_scope_s050",
    [string]$OriginalRoot = "$PSScriptRoot\..\test_data\real_scene\source_pool",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\real_scene\primitive_mesh_obj_s050"
)

$ErrorActionPreference = "Stop"
$replacementIds = [System.Collections.Generic.HashSet[int]]::new()
@(2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20) | ForEach-Object {
    [void]$replacementIds.Add($_)
}
$analysisPath = (Resolve-Path -LiteralPath $AnalysisRoot).Path
$originalPath = (Resolve-Path -LiteralPath $OriginalRoot).Path
$outputPath = [System.IO.Path]::GetFullPath($OutputRoot)
New-Item -ItemType Directory -Path $outputPath -Force | Out-Null

foreach ($modelId in 1..23) {
    if ($replacementIds.Contains($modelId)) {
        $source = Join-Path $analysisPath "models\$modelId\proxy.obj"
    } else {
        $source = Join-Path $originalPath "$modelId.obj"
    }
    if (-not (Test-Path -LiteralPath $source)) { throw "Missing model: $source" }
    Copy-Item -LiteralPath $source -Destination (Join-Path $outputPath "$modelId.obj") -Force
}

Write-Output "pool=$outputPath"
