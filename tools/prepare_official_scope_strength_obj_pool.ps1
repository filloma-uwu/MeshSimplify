param(
    [int]$StrengthIndex = 41,
    [string]$SweepRoot = "$PSScriptRoot\..\outputs\primitive_mesh_cpp_official_scope_strengths_003",
    [string]$OriginalRoot = "$PSScriptRoot\..\test_data\real_scene\source_pool",
    [string]$OutputRoot = "$PSScriptRoot\..\outputs\real_scene\primitive_mesh_obj_s041"
)

$ErrorActionPreference = "Stop"
if ($StrengthIndex -lt 0 -or $StrengthIndex -gt 100) {
    throw "StrengthIndex must be in [0, 100]"
}
$replacementIds = [System.Collections.Generic.HashSet[int]]::new()
@(2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20) | ForEach-Object {
    [void]$replacementIds.Add($_)
}
$externalRoots = @{
    2 = (Resolve-Path "$PSScriptRoot\..\outputs\primitive_mesh_cpp_model2_sweep_batched_atomic_boxes_probe").Path
    3 = (Resolve-Path "$PSScriptRoot\..\outputs\primitive_mesh_cpp_model3_sweep_batched_atomic_boxes_probe").Path
}
$sweepPath = (Resolve-Path -LiteralPath $SweepRoot).Path
$originalPath = (Resolve-Path -LiteralPath $OriginalRoot).Path
$outputPath = [IO.Path]::GetFullPath($OutputRoot)
$strengthDirectory = "s{0:D3}" -f $StrengthIndex
New-Item -ItemType Directory -Path $outputPath -Force | Out-Null

foreach ($modelId in 1..23) {
    if (-not $replacementIds.Contains($modelId)) {
        $source = Join-Path $originalPath "$modelId.obj"
    } elseif ($externalRoots.ContainsKey($modelId)) {
        $source = Join-Path $externalRoots[$modelId] "$strengthDirectory\proxy.obj"
    } else {
        $source = Join-Path $sweepPath "models\$modelId\$strengthDirectory\proxy.obj"
    }
    if (-not (Test-Path -LiteralPath $source)) { throw "Missing model: $source" }
    Copy-Item -LiteralPath $source -Destination (Join-Path $outputPath "$modelId.obj") -Force
}
Write-Output "strength=$($StrengthIndex / 100.0)"
Write-Output "pool=$outputPath"

