<#
Run with PowerShell's execution-policy bypass when the local policy blocks
scripts, for example:
  powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\download_geo.ps1 -Accession GSE229711 -IncludeSoft
#>

param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('GSE229711', 'GSE230808', 'GSE165722', 'GSE244889', 'GSE251686', 'GSE153066', 'GSE160756')]
    [string]$Accession,
    [switch]$IncludeSoft
)

$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'

$projectRoot = Split-Path -Parent $PSScriptRoot
$rawRoot = Join-Path $projectRoot 'data\raw'
$logsRoot = Join-Path $projectRoot 'logs'
New-Item -ItemType Directory -Force -Path $rawRoot, $logsRoot | Out-Null

$seriesPrefix = $Accession.Substring(0, 6) + 'nnn'
$baseUri = "https://ftp.ncbi.nlm.nih.gov/geo/series/$seriesPrefix/$Accession"
# Most selected series publish a single RAW TAR, but GSE153066 publishes its
# usable public single-cell matrix as one compressed supplementary TSV.  Keep
# the source asset name unchanged and do not substitute SRA data.
if ($Accession -eq 'GSE153066') {
    $primaryAssets = @(
        [pscustomobject]@{
            asset = 'GSE153066_AllSample.counts.tsv.gz'
            url = "$baseUri/suppl/GSE153066_AllSample.counts.tsv.gz"
        }
    )
}
elseif ($Accession -eq 'GSE244889') {
    # This series publishes raw 10x matrices plus a separate processed bulk
    # FPKM table.  Retain both assets with their original GEO names, but do
    # not imply that the FPKM table is suitable for count-based scRNA models.
    $primaryAssets = @(
        [pscustomobject]@{
            asset = 'GSE244889_RAW.tar'
            url = "$baseUri/suppl/GSE244889_RAW.tar"
        },
        [pscustomobject]@{
            asset = 'GSE244889_FPKM.txt.gz'
            url = "$baseUri/suppl/GSE244889_FPKM.txt.gz"
        },
        [pscustomobject]@{
            asset = 'GSE244889_filelist.txt'
            url = "$baseUri/suppl/filelist.txt"
        }
    )
}
else {
    $primaryAssets = @(
        [pscustomobject]@{
            asset = "${Accession}_RAW.tar"
            url = "$baseUri/suppl/${Accession}_RAW.tar"
        }
    )
}

function Get-FileSha256([string]$Path) {
    (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}

function Download-GeoFile([string]$Uri, [string]$Destination) {
    $partial = "$Destination.partial"
    if (Test-Path -LiteralPath $Destination) {
        return
    }

    $head = [Net.HttpWebRequest]::Create($Uri)
    $head.Method = 'HEAD'
    $head.UserAgent = 'IVDD-cross-cohort-reproducibility/1.0'
    $headResponse = $head.GetResponse()
    try {
        $expectedBytes = [int64]$headResponse.ContentLength
    }
    finally {
        $headResponse.Dispose()
    }

    if ($expectedBytes -le 0) {
        throw "The server did not provide a valid file size for $Uri"
    }

    $existingBytes = if (Test-Path -LiteralPath $partial) {
        (Get-Item -LiteralPath $partial).Length
    }
    else {
        [int64]0
    }

    if ($existingBytes -gt $expectedBytes) {
        throw "Partial file is larger than the expected server file: $partial"
    }

    $chunkBytes = 8MB
    while ($existingBytes -lt $expectedBytes) {
        $request = [Net.HttpWebRequest]::Create($Uri)
        $request.Method = 'GET'
        $request.UserAgent = 'IVDD-cross-cohort-reproducibility/1.0'
        $request.Timeout = 120000
        $request.ReadWriteTimeout = 120000
        if ($existingBytes -gt 0) {
            $request.AddRange($existingBytes)
        }

        $response = $request.GetResponse()
        try {
            if ($existingBytes -gt 0 -and [int]$response.StatusCode -ne 206) {
                throw "Server did not honor the resume range for $Uri"
            }

            $stream = $response.GetResponseStream()
            $file = [IO.File]::Open($partial, [IO.FileMode]::Append, [IO.FileAccess]::Write, [IO.FileShare]::None)
            try {
                $buffer = New-Object byte[] 1048576
                $downloadedThisRequest = [int64]0
                while ($downloadedThisRequest -lt $chunkBytes -and $existingBytes -lt $expectedBytes) {
                    $remaining = [int][Math]::Min($buffer.Length, $chunkBytes - $downloadedThisRequest)
                    $count = $stream.Read($buffer, 0, $remaining)
                    if ($count -le 0) {
                        break
                    }
                    $file.Write($buffer, 0, $count)
                    $downloadedThisRequest += $count
                    $existingBytes += $count
                }
            }
            finally {
                $file.Dispose()
                $stream.Dispose()
            }
        }
        finally {
            $response.Dispose()
        }

        Write-Host ("{0}: {1:N1}% ({2:N0}/{3:N0} bytes)" -f $Destination, (100 * $existingBytes / $expectedBytes), $existingBytes, $expectedBytes)
    }

    if ((Get-Item -LiteralPath $partial).Length -ne $expectedBytes) {
        throw "Completed byte count does not match expected file size for $Uri"
    }

    Move-Item -LiteralPath $partial -Destination $Destination -Force
}

$records = @()
foreach ($asset in $primaryAssets) {
    $assetPath = Join-Path $rawRoot $asset.asset
    Download-GeoFile -Uri $asset.url -Destination $assetPath
    $records += [pscustomobject]@{
        accession = $Accession
        asset = $asset.asset
        url = $asset.url
        retrieved_at_utc = (Get-Date).ToUniversalTime().ToString('o')
        bytes = (Get-Item -LiteralPath $assetPath).Length
        sha256 = Get-FileSha256 -Path $assetPath
    }
}

if ($IncludeSoft) {
    $softName = "${Accession}_family.soft.gz"
    $softUri = "$baseUri/soft/$softName"
    $softPath = Join-Path $rawRoot $softName
    Download-GeoFile -Uri $softUri -Destination $softPath
    $records += [pscustomobject]@{
        accession = $Accession
        asset = $softName
        url = $softUri
        retrieved_at_utc = (Get-Date).ToUniversalTime().ToString('o')
        bytes = (Get-Item -LiteralPath $softPath).Length
        sha256 = Get-FileSha256 -Path $softPath
    }
}

$manifestPath = Join-Path $rawRoot 'download_manifest.ndjson'
$records | ForEach-Object { $_ | ConvertTo-Json -Compress } | Add-Content -LiteralPath $manifestPath -Encoding utf8
$records | Format-Table -AutoSize
