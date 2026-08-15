param(
    [Parameter(Mandatory = $true)]
    [string]$OutputDirectory
)

$ErrorActionPreference = 'Stop'
$outputDirectory = (Resolve-Path -LiteralPath $OutputDirectory).Path
$logPath = Join-Path $outputDirectory 'pdf_export.log'
Set-Content -LiteralPath $logPath -Value ("PDF export started: " + (Get-Date -Format o)) -Encoding UTF8

$writer = New-Object -ComObject KWPS.Application
$writer.Visible = $false
$writer.DisplayAlerts = 0
$fixedFormatClass = $null

try {
    foreach ($language in @('EN', 'ZH')) {
        $docxPath = Join-Path $outputDirectory ("IVDD_cohort_aware_manuscript_{0}.docx" -f $language)
        $temporaryPdf = Join-Path $outputDirectory ("IVDD_cohort_aware_manuscript_{0}.exporting.pdf" -f $language)
        $finalPdf = Join-Path $outputDirectory ("IVDD_cohort_aware_manuscript_{0}.pdf" -f $language)
        if (Test-Path -LiteralPath $temporaryPdf) {
            Remove-Item -LiteralPath $temporaryPdf -Force
        }

        $document = $writer.Documents.Open($docxPath, $false, $true)
        try {
            $pageCount = $document.ComputeStatistics(2)
            Add-Content -LiteralPath $logPath -Value ("{0}: {1} pages" -f $language, $pageCount) -Encoding UTF8
            $document.ExportAsFixedFormat(
                $temporaryPdf, 17, $false, 0, 0, 1, 1, 0,
                $true, $true, 0, $true, $true, $false,
                $fixedFormatClass
            )
        }
        finally {
            $document.Close(0)
            [void][Runtime.InteropServices.Marshal]::ReleaseComObject($document)
        }

        if (-not (Test-Path -LiteralPath $temporaryPdf)) {
            throw "Word did not create $temporaryPdf"
        }
        Move-Item -LiteralPath $temporaryPdf -Destination $finalPdf -Force
        Add-Content -LiteralPath $logPath -Value ("{0}: exported" -f $language) -Encoding UTF8
    }
    Add-Content -LiteralPath $logPath -Value ("PDF export completed: " + (Get-Date -Format o)) -Encoding UTF8
}
catch {
    Add-Content -LiteralPath $logPath -Value ("PDF export failed: " + $_.Exception.ToString()) -Encoding UTF8
    throw
}
finally {
    $writer.Quit()
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($writer)
}
