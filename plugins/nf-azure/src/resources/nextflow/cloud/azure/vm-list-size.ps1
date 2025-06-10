$scriptPath = Split-Path -Parent $MyInvocation.MyCommand.Path

# Get all supported locations for VM sizes using the compute provider
$supportedLocations = az provider show --namespace Microsoft.Compute `
    --query "resourceTypes[?resourceType=='virtualMachines'].locations[]" | ConvertFrom-Json

Write-Host "Found $($supportedLocations.Count) supported locations"

foreach ($region in $supportedLocations) {
    $normalizedRegion = $region.Replace(' ', '').ToLower()
    Write-Host "Listing VM size for $normalizedRegion"
    try {
        $output = az vm list-sizes -l $normalizedRegion
        if ($output) {
            $output | Out-File -FilePath "$scriptPath\vm-list-size-$normalizedRegion.json"
            Write-Host "Successfully saved VM sizes for $normalizedRegion"
        }
    }
    catch {
        Write-Warning "Failed to get VM sizes for ${normalizedRegion}: $_"
    }
}