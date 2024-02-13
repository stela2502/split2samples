# Set the directory path for storing Rustody files in Windows
$rustodyFilesPath = "$env:LOCALAPPDATA\Rustody\files"
$rustodyBinPath = "$env:LOCALAPPDATA\Rustody\bin"

# Check if RustodyFiles environment variable is set in user profile

$profilePath = $PROFILE.CurrentUserAllHosts
if (!(Test-Path $profilePath)) {
    New-Item -ItemType File -Path $profilePath -Force
}

if (-not (Test-Path -Path $PROFILE)) {
    New-Item -Type File -Path $PROFILE -Force
}

if (-not (Select-String -Path $PROFILE -Pattern "RustodyFiles")) {
    # Add RustodyFiles environment variable to user profile
    Add-Content -Path $PROFILE -Value "Set-Item -Path Env:RustodyFiles -Value '$rustodyFilesPath'"
}

# Add Rustody bin directory to the user's PATH environment variable in the registry
# Check if the directory is already in the PATH
if (!$env:PATH.ToLower().Contains($rustodyBinPath.ToLower())) {
    # Add Rustody bin directory to the user's profile script
    $profilePath = $PROFILE.CurrentUserAllHosts
    if (!(Test-Path $profilePath)) {
        New-Item -ItemType File -Path $profilePath -Force
    }
    Add-Content -Path $profilePath -Value '$env:PATH += ";$rustodyBinPath"'
    $env:PATH += ";$rustodyBinPath"
} else {
    Write-Host "Rustody bin directory is already in the PATH."
}


# activate the variables in this instance, too
$env:RustodyFiles = $rustodyFilesPath

# Create directory and copy files if not already exists
New-Item -ItemType Directory -Path $rustodyFilesPath -Force
Copy-Item -Path ./resources/CellRanger/*.gz -Destination $rustodyFilesPath -Force

New-Item -ItemType Directory -Path $rustodyBinPath -Force
Copy-Item -Path ./target/release/*.exe -Destination $rustodyBinPath -Force


