# TODO create an array of arguments so we can backup multiple env variables and create text file names dynamically. They are currently hardcoded for our scenario.
function Backup-EnvironmentVariables($userEnv, $systemEnv) {
    $backupDir = "C:\Users\$Env:UserName\Desktop\Environment-Variable-Backup"

    if (-not (Test-Path -Path $backupDir)) {
        New-Item -ItemType Directory -Path $backupDir
    }
        
    # Output the environment variable data to the files
    $machinePath | Out-File -FilePath "$backupDir\machine-env-backup.txt"
    $userPath | Out-File -FilePath "$backupDir\user-env-backup.txt"
    

}
function Install-BuildToolsAndLibraries {
    # Install swig and set environment variable
    $swigUrl = "https://drive.google.com/uc?id=137YbZdSDL8YAjvbWYfdTBj244pgmPcmN&export=download"
    $swigInstallationPath = "C:\swigwin-4.1.1"
    Write-Output "Downloading SWIG..."
    Invoke-WebRequest -Uri $swigUrl -OutFile "$env:TEMP\\swigwin-4.1.1.zip"
    Expand-Archive -Path "$env:TEMP\\swigwin-4.1.1.zip" -DestinationPath "C:\\" -Force
    [System.Environment]::SetEnvironmentVariable("Path", $env:Path + ";$swigInstallationPath", [System.EnvironmentVariableTarget]::Machine)
    # Install VS Build Tools    
    $buildToolsArgumentList = "--add Microsoft.VisualStudio.Workload.VCTools;includeRecommended --passive"
    Write-Output "Downloading VS Build Tools..."
    Invoke-WebRequest -Uri "https://aka.ms/vs/17/release/vs_BuildTools.exe" -OutFile "$env:TEMP\\vs_buildtools.exe"
    Start-Process -FilePath "$env:TEMP\\vs_buildtools.exe" -ArgumentList $buildToolsArgumentList -Wait
}
function Install-RequiredPythonPackages {
    py -m ensurepip --upgrade
    py -m pip install --upgrade pip
    py -m pip install --upgrade setuptools
    py -m pip install SimpleItk
}
$pythonDownloadUrl = "https://www.python.org/ftp/python/3.10.0/python-3.10.0-amd64.exe"
function Install-PyPore3d {
    # Refresh path environment variable for current session, to be able to use "py" command
    $env:Path = [System.Environment]::GetEnvironmentVariable("Path", [System.EnvironmentVariableTarget]::Machine)
    $pythonPath = (Split-Path (Get-Command python).Path)
    $pythonSitePackagesPath = Join-Path $pythonPath "\Lib\site-packages"
    $pyPore3dUrl = "https://gitlab.elettra.eu/aboulhassan.amal/PyPore3D/-/wikis/uploads/22d2cb768bc4934a88685f80810470f8/PyPore3D_Win.zip"
    
    Write-Output "Downloading PyPore3D..."
    
    Invoke-WebRequest -Uri $pyPore3dUrl -OutFile "$env:TEMP\\PyPore3D_Win.zip"
    Expand-Archive -Path "$env:TEMP\\PyPore3D_Win.zip" -DestinationPath "C:\\Program Files" -Force    
    Set-Location "C:\Program Files\PyPore3D_Win"
    py setup.py build_ext --inplace

    Move-Item -Path "C:\Program Files\PyPore3D_Win\pypore3d" -Destination $pythonSitePackagesPath -Force
    Read-Host -Prompt "Press Enter to continue" 
}

$IsAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltinRole] "Administrator")

if (-NOT $IsAdmin) {
    # This will automatically show windows default prompt for admin privileges, no need to add CLI prompt
    Start-Process -FilePath "powershell.exe" -ArgumentList "-File `"$PSCommandPath`" -NoNewWindow" -Verb RunAs
}
else {
    if (Get-Command 'py' -ErrorAction SilentlyContinue) {
        Write-Output "Launcher found. Installing PyPore3D..."
    }
    else {
        # Check if 'python' command exists
        if (Get-Command 'python' -ErrorAction SilentlyContinue) {
            $pythonCheck = & python --version 2>&1
            if ($pythonCheck -notmatch 'Microsoft Store') {
                $pythonVersion = (python --version).Split(" ")[1]
                Write-Output "Python is installed. Current version:$pythonVersion. `nPython launcher is not found `nDownloading and installing py launcher..."
                Invoke-WebRequest -Uri $pythonDownloadUrl -OutFile "$env:TEMP\\python-3.12.0-amd64.exe"
                Start-Process -FilePath "$env:TEMP\\python-3.12.0-amd64.exe" -ArgumentList "/passive Include_launcher=1 PrependPath=1" -Wait
            }
            else {
                Write-Output "Python and launcher not found. Installing Python and launcher..."
                Invoke-WebRequest -Uri $pythonDownloadUrl -OutFile "$env:TEMP\\python-3.12.0-amd64.exe"
                Start-Process -FilePath "$env:TEMP\\python-3.12.0-amd64.exe" -ArgumentList "/passive PrependPath=1" -Wait
            }
            # Could've gone for this condition directly, however the extra if statement is added in order to handle the case of missing python path in the future.
        }
        else {
            Write-Output "Python and launcher not found. Installing Python and launcher..."
            Invoke-WebRequest -Uri $pythonDownloadUrl -OutFile "$env:TEMP\\python-3.12.0-amd64.exe"
            Start-Process -FilePath "$env:TEMP\\python-3.12.0-amd64.exe" -ArgumentList "/passive PrependPath=1" -Wait
        }
    }
    Install-BuildToolsAndLibraries
    Install-RequiredPythonPackages
    Install-PyPore3d
}
