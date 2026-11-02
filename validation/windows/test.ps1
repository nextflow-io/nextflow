#
#  Copyright 2013-2026, Seqera Labs
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

# End-to-end test of Nextflow running natively on Windows with Docker.
#
# Usage: test.ps1 -Jar <path to nextflow one jar> [-WorkDir <launch dir>]
#
# Nextflow is started with the `nextflow.ps1` launcher using Windows PowerShell,
# and the `docker` command must run Linux containers (Docker Desktop, or the
# docker shim used by the CI on GitHub-hosted runners).

param(
    [Parameter(Mandatory = $true)][string]$Jar,
    [string]$WorkDir = (Join-Path ([Environment]::GetFolderPath('UserProfile')) 'nxf-windows-test')
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$launcher = Join-Path (Resolve-Path (Join-Path $root '..\..')) 'nextflow.ps1'
$pipeline = Join-Path $root 'main.nf'
$env:NXF_BIN = (Resolve-Path $Jar).Path
$failures = 0

function Check([bool]$condition, [string]$message) {
    if( $condition ) {
        Write-Host "OK:   $message"
    }
    else {
        Write-Host "FAIL: $message"
        $script:failures++
    }
}

function Invoke-Nextflow([string[]]$arguments) {
    Write-Host "`n> nextflow $($arguments -join ' ')"
    $ErrorActionPreference = 'Continue'
    $output = & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $launcher @arguments 2>&1 | ForEach-Object { "$_" }
    $code = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    $output | ForEach-Object { Write-Host "  $_" }
    return @{ code = $code; text = ($output -join "`n") }
}

if( Test-Path $WorkDir ) { Remove-Item -Recurse -Force $WorkDir }
New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
Set-Location $WorkDir

# -- the launcher
$result = Invoke-Nextflow @('-version')
Check ($result.code -eq 0) 'launcher runs with Windows PowerShell'
Check ($result.text -match 'N E X T F L O W') 'launcher prints the version'

# -- without Docker the tasks cannot run on Windows
$noDocker = Join-Path $WorkDir 'no-docker.config'
Set-Content -Path $noDocker -Value 'docker.enabled = false'
$result = Invoke-Nextflow @('run', $pipeline, '-ansi-log', 'false', '-c', $noDocker)
Check ($result.code -ne 0) 'pipeline fails when Docker is not enabled'
Check ($result.text -match 'On Windows, tasks can only be executed with Docker') 'the error explains that Docker is required'
Remove-Item -Recurse -Force (Join-Path $WorkDir 'work'), (Join-Path $WorkDir '.nextflow') -ErrorAction SilentlyContinue

# -- first run
$result = Invoke-Nextflow @('run', $pipeline, '-ansi-log', 'false')
Check ($result.code -eq 0) "pipeline completes successfully (exit code $($result.code))"
Check ($result.text -match 'kernel: Linux; workdir: /[a-z]/') 'tasks run in a Linux container with mapped paths'
Check ($result.text -match 'lines: 2') 'env output is captured'
Check ($result.text -match 'attempt 2 ok') 'failed task is retried'

$upper = Join-Path $WorkDir 'results\upper.txt'
Check ((Test-Path $upper) -and ((Get-Content -Raw $upper) -ceq "HELLO WORLD`nFROM WINDOWS`n")) 'output is published with mode copy'
$count = Join-Path $WorkDir 'results\count.txt'
Check ((Test-Path $count) -and ((Get-Content -Raw $count).Trim() -eq '3')) 'output is published with the default mode'

$log = Get-Content -Raw (Join-Path $WorkDir '.nextflow.log')
Check ($log -match 'Process exceeded running time limit') 'task exceeding the time limit is stopped'

$wrappers = @(Get-ChildItem -Path (Join-Path $WorkDir 'work') -Recurse -Force -Filter '.command.run')
Check ($wrappers.Count -gt 0) 'task wrapper scripts are created'
$invalid = @($wrappers | Where-Object { $text = [IO.File]::ReadAllText($_.FullName); $text.Contains("`r") -or ($text -match '[A-Za-z]:\\') })
Check ($invalid.Count -eq 0) 'task wrapper scripts have Unix line endings and no Windows paths'

Start-Sleep -Seconds 5
$ErrorActionPreference = 'Continue'
$containers = (& docker ps -aq 2>&1 | Out-String).Trim()
$ErrorActionPreference = 'Stop'
Check ($containers -eq '') "no task container is left running [$containers]"

# -- resumed run
$result = Invoke-Nextflow @('run', $pipeline, '-ansi-log', 'false', '-resume')
Check ($result.code -eq 0) "resumed pipeline completes successfully (exit code $($result.code))"
Check ($result.text -match 'Cached process > UPPER') 'resumed run uses the task cache'
Check ($result.text -match 'Cached process > COUNT') 'resumed run uses the task cache for scratch tasks'

if( $failures -gt 0 ) {
    Write-Host "`n$failures check(s) failed"
    exit 1
}
Write-Host "`nAll checks passed"
