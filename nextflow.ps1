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

# Nextflow launcher for Windows (Windows PowerShell 5.1 or PowerShell 7+).
#
# On Windows, tasks are executed with Docker: enable it with `docker.enabled = true`
# and define a container for each process. See https://docs.seqera.io/nextflow/install

$ErrorActionPreference = 'Stop'

function Get-EnvOrDefault([string]$name, [string]$default) {
    $value = [Environment]::GetEnvironmentVariable($name)
    if( $value ) { return $value } else { return $default }
}

function Fail([string]$message) {
    [Console]::Error.WriteLine($message)
    exit 1
}

$NXF_VER = Get-EnvOrDefault 'NXF_VER' '26.10.3'
$NXF_ORG = Get-EnvOrDefault 'NXF_ORG' 'nextflow-io'
$NXF_HOME = Get-EnvOrDefault 'NXF_HOME' (Join-Path ([Environment]::GetFolderPath('UserProfile')) '.nextflow')
$NXF_BASE = Get-EnvOrDefault 'NXF_BASE' 'https://www.nextflow.io/releases'
$NXF_DIST = Get-EnvOrDefault 'NXF_DIST' (Join-Path $NXF_HOME 'framework')
$NXF_TEMP = Get-EnvOrDefault 'NXF_TEMP' ([IO.Path]::GetTempPath())
$NXF_JAR = Get-EnvOrDefault 'NXF_JAR' "nextflow-$NXF_VER-one.jar"
$NXF_BIN = Get-EnvOrDefault 'NXF_BIN' (Join-Path (Join-Path $NXF_DIST $NXF_VER) $NXF_JAR)
$NXF_URL = Get-EnvOrDefault 'NXF_URL' "$NXF_BASE/v$NXF_VER/$NXF_JAR"

if( $NXF_VER -notmatch '^(\d+)\.(\d+)\.(\d.*)$' ) {
    Fail "Not a valid Nextflow version: $NXF_VER"
}

# -- parse the launcher options, everything is passed to Nextflow as it is
$jvmOpts = @()
$command = $null
foreach( $arg in $args ) {
    if( !$command -and $arg -like '-D*' ) {
        $jvmOpts += $arg
    }
    elseif( !$command -and ($arg -eq '-bg' -or $arg -eq '-remote-debug') ) {
        Fail "The option '$arg' is not supported by the Nextflow launcher for Windows"
    }
    elseif( !$command -and ($arg -eq '-self-update' -or $arg -eq 'self-update') ) {
        Fail "Self update is not supported by the Nextflow launcher for Windows -- download the latest 'nextflow.ps1' or set the NXF_VER variable to the version to use"
    }
    elseif( !$command -and $arg -eq '-download' ) {
        $versionDir = Split-Path $NXF_BIN -Parent
        if( Test-Path $versionDir ) { Remove-Item -Recurse -Force $versionDir }
        $command = '-download'
    }
    elseif( !$command -and $arg -notlike '-*' ) {
        $command = $arg
    }
}

# -- find the Java VM
$javaHome = Get-EnvOrDefault 'NXF_JAVA_HOME' $env:JAVA_HOME
if( $javaHome -and (Test-Path (Join-Path $javaHome 'bin\java.exe')) ) {
    $JAVA_CMD = Join-Path $javaHome 'bin\java.exe'
}
elseif( $javaHome -and (Test-Path (Join-Path $javaHome 'bin/java')) ) {
    $JAVA_CMD = Join-Path $javaHome 'bin/java'
}
else {
    $found = Get-Command java -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
    if( !$found ) {
        Fail "Cannot find Java -- Nextflow requires Java 17 or later (up to 26); install it and make sure the 'java' command is available or JAVA_HOME is defined"
    }
    $JAVA_CMD = $found.Source
}

# -- check the Java version, using the JDK `release` file when available
$JAVA_VER = $null
$releaseFile = Join-Path (Split-Path (Split-Path $JAVA_CMD -Parent) -Parent) 'release'
if( Test-Path $releaseFile ) {
    $line = Select-String -Path $releaseFile -Pattern '^JAVA_VERSION="?([^"]+)' | Select-Object -First 1
    if( $line ) { $JAVA_VER = $line.Matches[0].Groups[1].Value }
}
if( !$JAVA_VER ) {
    # note: native command stderr must not be turned into a terminating error
    $ErrorActionPreference = 'Continue'
    $output = (& $JAVA_CMD -version 2>&1 | ForEach-Object { "$_" }) -join "`n"
    $ErrorActionPreference = 'Stop'
    if( $output -match 'version "?([0-9][^" ]*)' ) { $JAVA_VER = $Matches[1] }
}
if( $JAVA_VER -notmatch '^(17|18|19|20|21|22|23|24|25|26)' ) {
    Fail "Cannot find Java or it's a wrong version -- please make sure that Java 17 or later (up to 26) is installed`nJava command: $JAVA_CMD; detected version: $JAVA_VER"
}

# -- download the Nextflow runtime when missing
if( !(Test-Path $NXF_BIN) ) {
    $target = Split-Path $NXF_BIN -Parent
    New-Item -ItemType Directory -Force -Path $target | Out-Null
    $tmpFile = Join-Path $target "$NXF_JAR.download"
    Write-Host 'Downloading nextflow dependencies. It may require a few seconds, please wait ..'
    try {
        [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
        $progress = $ProgressPreference
        $ProgressPreference = 'SilentlyContinue'
        Invoke-WebRequest -Uri $NXF_URL -OutFile $tmpFile -UseBasicParsing
        $ProgressPreference = $progress
        Move-Item -Force $tmpFile $NXF_BIN
    }
    catch {
        if( Test-Path $tmpFile ) { Remove-Item -Force $tmpFile }
        Fail "ERROR: Cannot download nextflow required file -- make sure you can connect to the internet`n`nAlternatively you can try to download this file:`n    $NXF_URL`n`nand save it as:`n    $NXF_BIN"
    }
}
if( $command -eq '-download' ) {
    exit 0
}

# -- clean up and export the environment
# see https://github.com/nextflow-io/nextflow/issues/1716
Remove-Item Env:JAVA_TOOL_OPTIONS -ErrorAction SilentlyContinue
$env:NXF_HOME = $NXF_HOME
$env:NXF_ORG = $NXF_ORG
$env:NXF_CLI = "nextflow $($args -join ' ')"
$env:NXF_PACK = 'one'
$env:JAVA_CMD = $JAVA_CMD
if( !$env:COLUMNS ) {
    try { $env:COLUMNS = [string]$Host.UI.RawUI.WindowSize.Width } catch { }
}

# -- the Java VM options
$javaOpts = @('-Dfile.encoding=UTF-8', '-Dcom.sun.security.enableAIAcaIssuers=true', '-Djava.awt.headless=true')
$javaOpts += @('-XX:TieredStopAtLevel=1', '-XX:ReservedCodeCacheSize=240m')
if( "$env:NXF_OPTS $env:NXF_JVM_ARGS" -notmatch '-XX:\+Use[A-Za-z0-9]+GC' ) {
    $javaOpts += '-XX:+UseSerialGC'
}
# the task cache database deletes files that are memory mapped, which is not allowed on Windows
$javaOpts += '-Dleveldb.mmap=false'
if( $env:NXF_OPTS ) { $javaOpts += ($env:NXF_OPTS -split '\s+' | Where-Object { $_ }) }
if( $env:NXF_JVM_ARGS ) { $javaOpts += ($env:NXF_JVM_ARGS -split '\s+' | Where-Object { $_ }) }
if( $env:NXF_TEMP ) { $javaOpts += "-Djava.io.tmpdir=$NXF_TEMP" }
$javaOpts += $jvmOpts

$javaOpts += @(
    '--add-opens=java.base/java.lang=ALL-UNNAMED',
    '--add-opens=java.base/java.io=ALL-UNNAMED',
    '--add-opens=java.base/java.nio=ALL-UNNAMED',
    '--add-opens=java.base/java.net=ALL-UNNAMED',
    '--add-opens=java.base/java.util=ALL-UNNAMED',
    '--add-opens=java.base/java.util.concurrent.locks=ALL-UNNAMED',
    '--add-opens=java.base/java.util.concurrent.atomic=ALL-UNNAMED',
    '--add-opens=java.base/jdk.internal.misc=ALL-UNNAMED',
    '--add-opens=java.base/java.util.regex=ALL-UNNAMED'
)
if( $JAVA_VER -match '^(24|25|26)' ) {
    $javaOpts += @('--enable-native-access=ALL-UNNAMED', '--sun-misc-unsafe-memory-access=allow')
}
if( $env:NXF_ENABLE_VIRTUAL_THREADS -eq 'true' -and $JAVA_VER -match '^(19|20)' ) {
    $javaOpts += '--enable-preview'
}

# -- finally run it
$ErrorActionPreference = 'Continue'
& $JAVA_CMD @javaOpts -jar $NXF_BIN @args
exit $LASTEXITCODE
