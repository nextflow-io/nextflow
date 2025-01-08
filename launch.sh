#!/bin/bash
#
#  Copyright 2013-2024, Seqera Labs
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

#set -e
#set -u
#set -o errexit

function resolve_link() {
    [[ ! -f $1 ]] && exit 1
    if command -v realpath &>/dev/null; then
      realpath "$1"
    elif command -v readlink &>/dev/null; then
      local target="$1"
      cd $(dirname $target); target=$(basename $target)
      while [ -L "$target" ]; do
        target="$(readlink "$target")"
        cd $(dirname $target); target=$(basename $target)
      done
      echo "$(cd "$(dirname "$target")"; pwd -P)/$target"
    else
      echo_yellow "WARN: Neither \`realpath\` nor \`readlink\` command can be found"
      exit 1
    fi
}

# the application 'base' folder
bin_dir=`dirname $(resolve_link $0)`
base_dir=$bin_dir

# define the java env
JAVA_BIN=java
if test -n "$JAVA_HOME"; then
	JAVA_BIN="$JAVA_HOME/bin/java"
fi

#
# debugging JVM options
#  export _JAVA_OPTIONS='-Xms1g -Xmx2g -verbose:gc -XX:+AggressiveOpts -XX:+UnlockDiagnosticVMOptions -XX:+UnlockExperimentalVMOptions -XX:+PrintFlagsFinal -XX:+TraceClassLoading -XX:+TraceClassUnloading -verbose:class'
#

#
# Variable definition
#
declare -a args=()
DEBUG=''
COLUMNS=${COLUMNS:-`tput cols 2> /dev/tty`}
MAIN_CLASS=${MAIN_CLASS:-'nextflow.cli.Launcher'}
JAVA_VER="$($JAVA_BIN -version 2>&1)"
if [ $? -ne 0 ]; then
    echo "${JAVA_VER:-Failed to launch the Java virtual machine}"
    exit 1
fi
JAVA_VER=$(echo "$JAVA_VER" | awk '/version/ {gsub(/"/, "", $3); print $3}')
major=${BASH_REMATCH[1]}
minor=${BASH_REMATCH[2]}
version_check="^(17|18|19|20|21|23)"
if [[ ! $JAVA_VER =~ $version_check ]]; then
    echo "Error: cannot find Java or it's a wrong version -- please make sure that Java 17 or higher is installed"
    exit 1
fi
JVM_ARGS+=" -Dfile.encoding=UTF-8 -XX:+TieredCompilation -XX:TieredStopAtLevel=1"
JVM_ARGS+=" --add-opens=java.base/java.lang=ALL-UNNAMED --add-opens=java.base/java.io=ALL-UNNAMED --add-opens=java.base/java.nio=ALL-UNNAMED --add-opens=java.base/java.net=ALL-UNNAMED --add-opens=java.base/java.nio.file.spi=ALL-UNNAMED --add-opens=java.base/java.util=ALL-UNNAMED --add-opens=java.base/java.util.concurrent.locks=ALL-UNNAMED --add-opens=java.base/java.util.concurrent.atomic=ALL-UNNAMED --add-opens=java.base/sun.nio.ch=ALL-UNNAMED --add-opens=java.base/sun.nio.fs=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.http=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.https=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.ftp=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.file=ALL-UNNAMED --add-opens=java.base/jdk.internal.misc=ALL-UNNAMED --add-opens=java.base/jdk.internal.vm=ALL-UNNAMED --add-opens=java.base/java.util.regex=ALL-UNNAMED"
[[ $NXF_ENABLE_VIRTUAL_THREADS == 'true' ]] && [[ "$JAVA_VER" =~ ^(19|20) ]] && JVM_ARGS+=" --enable-preview"
[[ "$JAVA_VER" =~ ^(21) ]] && [[ ! "$NXF_ENABLE_VIRTUAL_THREADS" ]] && NXF_ENABLE_VIRTUAL_THREADS=true

## flight recorded -- http://docs.oracle.com/javacomponents/jmc-5-4/jfr-runtime-guide/run.htm
##JVM_ARGS+=" -XX:+UnlockCommercialFeatures -XX:+FlightRecorder -XX:StartFlightRecording=duration=60s,filename=myrecording.jfr"
NXF_PLUGINS_DIR=${NXF_PLUGINS_DIR:-$base_dir/plugins}
NXF_PLUGINS_MODE=${NXF_PLUGINS_MODE:-dev}
NXF_PLUGINS_DEFAULT=${NXF_PLUGINS_DEFAULT:-true}
NXF_HOME=${NXF_HOME:-$HOME/.nextflow}
NXF_OPTS=${NXF_OPTS:-}
NXF_CLI="$0 $@"
NXF_REMOTE_DEBUG_PORT=${NXF_REMOTE_DEBUG_PORT:-5005}
export NXF_CLI
export COLUMNS
export NXF_PLUGINS_DIR
export NXF_PLUGINS_MODE
export NXF_PLUGINS_DEFAULT
export NXF_ENABLE_VIRTUAL_THREADS

# Yourkit profiling library
YOURKIT_AGENT=${YOURKIT_AGENT:-/Applications/YourKit-Java-Profiler-2021.11.app/Contents/Resources/bin/mac/libyjpagent.dylib}

#
# classpath when the application is compiled with gradle
#
[[ ! -f "$base_dir/.launch.classpath" ]] && echo "Missing '.launch.classpath' file -- create it by running: ./gradlew exportClasspath" && exit 1
CLASSPATH+=":`cat $base_dir/.launch.classpath`"

#
# Handle special program cli options
#
while [ "$*" != "" ]; do
  if [[ "$1" == '-debug' || "$1" == '-trace' ]]; then
    args+=("$1")

  elif [ "$1" == --yourkit ]; then
    JVM_ARGS+=" -agentpath:$YOURKIT_AGENT "
  elif [ "$1" == '--with-jrebel' ]; then
    if [ "$JREBEL_HOME" ]; then
    JVM_ARGS+=" -javaagent:$JREBEL_HOME/jrebel.jar -Drebel.log.file=./jrebel-client.log"
    else
    echo "WARN: To use JRebel define the JREBEL_HOME variable in environment"
    fi

  elif [ "$1" == '-remote-debug' ]; then
    DEBUG="-agentlib:jdwp=transport=dt_socket,server=y,suspend=y,address=$NXF_REMOTE_DEBUG_PORT"
    args+=("$1")
  elif [ "$1" == '-enable-checkpoint' ]; then
    mkdir -p crac-files
    JVM_ARGS+=" -XX:CRaCCheckpointTo=$PWD/crac-files"
  elif [ "$1" == '-checkpoint' ]; then
    jcmd $CLASSPATH JDK.checkpoint
    exit 0
  else
   args+=("$1")
  fi
  # move to the next option
  shift
done

# Show some variable when in DEBUG mode
if [ "$DEBUG" != '' ]; then
  echo Launch environment
  echo ------------------
  echo base_dir: $base_dir
  echo jvmargs: $JVM_ARGS
  echo debug: $DEBUG
  echo classpath:
  echo $CLASSPATH | tr ":" "\n" | sort
  echo ''
  echo Launching it!
  echo ------------------
fi

# Launch the APP
exec $JAVA_BIN $JVM_ARGS $DEBUG $NXF_OPTS -cp "$CLASSPATH" "$MAIN_CLASS" "${args[@]}"
