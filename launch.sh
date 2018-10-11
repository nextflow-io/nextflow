#!/bin/bash
#
# Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
# Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
#
#   This file is part of 'Nextflow'.
#
#   Nextflow is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Nextflow is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
#

#set -e
#set -u
#set -o errexit

# the application 'base' folder
bin_dir=`dirname "$0"`
bin_dir=`cd "$bin_dir"; pwd`
base_dir=$bin_dir

# define the java env
java=java
if test -n "$JAVA_HOME"; then
	java="$JAVA_HOME/bin/java"
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
JAVA_VER="$(java -version 2>&1)"
if [ $? -ne 0 ]; then
    echo "${JAVA_VER:-Failed to launch the Java virtual machine}"
    exit 1
fi
JAVA_VER=$(echo "$JAVA_VER" | awk '/version/ {gsub(/"/, "", $3); print $3}')
major=${BASH_REMATCH[1]}
minor=${BASH_REMATCH[2]}
version_check="^(1.8|9|10)"
if [[ ! $JAVA_VER =~ $version_check ]]; then
    echo "Error: cannot find Java or it's a wrong version -- please make sure that Java 8 or higher is installed"
    exit 1
fi
JVM_ARGS+=" -Dfile.encoding=UTF-8 -noverify -XX:+TieredCompilation -XX:TieredStopAtLevel=1"
[[ $JAVA_VER =~ ^9|10 ]] && JVM_ARGS+=" --add-opens=java.base/java.lang=ALL-UNNAMED --add-opens=java.base/java.io=ALL-UNNAMED --add-opens=java.base/java.nio=ALL-UNNAMED --add-opens=java.base/java.util=ALL-UNNAMED --add-opens=java.base/sun.nio.ch=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.http=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.https=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.ftp=ALL-UNNAMED --add-opens=java.base/sun.net.www.protocol.file=ALL-UNNAMED --add-opens=java.base/jdk.internal.misc=ALL-UNNAMED --add-modules=java.xml.bind --illegal-access=deny"

## flight recorded -- http://docs.oracle.com/javacomponents/jmc-5-4/jfr-runtime-guide/run.htm
##JVM_ARGS+=" -XX:+UnlockCommercialFeatures -XX:+FlightRecorder -XX:StartFlightRecording=duration=60s,filename=myrecording.jfr"
NXF_HOME=${NXF_HOME:-$HOME/.nextflow}
NXF_OPTS=${NXF_OPTS:-}
EXTRAE_CONFIG_FILE=${EXTRAE_CONFIG_FILE:-$NXF_HOME/extrae/config}
NXF_CLI="$0 $@"
export NXF_CLI
export COLUMNS
export EXTRAE_CONFIG_FILE

DRIP_INIT_CLASS=nextflow.cli.DripMain
DRIP_INIT=''

EXTRAE_CONFIG_FILE=${EXTRAE_CONFIG_FILE:-$NXF_HOME/extrae/config}

#
# classpath when the application is compiled with gradle
#
if [ -e "$base_dir/build/libs" ]; then
  CLASSPATH=`ls $base_dir/build/libs/* | egrep 'nextflow-[0-9]+\.[0-9]+\.[0-9a-z]+(-[A-Z0-9]+)?.jar$'`

  # -- append runtime libraries
  [[ ! -f "$base_dir/.launch.classpath" ]] && echo "Missing '.launch.classpath' file -- create it by running: ./gradlew exportClasspath" && exit 1
  CLASSPATH+=":`cat $base_dir/.launch.classpath`"
  
else
  echo "Missing application libraries -- Compile Nextflow by using the 'compile.sh' script"
  exit 1
fi


#
# Handle special program cli options
#
while [ "$*" != "" ]; do
  if [[ "$1" == '-debug' || "$1" == '-trace' ]]; then
    args+=("$1")

  elif [ "$1" == --with-yourkit ]; then 
    JVM_ARGS+=" -agentpath:/Applications/YourKit_Java_Profiler_2014_build_14104.app/bin/mac/libyjpagent.jnilib=onexit=snapshot,tracing,dir=$PWD/snapshot "
  elif [ "$1" == '--with-jrebel' ]; then
    if [ "$JREBEL_HOME" ]; then
    JVM_ARGS+=" -javaagent:$JREBEL_HOME/jrebel.jar -Drebel.log.file=./jrebel-client.log"
    else
    echo "WARN: To use JRebel define the JREBEL_HOME variable in environment"
    fi

  elif [ "$1" == '-remote-debug' ]; then
    DEBUG='-agentlib:jdwp=transport=dt_socket,server=y,suspend=y,address=8010'

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
exec java $JVM_ARGS $DEBUG $NXF_OPTS -cp "$CLASSPATH" "$MAIN_CLASS" "${args[@]}"
