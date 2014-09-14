#!/bin/bash

#
#  Copyright (c) 2012, the authors.
#
#    This file is part of 'Nextflow'.
#
#    Nextflow is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Nextflow is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
# Variable definition
#
declare -a args=()
DEBUG=''
COLUMNS=${COLUMNS:-`tput cols 2> /dev/tty`}
MAIN_CLASS='nextflow.cli.Launcher'
JVM_ARGS+=" -Djava.awt.headless=true -noverify"
NXF_HOME=${NXF_HOME:-$HOME/.nextflow}
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
  CLASSPATH=`ls $base_dir/build/libs/nextflow-*.jar`

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
  if [[ "$1" == '--debug' || "$1" == '--trace' ]]; then
    args+=("$1")

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
exec java $JVM_ARGS $DEBUG -noverify -cp "$CLASSPATH" "$MAIN_CLASS" "${args[@]}"
