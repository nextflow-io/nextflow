#!/usr/bin/env bash
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

[[ "$NXF_DEBUG" == 'x' ]] && set -x
NXF_VER=${NXF_VER:-'24.10.5'}
NXF_ORG=${NXF_ORG:-'nextflow-io'}
NXF_HOME=${NXF_HOME:-$HOME/.nextflow}
NXF_PROT=${NXF_PROT:-'https'}
NXF_BASE=${NXF_BASE:-$NXF_PROT://www.nextflow.io/releases}
NXF_TEMP=${NXF_TEMP:-$TMPDIR}
NXF_DIST=${NXF_DIST:-$NXF_HOME/framework}
NXF_CLI="$0 $@"
NXF_CLI_OPTS=${NXF_CLI_OPTS:-}
NXF_REMOTE_DEBUG_PORT=${NXF_REMOTE_DEBUG_PORT:-5005}

function cmp_ver() {
    local IFS=.
    local i ver1 ver2
    read -r -a ver1 <<< "${1//./ }"
    read -r -a ver2 <<< "${2//./ }"

    # Compare major, minor, and patch numbers
    for ((i=0; i<3; i++)); do
        ver1[i]=${ver1[i]//[^0-9]}
        ver2[i]=${ver2[i]//[^0-9]}
        [[ ${ver1[i]:-0} -lt ${ver2[i]:-0} ]] && echo -1 && return
        [[ ${ver1[i]:-0} -gt ${ver2[i]:-0} ]] && echo 1 && return
    done

    # Extract suffixes for comparison
    local suffix1="${1##*.}"
    local suffix2="${2##*.}"
    suffix1="${suffix1//[0-9]}"
    suffix2="${suffix2//[0-9]}"

    # Compare suffixes
    [[ -z $suffix1 && -n $suffix2 ]] && echo 1 && return
    [[ -n $suffix1 && -z $suffix2 ]] && echo -1 && return
    [[ $suffix1 < $suffix2 ]] && echo -1 && return
    [[ $suffix1 > $suffix2 ]] && echo 1 && return

    # Versions are equal
    echo 0
}

# if the nextflow version is greater or equals to "24.07.0-edge" the new shadow jar launcher
# should be used, otherwise fallback on the legacy behavior setting the variable NXF_LEGACY_LAUNCHER
NXF_LEGACY_LAUNCHER=1
if [[ $(cmp_ver "$NXF_VER" "24.07.0-edge") -ge 0 ]]; then
  unset NXF_LEGACY_LAUNCHER
fi

export NXF_CLI
export NXF_ORG
export NXF_HOME

if [[ $TERM && $TERM != 'dumb' ]]; then
if command -v tput &>/dev/null; then
GREEN=$(tput setaf 2; tput bold)
YELLOW=$(tput setaf 3)
RED=$(tput setaf 1)
NORMAL=$(tput sgr0)
fi
fi

function echo_red() {
    >&2 echo -e "$RED$*$NORMAL"
}

function echo_green() {
    echo -e "$GREEN$*$NORMAL"
}

function echo_yellow() {
    >&2 echo -e "$YELLOW$*$NORMAL"
}

function die() {
  echo_red "$*"
  exit 1
}

function get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

function get() {
    if command -v curl &>/dev/null; then
        GET="curl -fsSL '$1' -o '$2'"
    elif command -v wget &>/dev/null; then
        GET="wget '$1' -O '$2' >/dev/null 2>&1"
    else
        echo_red "ERROR: Cannot find 'curl' nor 'wget' utility --  please install one of them"
        exit 1
    fi

    printf "Downloading nextflow dependencies. It may require a few seconds, please wait .. "
    eval $GET; status=$?
    printf "\r\033[K"
    if [ $status -ne 0 ]; then
        echo_red "ERROR: Cannot download nextflow required file -- make sure you can connect to the internet"
        echo ""
        echo "Alternatively you can try to download this file:"
        echo "    $1"
        echo ""
        echo "and save it as:"
        echo "    ${3:-$2}"
        echo ""
        exit 1
    fi
}

function get_ver() {
    if command -v curl &>/dev/null; then
        curl -fsSL "$1"
    elif command -v wget &>/dev/null; then
        wget "$1" >/dev/null 2>&1
    else
        echo_red "ERROR: Cannot find 'curl' nor 'wget' utility -- please install one of them"
        exit 1
    fi
}

function make_temp() {
    local base=${NXF_TEMP:=$PWD}
    if [ "$(uname)" = 'Darwin' ]; then mktemp "${base}/nxf-tmp.XXXXXX" || exit $?
    else mktemp -t nxf-tmp.XXXXXX -p "${base}" || exit $?
    fi
}

function resolve_link() {
    [[ ! -f $1 ]] && exit 1
    if command -v realpath &>/dev/null; then
      realpath "$1"
    elif command -v readlink &>/dev/null; then
      local target="$1"
      cd "$(dirname "$target")"; target="$(basename "$target")"
      while [ -L "$target" ]; do
        target="$(readlink "$target")"
        cd "$(dirname "$target")"; target="$(basename "$target")"
      done
      echo "$(cd "$(dirname "$target")"; pwd -P)/$target"
    else
      echo_yellow "WARN: Neither \`realpath\` nor \`readlink\` command can be found"
      exit 1
    fi
}

function current_ver() {
  [[ $NXF_EDGE == 1 || $NXF_VER == *"-edge" ]] && printf 'edge' || printf 'latest'
}

function install() {
    local tmpfile=$(make_temp)
    local version=$(set +u; [[ $NXF_VER ]] && printf "v$NXF_VER" || current_ver)
    local action="a=${2:-default}"
    get "$NXF_BASE/$version/nextflow?$action" "$tmpfile" "$1" || exit $?
    mv "$tmpfile" "$1" || exit $?
    chmod +x "$1" || exit $?
    bash "$1" -download || exit $?
    echo ''
    echo -e $'Nextflow installation completed. Please note:'
    echo -e $'- the executable file `nextflow` has been created in the folder:' $(dirname $1)
    if [[ ! "$PATH" =~ (^|:)"$(dirname $1)"(:|$) ]]; then
    echo -e $'- you may complete the installation by moving it to a directory in your $PATH'
    fi
    echo ''
}

function check_latest() {
    [[ $cmd != run ]] && return 0
    [[ $NXF_OFFLINE == true || $NXF_DISABLE_CHECK_LATEST == true ]] && return 0
    local latest=$(get_ver "$NXF_BASE/$(current_ver)/version?current=$NXF_VER")
    if [[ -n "$latest" && $(cmp_ver "$latest" "$NXF_VER") -gt 0 ]]; then
      echo_yellow "Nextflow $latest is available - Please consider updating your version to it"
    fi
}

function launch_nextflow() {
    # the launch command line
    local cmdline=()
    # remove leading and trailing double-quotes
    for x in "${launcher[@]}"; do
        x="${x%\"}"
        x="${x#\"}"
        cmdline+=("$x") 
    done 

    if [[ "$bg" ]]; then
      local pid_file="${NXF_PID_FILE:-.nextflow.pid}"
      cmdline+=("${args[@]}")
      exec "${cmdline[@]}" &
      disown
      echo $! > "$pid_file"
      exit 0
    fi

    cmdline+=("${args[@]}")
    exec "${cmdline[@]}"
    exit 1
}

# check self-install
if [ "$0" = "bash" ] || [[ "$0" =~ .*/bash ]]; then
    if [ -d nextflow ]; then
        echo 'Please note:'
        echo "- The install procedure needs to create a file named 'nextflow' in this folder, but a directory with this name already exists."
        echo "- Please renamed/delete that directory, or execute the Nextflow install procedure in another folder."
        echo ''
        exit 1
    fi
    install "$PWD/nextflow" install
    exit 0
fi

# clean up env
# see https://github.com/nextflow-io/nextflow/issues/1716
unset JAVA_TOOL_OPTIONS

# parse the command line
bg=''
declare -a jvmopts=()
declare -a args=("$@")
declare -a commands=(clone config drop help history info ls pull run view node console kuberun)
# $NXF_CLI_OPTS allow to pass arbitrary cli opts via the environment
# note: do not wrap with quotes because the variable can be used to specify more than on option separate by blanks
[ "$NXF_CLI_OPTS" ] && args+=($NXF_CLI_OPTS)

cmd=''
while [[ $# != 0 ]]; do
    case $1 in
    -D*)
      if [[ ! "$cmd" ]]; then
      jvmopts+=("$1")
      fi
      ;;
    -bg)
      bg=1
      ;;
    -remote-debug)
      echo_yellow "Enabling script debugging - continue the execution launching the remote VM debugger in your favourite IDE using port $NXF_REMOTE_DEBUG_PORT"
      remote_debug=1
      ;;
    -download)
      if [[ ! "$cmd" ]]; then
      rm -rf "$NXF_DIST/$NXF_VER" || exit $?
      bash "$0" -version || exit $?
      exit 0
      fi
      ;;
    -self-update|self-update)
      if [[ ! "$cmd" ]]; then
      [[ -z $NXF_EDGE && $NXF_VER = *-edge ]] && NXF_EDGE=1
      unset NXF_VER
      install "$0" update
      exit 0
      fi
      ;;
    *)
      [[ $1 && $1 != -* && ! "$cmd" && ${commands[*]} =~ $1 ]] && cmd=$1
      ;;
    esac
    shift
done

CAPSULE_LOG=${CAPSULE_LOG:=''}
CAPSULE_RESET=${CAPSULE_RESET:=''}
CAPSULE_CACHE_DIR=${CAPSULE_CACHE_DIR:="$NXF_HOME/capsule"}

NXF_PACK=one
NXF_MODE=${NXF_MODE:-''}
NXF_JAR=${NXF_JAR:-nextflow-$NXF_VER-$NXF_PACK.jar}
NXF_BIN=${NXF_BIN:-$NXF_DIST/$NXF_VER/$NXF_JAR}
NXF_PATH=$(dirname "$NXF_BIN")
NXF_URL=${NXF_URL:-$NXF_BASE/v$NXF_VER/$NXF_JAR}
NXF_HOST=${HOSTNAME:-localhost}
[[ $NXF_LAUNCHER ]] || NXF_LAUNCHER=${NXF_HOME}/tmp/launcher/nextflow-${NXF_PACK}_${NXF_VER}/${NXF_HOST}
# both NXF_GRAB and NXF_CLASSPATH are not supported any more as of version 24.04.7-edge
NXF_GRAB=${NXF_GRAB:-''}
NXF_CLASSPATH=${NXF_CLASSPATH:-''}

# Determine the path to this file
if [[ $NXF_PACK = dist ]]; then
    NXF_BIN=$(which "$0" 2>/dev/null)
    [ $? -gt 0 -a -f "$0" ] && NXF_BIN="./$0"
fi

# use nextflow custom java home path
if [[ "$NXF_JAVA_HOME" ]]; then
  JAVA_HOME="$NXF_JAVA_HOME"
  unset JAVA_CMD
fi
# Determine the Java command to use to start the JVM.
if [ ! -x "$JAVA_CMD" ] ; then
    if [ -d "$JAVA_HOME" ] ; then
        if [ -x "$JAVA_HOME/jre/sh/java" ] ; then
            # IBM's JDK on AIX uses strange locations for the executables
            JAVA_CMD="$JAVA_HOME/jre/sh/java"
        else
            JAVA_CMD="$JAVA_HOME/bin/java"
        fi
    elif [ -x /usr/libexec/java_home ]; then
        JAVA_CMD="$(/usr/libexec/java_home -v 11+ 2>/dev/null)/bin/java" || JAVA_CMD=java
    else
        JAVA_CMD="$(which java)" || JAVA_CMD=java
    fi
fi

# Retrieve the java version from a NF local file
JAVA_KEY="$NXF_HOME/tmp/ver/$(resolve_link "$JAVA_CMD" | sed 's@/@.@g')"
if [ -f "$JAVA_KEY" ]; then
  JAVA_VER="$(cat "$JAVA_KEY")"
else
  JAVA_VER="$("$JAVA_CMD" $NXF_OPTS -version 2>&1)"
  if [ $? -ne 0 ]; then
      getstarted_web="https://www.nextflow.io/docs/latest/getstarted.html"
      echo_red "${JAVA_VER:-Failed to launch the Java virtual machine}"
      echo_red "NOTE: Nextflow needs a Java virtual machine to run. To this end:
 - make sure a \`java\` command can be found; or
 - manually define the variables JAVA_HOME to point to an existing installation; or
 - install a Java virtual machine, for instance through https://sdkman.io (read the docs);
 - for more details please refer to the Nextflow Get Started page at http://docs.nextflow.io."
      echo_yellow "NOTE: Nextflow is trying to use the Java VM defined by the following environment variables:\n JAVA_CMD: $JAVA_CMD\n NXF_OPTS: $NXF_OPTS\n"
      exit 1
  fi
  JAVA_VER=$(echo "$JAVA_VER" | awk '/version/ {gsub(/"/, "", $3); print $3}')
  # check NF version
  if [[ ! $NXF_VER =~ ([0-9]+)\.([0-9]+)\.([0-9].*) ]]; then
    echo_red "Not a valid Nextflow version: $NXF_VER"
    exit 1
  fi
  major=${BASH_REMATCH[1]}
  minor=${BASH_REMATCH[2]}
  # legacy version - Java 7/8 only
  if [ $major -eq 0 ] && [ $minor -lt 26 ]; then
    version_check="^(1.7|1.8)"
    version_message="Java 7 or 8"
  else
    version_check="^(1.8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23)"
    version_message="Java 8 or later (up to 22)"
  fi
  if [[ ! $JAVA_VER =~ $version_check ]]; then
      echo_red "ERROR: Cannot find Java or it's a wrong version -- please make sure that $version_message is installed"
      if [[ "$NXF_JAVA_HOME" ]]; then
      echo_yellow "NOTE: Nextflow is trying to use the Java VM defined by the following environment variables:\n JAVA_CMD: $JAVA_CMD\n NXF_JAVA_HOME: $NXF_JAVA_HOME\n"
      else
      echo_yellow "NOTE: Nextflow is trying to use the Java VM defined by the following environment variables:\n JAVA_CMD: $JAVA_CMD\n JAVA_HOME: $JAVA_HOME\n"
      fi
      exit 1
  fi
  if [[ ! $JAVA_VER =~ ^(11|12|13|14|15|16|17|18|19|20|21|22|23) ]]; then
      echo_yellow "NOTE: Nextflow is not tested with Java $JAVA_VER -- It's recommended the use of version 11 up to 23\n"
  fi
  mkdir -p "$(dirname "$JAVA_KEY")"
  [[ -f $JAVA_VER ]] && echo $JAVA_VER > "$JAVA_KEY"
fi

# Verify nextflow jar is available
if [ ! -f "$NXF_BIN" ]; then
    [ -f "$NXF_PATH" ] && rm "$NXF_PATH"
    mkdir -p "$NXF_PATH" || exit $?
    tmpfile=$(make_temp)
    get "$NXF_URL" "$tmpfile" "$NXF_BIN"
    mv "$tmpfile" "$NXF_BIN"
fi

COLUMNS=${COLUMNS:-`tty -s && tput cols 2>/dev/null || true`}
declare -a JAVA_OPTS=()
JAVA_OPTS+=(-Dfile.encoding=UTF-8 -Dcapsule.trampoline -Dcapsule.java.cmd="$JAVA_CMD" -Dcom.sun.security.enableAIAcaIssuers=true)
if [[ $cmd == console ]]; then bg=1;
else JAVA_OPTS+=(-Djava.awt.headless=true)
fi

if [[ $NXF_LEGACY_LAUNCHER ]]; then
  [[ "$JAVA_HOME" ]] && JAVA_OPTS+=(-Dcapsule.java.home="$JAVA_HOME")
  [[ "$CAPSULE_LOG" ]] && JAVA_OPTS+=(-Dcapsule.log=$CAPSULE_LOG)
  [[ "$CAPSULE_RESET" ]] && JAVA_OPTS+=(-Dcapsule.reset=true)
fi
[[ "$cmd" != "run" && "$cmd" != "node" ]] && JAVA_OPTS+=(-XX:+TieredCompilation -XX:TieredStopAtLevel=1)
[[ "$NXF_OPTS" ]] && JAVA_OPTS+=($NXF_OPTS)
[[ "$NXF_CLASSPATH" ]] && export NXF_CLASSPATH
[[ "$NXF_GRAB" ]] && export NXF_GRAB
[[ "$COLUMNS" ]] && export COLUMNS
[[ "$NXF_TEMP" ]] && JAVA_OPTS+=(-Djava.io.tmpdir="$NXF_TEMP")
[[ "${jvmopts[@]}" ]] && JAVA_OPTS+=("${jvmopts[@]}")
export JAVA_CMD
export CAPSULE_CACHE_DIR
export NXF_PLUGINS_DIR
export NXF_PLUGINS_MODE
export NXF_PLUGINS_DEFAULT
export NXF_PACK
export NXF_ENABLE_VIRTUAL_THREADS

# lookup the a `md5` command
if hash md5sum 2>/dev/null; then MD5=md5sum;
elif hash gmd5sum 2>/dev/null; then MD5=gmd5sum;
elif hash md5 2>/dev/null; then MD5=md5;
else MD5=''
fi

# when no md5 command is available fallback on default execution
if [ ! "$MD5" ] || [ "$CAPSULE_RESET" ]; then
    launcher=($("$JAVA_CMD" "${JAVA_OPTS[@]}" -jar "$NXF_BIN"))
    launch_nextflow
    exit 1
fi

# creates a md5 unique for the given variables
env_md5() {
cat <<EOF | $MD5 | cut -f1 -d' '
$JAVA_CMD
$JAVA_VER
${JAVA_OPTS[@]}
$NXF_HOME
$NXF_VER
$NXF_OPTS
$NXF_GRAB
$NXF_CLASSPATH
$NXF_JVM_ARGS
$NXF_ENABLE_VIRTUAL_THREADS
EOF
}

# checked if a cached classpath file exists and it newer that the nextflow boot jar file
LAUNCH_FILE="${NXF_LAUNCHER}/classpath-$(env_md5)"

if [ -s "$LAUNCH_FILE" ] && [ "$LAUNCH_FILE" -nt "$NXF_BIN" ] && [[ "$remote_debug" -ne 1 ]]; then
    declare -a launcher="($(cat "$LAUNCH_FILE"))"
else
    if [[ $NXF_LEGACY_LAUNCHER ]]; then
      # otherwise run the capsule and get the result classpath in the 'launcher' and save it to a file
      cli=($("$JAVA_CMD" "${JAVA_OPTS[@]}" -jar "$NXF_BIN"))
      [[ $? -ne 0 ]] && echo_red 'Unable to initialize nextflow environment' && exit 1
    else
      # otherwise parse the command and get the result classpath in the 'launcher' and save it to a file
      cli=("\"$JAVA_CMD\"" "${JAVA_OPTS[@]}" -jar "$NXF_BIN")
    fi

    # first string between double quotes is the full path to java, also blank spaces are included
    # remainder string are arguments
    # we extract first part into `cmd_base`` and remainder into `cmd_tail`` and convert them to array as previous version
    cmd_pattern='"([^"]*)"(.*)'
    [[ "${cli[@]}" =~ $cmd_pattern ]]
    declare -a cmd_base="(${BASH_REMATCH[1]})"
    declare -a cmd_tail="(${BASH_REMATCH[2]})"

    launcher="${cmd_base[@]}"
    [[ "$NXF_JVM_ARGS" ]] && launcher+=($NXF_JVM_ARGS)
    [[ "$remote_debug" ]] && launcher+=(-agentlib:jdwp=transport=dt_socket,server=y,suspend=y,address=$NXF_REMOTE_DEBUG_PORT)

    if [[ "$JAVA_VER" =~ ^(9|10|11|12|13|14|15|16|17|18|19|20|21|22|23) ]]; then
      launcher+=(--add-opens=java.base/java.lang=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.io=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.nio=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.net=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.util=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.util.concurrent.locks=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.util.concurrent.atomic=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.nio.file.spi=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/sun.nio.ch=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/sun.nio.fs=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/sun.net.www.protocol.http=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/sun.net.www.protocol.https=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/sun.net.www.protocol.ftp=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/sun.net.www.protocol.file=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/jdk.internal.misc=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/jdk.internal.vm=ALL-UNNAMED)
      launcher+=(--add-opens=java.base/java.util.regex=ALL-UNNAMED)
      if [[ "$NXF_ENABLE_VIRTUAL_THREADS" == 'true' ]]; then
        if [[ "$JAVA_VER" =~ ^(19|20) ]]; then launcher+=(--enable-preview)
        elif [[ ! "$JAVA_VER" =~ ^(21|22|23) ]]; then die "Virtual threads require Java 19 or later - current version $JAVA_VER"
        fi
      fi
      launcher+=("${cmd_tail[@]}")
    else
      launcher+=("${cmd_tail[@]}")
    fi

    # create the launch file only when using the legacy launcher (capsule)
    if [[ $NXF_LEGACY_LAUNCHER ]]; then
    # Don't show errors if the LAUNCH_FILE can't be created
    if mkdir -p "${NXF_LAUNCHER}" 2>/dev/null; then
        STR=''
        for x in "${launcher[@]}"; do
        [[ "$x" != "\"-Duser.dir=$PWD\"" ]] && [[ ! "$x" == *"-agentlib:jdwp"* ]] && STR+=$(printf '%q ' "$x")
        done
        printf "$STR">"$LAUNCH_FILE"
    else
        echo_yellow "Warning: Couldn't create cached classpath folder: $NXF_LAUNCHER -- Maybe NXF_HOME is not writable?"
    fi
    fi

fi

# check for latest version
check_latest
# finally run it
launch_nextflow
