#!/bin/bash
#BSUB -x 1
#BSUB -y 2
### ---
### name: 'Hello 1'
### ...
set -e
set -u
NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
NXF_ENTRY=${1:-nxf_main}


nxf_sleep() {
  sleep $1 2>/dev/null || sleep 1;
}

nxf_date() {
    local ts=$(date +%s%3N);
    if [[ ${#ts} == 10 ]]; then echo ${ts}000
    elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
    elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
    elif [[ ${#ts} == 13 ]]; then echo $ts
    else echo "Unexpected timestamp value: $ts"; exit 1
    fi
}

nxf_env() {
    echo '============= task environment ============='
    env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
    echo '============= task output =================='
}

nxf_kill() {
    declare -a children
    while read P PP;do
        children[$PP]+=" $P"
    done < <(ps -e -o pid= -o ppid=)

    kill_all() {
        [[ $1 != $$ ]] && kill $1 2>/dev/null || true
        for i in ${children[$1]:=}; do kill_all $i; done
    }

    kill_all $1
}

nxf_mktemp() {
    local base=${1:-/tmp}
    mkdir -p "$base"
    if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
    else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
    fi
}

nxf_fs_copy() {
  local source=$1
  local target=$2
  local basedir=$(dirname $1)
  mkdir -p $target/$basedir
  cp -fRL $source $target/$basedir
}

nxf_fs_move() {
  local source=$1
  local target=$2
  local basedir=$(dirname $1)
  mkdir -p $target/$basedir
  mv -f $source $target/$basedir
}

nxf_fs_rsync() {
  rsync -rRl $1 $2
}

nxf_fs_rclone() {
  rclone copyto $1 $2/$1
}

nxf_fs_fcp() {
  fcp $1 $2/$1
}

on_exit() {
    exit_status=${nxf_main_ret:=$?}
    printf -- $exit_status > {{folder}}/.exitcode
    set +u
    exit $exit_status
}

on_term() {
    set +e
    [[ "$pid" ]] && nxf_kill $pid
}

nxf_launch() {
    /bin/bash -ue {{folder}}/.command.sh
}

nxf_stage() {
    true
}

nxf_unstage() {
    true
    [[ ${nxf_main_ret:=0} != 0 ]] && return
}

nxf_main() {
    trap on_exit EXIT
    trap on_term TERM INT USR2
    trap '' USR1

    [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
    NXF_SCRATCH=''
    [[ $NXF_DEBUG > 0 ]] && nxf_env
    touch {{folder}}/.command.begin
    set +u
    # beforeScript directive
    echo Before
    # conda environment
    source $(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", $2); print $2 }')/bin/activate /conda/env/path
    # spack environment
    spack env activate -d /spack/env/path
    set -u
    [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
    export NXF_TASK_WORKDIR="$PWD"
    nxf_stage

    set +e
    (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
    pid=$!
    wait $pid || nxf_main_ret=$?
    nxf_unstage
}

$NXF_ENTRY
