#!/bin/bash
#BSUB -x 1
#BSUB -y 2
# NEXTFLOW TASK: Hello 1
set -e
set -u
NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
NXF_ENTRY=${1:-nxf_main}


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
    if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
    else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
    fi
}

on_exit() {
    exit_status=${nxf_main_ret:=$?}
    printf $exit_status > {{folder}}/.exitcode
    set +u
    [[ "$tee1" ]] && kill $tee1 2>/dev/null
    [[ "$tee2" ]] && kill $tee2 2>/dev/null
    [[ "$ctmp" ]] && rm -rf $ctmp || true
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
    trap on_term TERM INT USR1 USR2

    [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
    NXF_SCRATCH=''
    [[ $NXF_DEBUG > 0 ]] && nxf_env
    touch {{folder}}/.command.begin
    set +u
    # beforeScript directive
    echo Before
    # conda environment
    source $(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", $2); print $2 }')/bin/activate /conda/env/path
    set -u
    [[ $NXF_SCRATCH ]] && echo "nxf-scratch-dir $HOSTNAME:$NXF_SCRATCH" && cd $NXF_SCRATCH
    nxf_stage

    set +e
    local ctmp=$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp $TMPDIR)
    local cout=$ctmp/.command.out; mkfifo $cout
    local cerr=$ctmp/.command.err; mkfifo $cerr
    tee .command.out < $cout &
    tee1=$!
    tee .command.err < $cerr >&2 &
    tee2=$!
    ( nxf_launch ) >$cout 2>$cerr &
    pid=$!
    wait $pid || nxf_main_ret=$?
    wait $tee1 $tee2
    nxf_unstage
}

$NXF_ENTRY
