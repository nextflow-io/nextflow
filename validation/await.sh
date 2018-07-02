#!/bin/bash
set -eu

TIMEOUT=3600    # terminate after one hour
KEEPALIVE=480   # time after which print an output to keep alive travis

# launch main command to await  
( "$@" | tee stdout.log ) &
pid=$!

# make sure to terminate target command on exit
trap "kill $pid" TERM INT USR1 USR2

# monitor task execution 
begin=$(date +%s)
size=0
while kill -0 $pid > /dev/null 2>&1; do 
    sleep 1

    # continue to await if there's a change 
    # in the stdout file 
    ### BSD current=$(stat -f%z stdout.log)
    current=$(stat -c '%s' stdout.log)
    if [[ $current != $size ]]; then
      size=$current
      begin=$(date +%s)
      continue
    fi

    # kill the execution if it's taking too 
    # much time without producing any output
    now=$(date +%s)
    delta=$((now-begin))
    if ((delta>=$TIMEOUT)); then
        echo Taking too long ... killing it!
        kill $pid
        exit 1
    elif ((delta>=$KEEPALIVE)); then    
        echo "[keep alive]"
    fi  
done