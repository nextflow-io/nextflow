#!/bin/bash
set -e

# Main entry point for this app.
main() {
    # log file name
    LOG=${name}.log
    export NXF_HOME='/opt/nextflow'
    export NXF_PACKAGE='dx'

    # Launch it !
    set +e
    bash /usr/bin/nextflow \
         -log $PWD/$LOG \
         run $script \
         -process.executor dnanexus \
         -cache false \
         -qs $queue \
         -w dxfs://${DX_WORKSPACE_ID}:/workspace/ \
         $params
    exit_status=$?
    set -e

    echo "nextflow exitstatus > ${exit_status}"
    if [ "$exit_status" -ne "0" ]; then tail -n 100 $PWD/$LOG >&2; fi

    # upload the nextflow log
    log_file_id=$(dx upload $PWD/$LOG --brief --path ${DX_PROJECT_CONTEXT_ID}:$LOG )
    echo "$LOG > $log_file_id"

    # report the application exit status
    exit ${exit_status}
}

# Entry point for parallel sub-tasks.
process() {
    echo "Sys info: cpus `nproc`; mem `cat /proc/meminfo | grep MemTotal | awk '{print $2$3}'`"
    echo "Starting: ${task_name}"

    # stage scripts and input files to current folder
    dx download --no-progress $task_script -o .command.sh
    [ ! -z "$task_env" ] && dx download --no-progress $task_env -o .command.env
    [ ! -z "$task_input" ] && dx download --no-progress $task_input -o .command.in
    [ ! -z "$stage_inputs" ] && eval "$stage_inputs"
    [ -f .command.env ] && source .command.env

    # Change permissions to '.command.sh' to execute
    chmod +x .command.sh

    # Execution of the task's script
    set +e
    if [ -f .command.in ]; then
      ( ./.command.sh &> .command.out < .command.in )
    else
      ( ./.command.sh &> .command.out )
    fi
    exit_status=$?
    set -e

    # Set as output the exit code of the script's execution
    dx-jobutil-add-output exit_code "$exit_status" --class=int

    # Upload and set as outputs the files which matches the
    # names or structures of the files declared in "outputs[]"
    dx upload .command.out --path ${target_dir}/.command.out --brief --no-progress --wait
    for item in "${output_files[@]}"; do
        for name in `ls $item 2>/dev/null`; do
            dx upload -r $name --path "${target_dir}/$name" --brief --no-progress --wait
        done
    done

    if [ "$exit_status" -ne "0" ]; then
        ls -la >&2
        exit
    fi

}
