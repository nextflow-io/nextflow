#!/bin/bash
set -e
#set -x
#mkdir -p /cloud; dx-mount --project-id $DX_PROJECT_CONTEXT_ID /cloud


# Main entry point for this app.
main() {
    #echo 'Install Java 6 runtime'
    #apt-get install -y openjdk-6-jre-headless > /dev/null

    # create workspace path
    dx mkdir -p /workspace

    # Launch it !
    set +e
    java -jar /usr/bin/nextflow.jar -log $PWD/nextflow.log -process.executor dnanexus -w dxfs:///workspace/ $script $params
    exit_status=$?
    set -e

    echo "nextflow exitstatus > ${exit_status}"
    if [ "$exit_status" -ne "0" ]; then tail -n 100 $PWD/nextflow.log >&2; fi

    # upload the nextflow log
    log_file_id=$(dx upload $PWD/nextflow.log --brief --path ${DX_PROJECT_CONTEXT_ID}:nextflow.log )
    echo "nextflow.log > $log_file_id"

    # report the application exit status
    exit ${exit_status}
}

# Entry point for parallel sub-tasks.
process() {
    echo "Starting PROCESS subtask ${task_name}"

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
    dx-jobutil-add-output exit_code "$exit_status" --class int

    # Upload and set as outputs the files which matches the
    # names or structures of the files declared in "outputs[]"
    PRJ=$(dx describe $task_script | grep Project | awk '{print $2}')
    TARGET=$(dx describe $task_script | grep Folder | awk '{print $2}')
    dx upload .command.out --path $PRJ:$TARGET/.command.out --brief --no-progress --wait
    for item in "${output_files[@]}"; do
        for name in `ls $item 2>/dev/null`; do
            dx upload -r $name --path "$PRJ:$TARGET/$name" --brief --no-progress --wait
        done
    done

    if [ "$exit_status" -ne "0" ]; then
        ls -la >&2
        exit
    fi

}
