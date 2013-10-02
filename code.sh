#!/bin/bash
set -e

# Main entry point for this app.
main() {
    #echo 'Install Java 6 runtime'
    #apt-get install -y openjdk-6-jre-headless > /dev/null

    # create workspace path
    dx mkdir -p /cloud/workspace
    echo CLASSPATH: $CLASSPATH

    # Download the script file
    script_id=$(dx-jobutil-parse-link "${script}")
    dx download ${script_id} -o script.nf

    # Launch it !
    set +e
    nextflow -log $PWD\nextflow.log -c ${DX_FS_ROOT}/nextflow.config -w dxfs:///cloud/workspace/ script.nf
    exit_status=$?
    set -e

    echo "nextflow exitstatus > ${exit_status}"
    if [ "$exit_status" -ne "0" ]; then tail -n 100 $PWD\nextflow.log >&2; fi

    # upload the nextflow log
    log_file_id=$(dx upload $PWD\nextflow.log --brief --path ${DX_PROJECT_CONTEXT_ID}:nextflow.log )
    echo "nextflow.log > $log_file_id"

    # report the application exit status
    exit ${exit_status}
}

# Entry point for parallel sub-tasks.
process() {
    echo "Starting PROCESS subtask ${taskName}"
    set -x

    # stage input files to current folder
    dx download --no-progress $task_script -o .command.sh
    [ ! -z $task_env ] && dx download --no-progress $task_env -o .command.env
    [ ! -z $task_input ] && download --no-progress $task_input -o .command.in

    # source the env file
    if [ -f .command.env ]; then
    source .command.env
    fi

    # Change permissions to '.command.sh' to execute
    chmod +x .command.sh

    # Execution of the task's script
    set +e
    if [ -f .command.in ]; then
      ( ./.command.sh > .command.out < .command.in )
    else
      ( ./.command.sh > .command.out )
    fi
    exit_status=$?
    set -e

    # Set as output the exit code of the script's execution
    dx-jobutil-add-output exit_code "$exit_status" --class int

    if [ "$exit_status" -ne "0" ]; then
        ls -la >&2
        exit
    fi

    # Upload and set as outputs the files which matches the
    # names or structures of the files declared in "outputs[]"
    PRJ=$(dx describe $task_script | grep Project | awk '{print $2}')
    TARGET=$(dx describe $task_script | grep Folder | awk '{print $2}')
    dx upload .command.out --path $PRJ:$TARGET/.command.out --brief --no-progress --wait
    for item in "${outputs[@]}"; do
        for name in `ls $item 2>/dev/null`; do
            dx upload $name --path "$PRJ:$(dirname $TARGET)/$name" --brief --no-progress --wait
        done
    done

}


# This entry point is run after all 'process' tasks have finished.
postprocess() {
  echo "Starting POSTPROCESS subtask"

  # Replace this with whatever work is needed to combine the work from the
  # PROCESS stages to make the final result.
  #
  # In this case we just concatenate all the files, so we end up with the
  # original input file repeated numSubtasks times.
  for process_output in "${process_outputs[@]}"
  do
    process_output_id=$(dx-jobutil-parse-link "$process_output")
    dx cat $process_output_id >> combined_output
  done

  # Upload the output files and add the output.
  combined_output_id=`dx upload combined_output --brief --no-progress`
  dx-jobutil-add-output combined_output "$combined_output_id"
}
