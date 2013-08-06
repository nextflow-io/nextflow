#!/bin/bash
set -e

# Main entry point for this app.
main() {
    #echo 'Install Java 6 runtime'
    #apt-get install -y openjdk-6-jre-headless > /dev/null

    # Download the script file
    script_id=$(dx-jobutil-parse-link "${script}")
    dx download ${script_id} -o script.nf

    # Launch it !
    set +e
    nextflow -log $PWD\nextflow.log -c ${DX_FS_ROOT}/nextflow.config script.nf
    exit_status=$?
    set -e

    echo "nextflow exitstatus > ${exit_status}"
    if [ "$exit_status" -ne "0" ]; then cat $PWD\nextflow.log; fi

    # Returns the nextflow log
    log_file_id=$(dx upload $PWD\nextflow.log --brief)
    dx-jobutil-add-output 'log' "$log_file_id" --class=file

    exit ${exit_status}
}

# Entry point for parallel sub-tasks.
process() {
    echo "Starting PROCESS subtask ${taskName}"

    # Download all the files specified in "inputs"
    for input in "${inputs[@]}"
    do
        input_file_id=$(dx-jobutil-parse-link "${input}") ;
        dx download "$input_file_id" --no-progress  ;
    done

    # Download the script file
    input_file_id=$(dx-jobutil-parse-link "${taskScript}")
    dx download "$input_file_id" -o .command.sh --no-progress

    # Change permissions to '.command.sh' to execute
    # printf "$taskEnv" > task_env
    chmod +x .command.sh

    # Execution of the task's script
    set +e
    if [ -z "${taskInput}" ]; then
      ( ./.command.sh ) > .command.out
    else
      # Download the input file and save to .command.in
      input_taskfile_id=$(dx-jobutil-parse-link "${taskInput}")  ;
      dx download "$input_taskfile_id" -o .command.in --no-progress  ;
      ( cat .command.in | ./.command.sh ) > .command.out  ;
    fi
    exit_status=$?
    set -e

    # Set as output the exit code of the script's execution
    dx-jobutil-add-output exit_code "$exit_status" --class string

    # Upload the command result
    if [ -f .command.out ]; then
        output_file_id=`dx upload ".command.out" --brief --no-progress` ;
        dx-jobutil-add-output ".command.out" "$output_file_id" --class string ;
    fi

    if [ "$exit_status" -ne "0" ]; then
        ls -la >&2
        exit
    fi

    # Upload and set as outputs the files which matches the
    # names or structures of the files declared in "outputs[]"
    for item in "${outputs[@]}"; do
        for name in `ls $item 2>/dev/null`; do
            output_file_id=`dx upload "${name}" --brief --no-progress` ;
            dx-jobutil-add-output "${name}" "$output_file_id" --class string ;
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
