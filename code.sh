#!/bin/bash
set -e

# Main entry point for this app.
main() {
  nextflow -debug nextflow -c ${DX_FS_ROOT}/nextflow.config ${DX_FS_ROOT}/output_test.nf
}

# Entry point for parallel subtasks.
process() {
  echo "Starting PROCESS subtask ${taskName}"

  # Declare and set an associative array for inputs
  # Declare and set an array for outputs
  declare -A inputs="${inputs}"
  declare -a outputs="${outputs}"


  # Download all the files specified in "inputs"
  for i in "${!inputs[@]}";
      do
        input_file_id=$(dx-jobutil-parse-link "${inputs[$i]}") ;
        dx download "$input_file_id" --no-progress  ;
      done


  # Download the script file
  input_file_id=$(dx-jobutil-parse-link "${taskScript}")
  dx download "$input_file_id" -o task_script --no-progress


  # Change permissions to task_script to execute
  #printf "$taskEnv" > task_env
  chmod +x task_script


  # Execution of the task's script
  set +e
  #(source task_env; ./task_script) > my_output_file
  ./task_script
  exitcode=$?
  set -e


  # Set as output the exit code of the script's execution
  dx-jobutil-add-output exit_code "$exitcode" --class string


  # Upload and set as outputs the files which matches the
  # names or structures of the files declared in "outputs[]"
  for item in "${outputs[@]}"; do
    for name in `ls $item 2>/dev/null`; do
          output_file_id=`dx upload "${name}" --brief --no-progress` ;
          dx-jobutil-add-output "${name}" "$output_file_id" --class string ;
    done
  done

  echo "Finished PROCESS subtask ${taskName}"
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
