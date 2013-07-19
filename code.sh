#!/bin/bash
set -e

# Main entry point for this app.
main() {
  nextflow -debug nextflow -c ${DX_FS_ROOT}/nextflow.config ${DX_FS_ROOT}/test.nf
}

# Entry point for parallel subtasks.
process() {
  echo "Starting PROCESS subtask ${task_name}"

  # Download the input files.
  input_file_id=$(dx-jobutil-parse-link "$taskScript")
  dx download "$input_file_id" -o task_script --no-progress

  printf "$taskEnv" > task_env
  chmod +x task_script

  set +e
  (source task_env; ./task_script) > output_file
  exitcode=$?
  set -e

  # Upload the output files and add the output.
  output_file_id=`dx upload output_file --brief --no-progress`
  dx-jobutil-add-output exit_code "$exitcode" --class string
  dx-jobutil-add-output output_file "$output_file_id" --class file

  echo "Finished PROCESS subtask ${task_name}"
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
