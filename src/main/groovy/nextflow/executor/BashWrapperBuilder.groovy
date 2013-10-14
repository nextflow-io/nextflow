package nextflow.executor
import java.nio.file.Path

import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
/**
 * Builder to create the BASH script which is used to
 * wrap and launch the user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilder {

    final private TaskRun task

    def scratch

    def input

    Map<String,String> environment

    Closure<String> stagingScript

    Closure<String> unstagingScript

    String wrapperScript


    BashWrapperBuilder( TaskRun task )  {
        assert task
        assert task.workDirectory
        this.task = task
    }

    /* this constructor is used only for testing purpose */
    protected BashWrapperBuilder () { task = null }

    /**
     * @return The bash script fragment to change to the 'scratch' directory if it has been specified in the task configuration
     */
    protected String changeToScratchDirectory() {

        if( scratch == null || scratch == false ) {
            return null
        }

        /*
         * when 'scratch' is defined as a bool value
         * try to use the 'TMP' variable, if does not exist fallback to a tmp folder
         */
        if( scratch == true ) {
            return 'NF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NF_SCRATCH'
        }

        // convert to string for safety
        scratch = scratch.toString()

        // when it is defined by a variable, just use it
        if( scratch.startsWith('$') ) {
            return "NF_SCRATCH=$scratch && cd \$NF_SCRATCH"
        }

        if( scratch.toLowerCase() in ['ramdisk','ram-disk']) {
            return 'NF_SCRATCH=$(mktemp -d -p /dev/shm/) && cd $NF_SCRATCH'
        }


        return "NF_SCRATCH=\$(mktemp -d -p $scratch) && cd \$NF_SCRATCH"

    }

    /**
     * Build up the BASH wrapper script file which will launch the user provided script
     * @return The {@code Path} of the created wrapper script
     */

    Path build() {

        final scriptFile = task.getCmdScriptFile()
        final inputFile = task.getCmdInputFile()
        final environmentFile = task.getCmdEnvironmentFile()
        final startedFile = task.getCmdStartedFile()
        final outputFile = task.getCmdOutputFile()
        final exitedFile = task.getCmdExitFile()
        final wrapperFile = task.getCmdWrapperFile()

        /*
         * the script file
         */
        final taskScript = task.processor.normalizeScript(task.script.toString())
        scriptFile.text = taskScript

        /*
         * save the input when required
         */
        if( input != null ) {
            inputFile.text = input
        }

        /*
         * save the environment
         */
        if( environment ) {
            // create the *bash* environment script
            environmentFile.text = TaskProcessor.bashEnvironmentScript(environment)
        }

        /*
         * create a script wrapper which do the following
         * 1 - move the TMP directory provided by the sge/oge grid engine
         * 2 - pipe the input stream
         * 3 - launch the user script
         * 4 - un-stage e.g. copy back the result files to the working folder
         */

        def ENDL = '\n'
        def wrapper = new StringBuilder()
        wrapper << '#!/bin/bash -Eeu' << ENDL
        wrapper << 'trap onexit 1 2 3 15 ERR' << ENDL
        wrapper << 'function onexit() { local exit_status=${1:-$?}; printf $exit_status > ' << exitedFile.toString() << '; exit $exit_status; }' << ENDL
        wrapper << 'touch ' << startedFile.toString() << ENDL

        // source the environment
        wrapper << 'source ' << environmentFile.toString() << ENDL

        // whenever it has to change to the scratch directory
        def changeDir = changeToScratchDirectory()
        if( changeDir ) {
            wrapper << changeDir << ENDL
        }

        // staging input files when required
        def staging
        if( stagingScript && (staging=stagingScript.call()) ) {
            wrapper << staging << ENDL
        }

        // fetch the script interpreter
        final interpreter = task.processor.fetchInterpreter(taskScript)

        // execute the command script
        wrapper << '( '
        wrapper << interpreter << ' ' << scriptFile.toString()
        if( input != null ) wrapper << ' < ' << inputFile.toString()
        wrapper << ' ) &> ' << outputFile.toAbsolutePath() << ENDL

        def unstaging
        if( changeDir && unstagingScript && (unstaging=unstagingScript.call()) ) {
            wrapper << unstaging << ENDL
        }

        wrapper << 'onexit' << ENDL

        wrapperFile.text = wrapperScript = wrapper.toString()
        return wrapperFile
    }

}
