package nextflow.cloud.gce.pipelines

import groovy.transform.CompileStatic
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

import java.nio.file.Path

/**
 * Implements BASH launcher script for Google Pipeline
 */
//TODO: Need to use the scratch directory "more cleaner"
@CompileStatic
class GooglePipelinesScriptLauncher extends BashWrapperBuilder {

    GooglePipelinesScriptLauncher(TaskBean bean, GooglePipelinesTaskHandler handler) {
        super(bean, new GooglePipelinesFileCopyStrategy(bean, handler))

        // enable the copying of output file to the GS work dir
        scratch = "/work/scratch/"


        // include task script as an input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_SCRIPT] = bean.workDir.resolve(TaskRun.CMD_SCRIPT)
        // include the wrapper script as in input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_RUN] = bean.workDir.resolve(TaskRun.CMD_RUN)
        // add the wrapper file when stats are enabled
        if (bean.statsEnabled) {
            bean.inputFiles[TaskRun.CMD_STUB] = bean.workDir.resolve(TaskRun.CMD_STUB)
        }
        // include task stdin file
        if (bean.input != null) {
            bean.inputFiles[TaskRun.CMD_INFILE] = bean.workDir.resolve(TaskRun.CMD_INFILE)
        }
    }

    @Override
    protected String getScratchDirectoryCommand() {
        "NXF_SCRATCH=\$(mkdir $scratch)"
    }

    @Override
    String copyFile(String name, Path target) {
        "echo 'Google Pipelines file staging/unstaging happens in pre/post actions'"
    }
}