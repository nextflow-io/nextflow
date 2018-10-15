package nextflow.cloud.gce.pipelines

import groovy.transform.CompileStatic
import nextflow.executor.BashWrapperBuilder
import nextflow.file.FileHelper
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 * Implements BASH launcher script for Google Pipeline
 */
//TODO: This code is copied nearly 100% from AWS batch.  Need to grok it and rewrite it
//TODO: Scratch stuff and changing the targetDir are both unsuitable to trigger getUnstageOutputFilesScript
@CompileStatic
class GooglePipelinesScriptLauncher extends BashWrapperBuilder {

    GooglePipelinesScriptLauncher(TaskBean bean, GooglePipelinesTaskHandler handler) {
        super(bean, new GooglePipelinesFileCopyStrategy(bean, handler))

        // enable the copying of output file to the GS work dir
        //scratch = false

        //TODO: here lies a hackdragon, 'YARR!!!!
        bean.targetDir = FileHelper.asPath("/work")

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
}