package nextflow.cloud.gce.pipelines

import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean

import java.nio.file.Path

@Slf4j
@CompileStatic
class GooglePipelinesFileCopyStrategy extends SimpleFileCopyStrategy {

    GooglePipelinesTaskHandler handler
    TaskBean task

    GooglePipelinesFileCopyStrategy(TaskBean bean, GooglePipelinesTaskHandler handler) {
        super(bean)
        this.handler = handler
        this.task = bean
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {

        def stagingCommands = inputFiles.collect { stageName, storePath ->
            def absStorePath = storePath.toAbsolutePath()
            "gsutil -m  -q cp -P -c -r ${absStorePath.toUriString()} ${task.workDir}/${absStorePath.toString().endsWith(stageName) ? "" : stageName} || true".toString()
        }

        log.debug "[GOOGLE PIPELINE] Constructed the following file copy staging commands: $stagingCommands"

        handler.stagingCommands.addAll(stagingCommands)

        //Insert this comment into the task run script to note that the staging is done differently
        "# Google pipeline staging is done in a pipeline action step that is run prior to the main pipeline action"
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {

        def unstagingCommands = outputFiles.collect {
            "gsutil -m -q cp -P -r -c ${task.workDir}/$it ${task.workDir.toUriString()} || true".toString()
        }

        log.debug "[GOOGLE PIPELINE] Constructed the following file copy staging commands: $unstagingCommands"

        handler.unstagingCommands.addAll(unstagingCommands)

        //Insert this comment into the task run script to note that the unstaging is done differently
        ": # Google pipeline unstaging is done in a pipeline action step that is run after the main pipeline action"
    }

    //Although it seems like this is used as a general "touch" mechanism it's actually just used to create a .begin in the working directory
    @Override
    String touchFile(Path file) {
        if (file in CloudStoragePath) {
            handler.stagingCommands << "echo start | gsutil -q cp  -c - ${file.toUriString()} || true".toString()
            "# Google pipeline touchFile is done in a pipeline action step that is run prior to the main pipeline action"
        } else
            return super.touchFile(file)
    }
}