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
            "gsutil -m  -q cp -P -r ${absStorePath.toUriString()} ${task.workDir}/${absStorePath.toString().endsWith(stageName) ? "" : stageName}".toString()
        }

        log.debug "[GOOGLE PIPELINE] Constructed the following file copy staging commands: $stagingCommands"

        handler.stagingCommands.addAll(stagingCommands)

        //copy the remoteBinDir if it ia defined
        if(handler.pipelineConfiguration.remoteBinDir) {
            def createRemoteBinDir = "mkdir -p ${task.workDir}/nextflow-bin".toString()
            def remoteBinCopy = "gsutil -m -q cp -P -r ${handler.pipelineConfiguration.remoteBinDir.toUriString()}/* ${task.workDir}/nextflow-bin".toString()
            handler.stagingCommands << createRemoteBinDir << remoteBinCopy
        }

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

    /**
     * {@inheritDoc}
     */
    @Override
    String getEnvScript(Map environment, String wrapper) {
        if( wrapper )
            throw new IllegalArgumentException("Parameter `wrapHandler` not supported by ${this.class.simpleName}")

        final result = new StringBuilder()
        final copy = environment ? new HashMap<String,String>(environment) : Collections.<String,String>emptyMap()
        // remove any external PATH
        if( copy.containsKey('PATH') )
            copy.remove('PATH')
        // when a remote bin directory is provide managed it properly
        if( handler.pipelineConfiguration.remoteBinDir ) {
            result << "chmod +x $task.workDir/nextflow-bin/*\n"
            result << "export PATH=$task.workDir/nextflow-bin:\$PATH\n"
        }
        // finally render the environment
        final envSnippet = super.getEnvScript(copy,null)
        if( envSnippet )
            result << envSnippet
        return result.toString()
    }
}