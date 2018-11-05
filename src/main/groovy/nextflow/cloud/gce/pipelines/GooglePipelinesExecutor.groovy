package nextflow.cloud.gce.pipelines


import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Nextflow
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.executor.SupportedScriptTypes
import nextflow.extension.FilesEx
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration

import java.nio.file.Path

@Slf4j
@SupportedScriptTypes(ScriptType.SCRIPTLET)
@CompileStatic
class GooglePipelinesExecutor extends Executor {

    static GooglePipelinesConfiguration pipelineConfig
    final GooglePipelinesHelper helper

    /*
     * Used for testing purposes
     */
    GooglePipelinesExecutor(GooglePipelinesHelper helper) {
        this.helper = helper
    }

    GooglePipelinesExecutor() {
        this.helper = new GooglePipelinesHelper()
    }

    @Override
    final boolean isContainerNative() {
        return true
    }

    @Override
    void register() {
        super.register()

        pipelineConfig = validateConfiguration()

        log.debug "[GOOGLE PIPELINE] Finished registration for executor $name"
    }

    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 1000, Duration.of('10 sec'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new GooglePipelinesTaskHandler(task, this, pipelineConfig)
    }

    GooglePipelinesConfiguration validateConfiguration() {

        //Make sure that the workdir is a GS Bucket
        if (!(session.workDir instanceof CloudStoragePath)) {
            session.abort()
            throw new AbortOperationException("When using `$name` executor a GCE bucket must be provided as a working directory -- Add the option `-w gs://<your-bucket/path>` to your run command line or specify a workDir in your config file.")
        }

        //Check for the existence of all required configuration for our executor
        def requiredConfig = ["gce.project", "gce.zone"]

        requiredConfig.each {
            if (!session.config.navigate(it)) {
                session.abort()
                throw new AbortOperationException("Required config value '$it' for executor $name is not defined. Please add it to your process or nextflow configuration file.")
            }
        }

        def path = session.config.navigate('env.PATH')
        if( path ) {
            log.warn "Environment PATH defined in config file is ignored by AWS Batch executor"
        }

        /*
         * upload local binaries
         */
        def disableBinDir = session.getExecConfigProp(name, 'disableRemoteBinDir', false)
        Path remoteBinDir = null
        if( session.binDir && !disableBinDir ) {
            def cloudPath = Nextflow.tempDir() + "/" //need the ending slash to mark it as a directory since it doesn't exist yet
            log.info "Uploading local `bin` scripts folder to ${cloudPath.toUriString()}bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, cloudPath)
        }

        return new GooglePipelinesConfiguration(
                session.config.navigate("gce.project") as String,
                session.config.navigate("gce.zone") as String,
                session.config.navigate("cloud.instanceType") as String,
                remoteBinDir,
                session.config.navigate("cloud.preemptible") as boolean
        )
    }
}