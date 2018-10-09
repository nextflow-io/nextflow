package nextflow.cloud.gce.pipelines

import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.executor.SupportedScriptTypes
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration

@Slf4j
@SupportedScriptTypes(ScriptType.SCRIPTLET)
class GooglePipelinesExecutor extends Executor {


    @PackageScope
    static Genomics genomicsClient

    static GooglePipelinesConfiguration pipelineConfig

    @Override
    final boolean isContainerNative() {
        return true
    }

    @Override
    void register() {
        super.register()

        pipelineConfig = validateConfiguration()
        log.debug "[GOOGLE PIPELINE] Pipeline config: $pipelineConfig"

        genomicsClient = GooglePipelinesHelper.createGenomicClient()

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

        new GooglePipelinesConfiguration(
                session.config.navigate("gce.project") as String,
                session.config.navigate("gce.zone") as String,
                session.config.navigate("cloud.instanceType") as String,
                session.config.navigate("cloud.preemptible") as boolean
        )
    }
}