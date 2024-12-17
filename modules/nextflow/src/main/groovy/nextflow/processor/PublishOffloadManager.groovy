package nextflow.processor

import groovy.util.logging.Slf4j
import nextflow.Nextflow
import nextflow.Session
import nextflow.executor.Executor
import nextflow.executor.ExecutorFactory
import nextflow.extension.FilesEx
import nextflow.fusion.FusionHelper
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.TokenValCall
import nextflow.util.ArrayTuple

import java.nio.file.Path

@Slf4j
class PublishOffloadManager {
    Map<TaskRun, ArrayTuple> runningPublications= new HashMap<TaskRun, ArrayTuple>(10)
    static final Map SUPPORTED_SCHEMES = [awsbatch:['s3'], local:['file']]
    static final String S5CMD_CONTAINER = 'jorgeejarquea/s5cmd_aws:0.0.1'
    Session session
    PublishTaskProcessor copyProcessor
    PublishTaskProcessor moveProcessor

    PublishOffloadManager(Session session) {
        this.session = session
    }

    void init(){

        if (useS5cmd()){
            this.copyProcessor = createProcessor( "publish_dir_copy_process", new BodyDef({"s5cmd cp $source $target"},'copy data process') )
            this.moveProcessor = createProcessor( "publish_dir_move_process", new BodyDef({"s5cmd mv $source $target"},'move data process') )
        } else {
            this.copyProcessor = createProcessor( "publish_dir_copy_process", new BodyDef({"cp $source $target"},'copy data process') )
            this.moveProcessor = createProcessor( "publish_dir_move_process", new BodyDef({"mv $source $target"},'move data process') )
        }

    }

    private boolean checkOffload(Path source, Path destination, String executor){
        return session.publishOffload && source.scheme in SUPPORTED_SCHEMES[executor] && destination.scheme in SUPPORTED_SCHEMES[executor];
    }

    private synchronized boolean tryInvokeProcessor(TaskProcessor processor, Path origin, Path destination){
        if (checkOffload(origin, destination, processor.executor.name)) {
            final params = new TaskStartParams(TaskId.next(), processor.indexCount.incrementAndGet())
            final values = new ArrayList(1)
            log.debug("Creating task for file publication: ${origin.toUri().toString()} -> ${destination.toUri().toString()} " )
            values[0] = generateFileValues(origin, destination)
            final args = new ArrayList(2)
            args[0] = params
            args[1] = values
            assert args.size() == 2
            processor.invokeTask(args.toArray())
            runningPublications.put(processor.currentTask.get(), Nextflow.tuple(origin, destination))
            return true
        }
        return false
    }

    private isFusionEnabled(){
        return FusionHelper.isFusionEnabled(session)
    }

    private useS5cmd(){
        return ( (!isFusionEnabled()) && (ExecutorFactory.getDefaultExecutorName(session) == 'awsbatch') )
    }

    private ArrayTuple<String> generateFileValues(Path origin, Path destination){
        if ( isFusionEnabled() ){
            Nextflow.tuple(FusionHelper.toContainerMount(origin), FusionHelper.toContainerMount(destination))
        } else {
            Nextflow.tuple(FilesEx.toUriString(origin), FilesEx.toUriString(destination))
        }
    }

    boolean tryMoveOffload(Path origin, Path destination) {
        tryInvokeProcessor(moveProcessor, origin, destination)
    }

    boolean tryCopyOffload(Path origin, Path destination) {
        tryInvokeProcessor(copyProcessor, origin, destination)
    }

    private PublishTaskProcessor createProcessor( String name, BodyDef body){
        assert body
        assert session.script
        log.debug("Creating processor $name")
        // -- the config object
        final processConfig = new ProcessConfig(session.script, name)
        // Invoke the code block which will return the script closure to the executed.
        // As side effect will set all the property declarations in the 'taskConfig' object.
        if (useS5cmd()) {
            processConfig.put('container', S5CMD_CONTAINER);
        }
        processConfig._in_tuple(new TokenValCall('source'), new TokenValCall('target'))

        if ( !body )
            throw new IllegalArgumentException("Missing script in the specified process block -- make sure it terminates with the script string to be executed")

        // -- apply settings from config file to process config
        processConfig.applyConfig((Map)session.config.process, name, name, name)

        // -- get the executor for the given process config
        final execObj = session.executorFactory.getExecutor(name, processConfig, body, session)

        // -- create processor class
        new PublishTaskProcessor( name, execObj, session, session.script, processConfig, body, this )
    }

}

class PublishTaskProcessor extends TaskProcessor{

    PublishOffloadManager manager

    PublishTaskProcessor(String name, Executor executor, Session session, BaseScript baseScript, ProcessConfig processConfig, BodyDef bodyDef, PublishOffloadManager manager) {
        super(name, executor, session, baseScript, processConfig, bodyDef)
        this.manager = manager
    }

    @Override
    void finalizeTask0(TaskRun task){
        final tuple = manager.runningPublications.remove(task)
        session.notifyFilePublish((Path)tuple.get(0), (Path)tuple.get(1))
    }
}
