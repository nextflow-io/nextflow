/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.processor

import groovy.transform.PackageScope
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
import nextflow.util.ArrayTuple

import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger
/**
 * Manages the offload of publishing outputs to Nextflow tasks
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
class PublishOffloadManager {
    Map<Integer, ArrayTuple> runningPublications= new HashMap<Integer, ArrayTuple>(10)
    static final Map SUPPORTED_SCHEMES = [awsbatch:['s3'], local:['file']]
    static final String S5CMD_CONTAINER = 'jorgeejarquea/s5cmd_aws:0.0.1'
    static final String PUBLISH_FUNCTION = 'nxf_publish'
    private Session session
    private PublishTaskProcessor publishProcessor
    private List<String> commands = new LinkedList<String>()
    private boolean closed = false
    private PublishOffloadConfig config
    /**
     * Unique offloaded index number
     */
    final protected AtomicInteger indexCount = new AtomicInteger()

    PublishOffloadManager(Session session, PublishOffloadConfig config) {
        this.session = session
        this.config = config
    }
    @PackageScope
    TaskProcessor getPublishProcessor(){ publishProcessor }

    synchronized void init(){
        //Try with template
        BodyDef body = new BodyDef({
            def file = PublishOffloadManager.class.getResource('copy-group-template.sh')
            final engine = new TaskTemplateEngine()
            engine.setPlaceholder('!' as char)
            engine.eval(file.text, delegate)
            return engine.result
        },'publish file process', 'script')
        this.publishProcessor = createProcessor( "publish_process", body )
        log.debug("Publish Offload Manager initialized.")
    }

    private boolean checkOffload(Path source, Path destination, String executor){
        return this.config.enable && source.scheme in SUPPORTED_SCHEMES[executor] && destination.scheme in SUPPORTED_SCHEMES[executor]
    }

    private void invokeProcessor(inputValue) {
        final params = new TaskStartParams(TaskId.next(), publishProcessor.indexCount.incrementAndGet())
        final values = new ArrayList(2)
        values[0] = inputValue
        values[1] = this.config.maxParallel
        final args = new ArrayList(2)
        args[0] = params
        args[1] = values
        publishProcessor.invokeTask(args.toArray())
    }

    private synchronized boolean tryOffload(String command, Path origin, Path destination, PublishRetryConfig retryConfig, boolean failOnError){
        if (checkOffload(origin, destination, publishProcessor.executor.name)) {
            final id = indexCount.incrementAndGet()
            runningPublications.put(id, Nextflow.tuple(origin, destination, failOnError))
            commands.add(generateExecutionCommand(id, command, origin, destination, retryConfig))
            if (commands.size() == this.config.batchSize){
                invokeProcessor(commands.join(";"))
                commands.clear()
            }
            return true
        }
        return false
    }

    private isFusionEnabled(){
        return FusionHelper.isFusionEnabled(session) && config.useFusion
    }

    private useS5cmd(){
        return ( (!isFusionEnabled()) && (ExecutorFactory.getDefaultExecutorName(session) == 'awsbatch') )
    }

    private String generateExecutionCommand(Integer id, String command, Path origin, Path destination, PublishRetryConfig retryConfig){
        return "$PUBLISH_FUNCTION ${retryConfig.maxAttempts} ${retryConfig.delay.toMillis()} ${retryConfig.jitter} ${retryConfig.maxDelay.toMillis()} " +
                "$id $command ${convertFilePath(origin)} ${convertFilePath(destination)}"
    }

    private String convertFilePath(Path path) {
        if (isFusionEnabled()) {
            return FusionHelper.toContainerMount(path)
        } else {
            return FilesEx.toUriString(path)
        }
    }

    boolean tryMoveOffload(Path origin, Path destination, PublishRetryConfig retryConfig, boolean failOnError) {
        String command = 'mv'
        if ( useS5cmd() ) {
            command = 's5cmd mv'
        }
        tryOffload(command, origin, destination, retryConfig, failOnError)
    }

    boolean tryCopyOffload(Path origin, Path destination, PublishRetryConfig retryConfig, boolean failOnError) {
        String command = 'cp'
        if ( useS5cmd() ) {
            command = 's5cmd cp'
        }
        tryOffload(command, origin, destination, retryConfig, failOnError)
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
            processConfig.put('container', S5CMD_CONTAINER)
        }
        processConfig._in_val('executions')
        processConfig._in_val('max_parallel')
        processConfig._out_stdout()
        if ( !body )
            throw new IllegalArgumentException("Missing script in the specified process block -- make sure it terminates with the script string to be executed")

        // -- apply settings from config file to process config
        processConfig.applyConfig((Map)session.config.process, name, name, name)

        // -- get the executor for the given process config
        final execObj = session.executorFactory.getExecutor(name, processConfig, body, session)

        // -- create processor class
        new PublishTaskProcessor( name, execObj, session, session.script, processConfig, body, this )
    }

    synchronized void close() {
        closed=true
        if ( commands.size() ){
            invokeProcessor(commands.join(";"))
            commands.clear()
        }

    }
}

@Slf4j
class PublishTaskProcessor extends TaskProcessor{
    PublishOffloadManager manager

    PublishTaskProcessor(String name, Executor executor, Session session, BaseScript baseScript, ProcessConfig processConfig, BodyDef bodyDef, PublishOffloadManager manager) {
        super(name, executor, session, baseScript, processConfig, bodyDef)
        this.manager = manager
    }

    @Override
    void finalizeTask0(TaskRun task) {
        if( task.outputs.size() == 1 ){
            def value = task.outputs.values().first()
            value = value instanceof Path ? value.text : value?.toString()
            for ( String finishedCopy : value.split('\n') ){
                final result = finishedCopy.split(":")
                if (result.size() == 2) {
                    final id = result[0] as Integer
                    final tuple = manager.runningPublications.remove(id)
                    final exitCode = result[1] as Integer
                    if( exitCode == 0 ){
                        session.notifyFilePublish((Path) tuple.get(0), (Path) tuple.get(1))
                    } else {
                        if (tuple.get(2) as Boolean) {
                            log.error("Publication of file ${tuple.get(0)} -> ${tuple.get(1)} failed.")
                        } else {
                            log.warn("Publication of file ${tuple.get(0)} -> ${tuple.get(1)} failed.")
                        }
                        printPublishTaskError(task)
                    }
                }
            }
        } else {
            log.error("Incorrect number of outputs in the publish task")
        }
    }

    private void printPublishTaskError(TaskRun task){
        final List<String> message = []
        final max = 50
        final lines = task.dumpStderr(max)
        message << "Executed publish task:"
        if( lines ) {
                message << "\nCommand error:"
                for( String it : lines ) {
                    message << "  ${stripWorkDir(it, task.workDir)}"
                }
        }
        if( task?.workDir )
            message << "\nWork dir:\n  ${task.workDirStr}"
        log.debug(message.join('\n'))
    }
}
