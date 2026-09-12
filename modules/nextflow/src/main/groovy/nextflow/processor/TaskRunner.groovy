/*
 * Copyright 2013-2026, Seqera Labs
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

import static nextflow.processor.ErrorStrategy.*

import java.nio.file.Path
import java.util.concurrent.atomic.AtomicBoolean
import java.util.function.Consumer

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.CloudSpotTerminationException
import nextflow.exception.FailedGuardException
import nextflow.exception.IllegalArityException
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.exception.ProcessEvalException
import nextflow.exception.ProcessException
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessRetryableException
import nextflow.exception.ProcessSubmitTimeoutException
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.CachedTaskHandler
import nextflow.executor.Executor
import nextflow.executor.StoredTaskHandler
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.file.FilePorter
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ProcessConfigV1
import nextflow.script.ProcessConfigV2
import nextflow.script.ScriptType
import nextflow.script.dsl.Types
import nextflow.script.params.EachInParam
import nextflow.script.params.EnvInParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InParam
import nextflow.script.params.OutParam
import nextflow.script.params.StdInParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.TupleOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.ValueOutParam
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessTupleInput
import nextflow.script.types.Record
import nextflow.script.types.Tuple
import nextflow.trace.TraceRecord
import nextflow.util.HashBuilder
import nextflow.util.LockManager
import nextflow.util.RecordMap
/**
 * Implements the execution of a single task, that is the resolution of the task
 * inputs and outputs, the cache and store directory lookup, the delegation to
 * the executor and the task error strategy.
 *
 * A runner is scoped to a process, but it has no dependency on the dataflow
 * interface of the process: a task can be executed by invoking {@link #submit}
 * directly, without the process being part of a dataflow network. The completion
 * of a task is reported to the completion handler given at construction.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class TaskRunner {

    /**
     * The process for which this runner executes tasks. It provides the
     * process-scoped services (name, config, executor, session, script)
     * shared by all task executions.
     */
    private TaskProcessor processor

    /**
     * Handler invoked when a task execution has completed.
     */
    private Consumer<TaskRun> onComplete

    /**
     * Count the number of task failures for this process.
     */
    private volatile int errorCount

    /**
     * Make sure the error message is shown only the very first time across all processes.
     */
    @PackageScope
    static final AtomicBoolean errorShown = new AtomicBoolean()

    private static LockManager lockManager = new LockManager()

    TaskRunner(TaskProcessor processor, Consumer<TaskRun> onComplete) {
        this.processor = processor
        this.onComplete = onComplete
    }

    int getErrorCount() { errorCount }

    /*
     * Process-scoped services owned by the process.
     */
    private Session getSession() { processor.getSession() }

    private Executor getExecutor() { processor.getExecutor() }

    private ProcessConfig getConfig() { processor.getConfig() }

    private ProcessConfigV1 configV1() { getConfig() as ProcessConfigV1 }

    private ProcessConfigV2 configV2() { getConfig() as ProcessConfigV2 }

    private BodyDef getTaskBody() { processor.getTaskBody() }

    private ScriptType getScriptType() { processor.getScriptType() }

    private String getName() { processor.getName() }

    private String safeTaskName(TaskRun task) { processor.safeTaskName(task) }

    private void checkWarn(String msg, Map opts=null) { processor.checkWarn(msg, opts) }

    /**
     * Execute a task for the given input values.
     *
     * The task is executed asynchronously: this method returns once the task has
     * been submitted to the executor, or as soon as the task result has been
     * resolved from the store directory or the cache. The completion handler is
     * invoked when the task execution has completed.
     *
     * @param params
     *      The task id and index. They are supplied by the caller because the
     *      index determines the ordering of the task executions
     * @param values
     *      The list of task input values
     */
    void submit(TaskStartParams params, List values) {
        log.trace "Invoking task > $name with params=$params; values=$values"

        // -- create the task run instance
        final task = createTaskRun(params)
        // -- set the task instance as the current in this thread
        processor.currentTask.set(task)

        // -- map the inputs to a map and use to delegate closure values interpolation
        final foreignFiles = session.filePorter.newBatch(executor.getStageDir())

        resolveTaskInputs(task, values, foreignFiles)

        // verify that `when` guard, when specified, is satisfied
        if( !checkWhenGuard(task) )
            return

        // -- resolve the task command script
        task.resolve(taskBody)

        // -- verify if exists a stored result for this case,
        //    if true skip the execution and return the stored data
        if( checkStoredOutput(task) )
            return

        // -- download foreign files
        session.filePorter.transfer(foreignFiles)

        final hash = new TaskHasher(task).compute()
        checkCachedOrLaunchTask(task, hash, processor.isResumable())
    }

    /**
     * Create a new {@code TaskRun} instance, initializing the following properties :
     * <li>{@code TaskRun#id}
     * <li>{@code TaskRun#status}
     * <li>{@code TaskRun#index}
     * <li>{@code TaskRun#name}
     * <li>{@code TaskRun#process}
     *
     * @return The new newly created {@code TaskRun}
     */

    final protected TaskRun createTaskRun(TaskStartParams params) {
        final task = new TaskRun(
                id: params.id,
                index: params.index,
                processor: processor,
                type: scriptType,
                config: config.createTaskConfig(),
                context: new TaskContext(processor)
        )

        // setup config
        task.config.index = task.index
        task.config.process = task.processor.name
        task.config.executor = task.processor.executor.name

        if( config instanceof ProcessConfigV1 )
            initializeTaskRunV1(task)

        return task
    }

    private void initializeTaskRunV1(TaskRun task) {
        /*
         * initialize the inputs/outputs for this task instance
         */
        configV1().getInputs().each { InParam param ->
            if( param instanceof TupleInParam )
                param.inner.each { task.setInput(it)  }
            else if( param instanceof EachInParam )
                task.setInput(param.inner)
            else
                task.setInput(param)
        }

        configV1().getOutputs().each { OutParam param ->
            if( param instanceof TupleOutParam ) {
                param.inner.each { task.setOutput(it) }
            }
            else
                task.setOutput(param)
        }
    }

    final protected void resolveTaskInputs(TaskRun task, List values, FilePorter.Batch foreignFiles) {
        if( config instanceof ProcessConfigV1 )
            resolveTaskInputsV1(task, values, foreignFiles)
        else if( config instanceof ProcessConfigV2 )
            resolveTaskInputsV2(task, values, foreignFiles)
    }

    private void resolveTaskInputsV1(TaskRun task, List values, FilePorter.Batch foreignFiles) {

        // -- validate input lengths
        for( final param : configV1().getInputs().ofType(TupleInParam) ) {
            final value = values[param.index]
            final expected = param.inner.size()
            final actual = value instanceof Collection
                ? value.size()
                : (value instanceof Map ? value.size() : 1)

            if( actual != expected ) {
                final message = "Input tuple does not match tuple declaration in process `$name` -- offending value: $value"
                checkWarn(message, [firstOnly: true, cacheKey: processor])
            }
        }

        // -- add input params to task context
        final Map<FileInParam,?> fileParams = [:]

        task.inputs.keySet().each { InParam param ->

            // add the value to the task instance
            def val = param.decodeInputs(values)

            switch(param) {
                case ValueInParam:
                    task.context.put(param.name, val)
                    break

                case FileInParam:
                    fileParams[param] = val
                    return // <-- leave it, because we do not want to add this 'val' at this stage

                case StdInParam:
                    task.stdin = val
                    break

                case EnvInParam:
                    task.inputEnv.put(param.name, val?.toString())
                    break

                default:
                    throw new IllegalStateException("Unsupported input param type: ${param?.class?.simpleName}")
            }

            // add the value to the task instance context
            task.setInput(param, val)
        }

        // -- all file parameters are processed in a second pass
        //    so that we can use resolve the variables that eventually are in the file name
        final resolver = new TaskInputResolver(task, foreignFiles, executor)

        for( Map.Entry<FileInParam,?> entry : fileParams.entrySet() ) {
            final param = entry.getKey()
            final val = entry.getValue()
            final resolved = resolver.resolve(param, val)

            if( !param.isValidArity(resolved.size()) )
                throw new IllegalArityException("Incorrect number of input files for process `${safeTaskName(task)}` -- expected ${param.arity}, found ${resolved.size()}")

            // add the value to the task instance context
            task.setInput(param, resolved)
            task.inputFiles.addAll(resolved)
        }

        // -- set the delegate map as context in the task config
        //    so that lazy directives will be resolved against it
        task.config.context = task.context

        // check conflicting file names
        checkConflicts(task.inputFiles)
    }

    private void checkConflicts(List<FileHolder> allFiles) {
        final allNames = new HashMap<String,Integer>()
        for( final holder : allFiles ) {
            final num = allNames.getOrCreate(holder.stageName, 0) + 1
            allNames.put(holder.stageName, num)
        }

        final conflicts = allNames.findAll { name, num -> num > 1 }
        if( conflicts ) {
            log.debug("Process $name > collision check staging file names: $allNames")
            def message = "Process `$name` input file name collision -- There are multiple input files for each of the following file names: ${conflicts.keySet().join(', ')}"
            throw new ProcessUnrecoverableException(message)
        }
    }

    @CompileStatic
    private void resolveTaskInputsV2(TaskRun task, List values, FilePorter.Batch foreignFiles) {
        final declaredInputs = configV2().getInputs()
        final ctx = task.context

        // -- remove control input
        values = values.init()

        // -- validate input lengths
        if( declaredInputs.size() != values.size() )
            throw new ProcessUnrecoverableException("Process `$name` declares ${declaredInputs.size()} inputs but received ${values.size()} -- offending value: $values")

        // -- add input params to task context
        for( int i = 0; i < declaredInputs.getParams().size(); i++ ) {
            final param = declaredInputs.getParams()[i]
            final value = values[i]
            if( param instanceof ProcessTupleInput && param.getType() == Record.class )
                assignTaskRecordInput(task, param, value, i)
            else if( param instanceof ProcessTupleInput && param.getType() == Tuple.class )
                assignTaskTupleInput(task, param, value, i)
            else
                assignTaskInput(task, param, value, i)
        }

        // -- resolve environment vars
        for( final entry : declaredInputs.getEnv() ) {
            final value = ctx.resolveLazy(entry.value)
            task.inputEnv.put(entry.key, value?.toString())
        }

        // -- resolve stdin
        task.stdin = ctx.resolveLazy(declaredInputs.stdin)

        // -- resolve input files
        final resolver = new TaskInputResolver(task, foreignFiles, executor)
        final resolvedValues = new HashSet<>()

        for( final fileInput : declaredInputs.getFiles() ) {
            final value = fileInput.resolve(ctx)
            // allow nullable file inputs
            if( value == null )
                continue
            // user-defined file stagers take precedence over default file stagers
            if( resolvedValues.contains(value) )
                continue
            final resolved = resolver.resolve(fileInput, value)
            task.inputFiles.addAll(resolved)
            resolvedValues.add(value)
        }

        // -- normalize input values by replacing source paths with staged paths
        final Map<Path,FileHolder> holders = [:]
        for( final holder : task.inputFiles )
            holders.put(holder.getSourcePath(), holder)

        for( final param : declaredInputs.getParams() ) {
            if( param instanceof ProcessTupleInput ) {
                for( final innerParam : param.getComponents() ) {
                    final value = task.inputs[innerParam]
                    final normalizedValue = resolver.normalizeValue(value, holders)
                    task.context.put( innerParam.name, normalizedValue )
                }
            }
            else {
                final value = task.inputs[param]
                final normalizedValue = resolver.normalizeValue(value, holders)
                task.context.put( param.name, normalizedValue )
            }
        }

        // -- set the delegate map as context in the task config
        //    so that lazy directives will be resolved against it
        task.config.context = ctx
    }

    @CompileStatic
    private void assignTaskRecordInput(TaskRun task, ProcessTupleInput param, Object value, int index) {
        if( value == null && !param.optional ) {
            throw new ProcessUnrecoverableException("[${safeTaskName(task)}] input at index ${index} cannot be null")
        }
        if( value !instanceof RecordMap ) {
            throw new ProcessUnrecoverableException("[${safeTaskName(task)}] input at index ${index} expected a record but received: ${value} [${value.class.simpleName}]")
        }
        final recordParams = param.getComponents()
        final record = value as Map
        for( final recordParam : recordParams ) {
            assignTaskInput(task, recordParam, record[recordParam.getName()], index)
        }
    }

    @CompileStatic
    private void assignTaskTupleInput(TaskRun task, ProcessTupleInput param, Object value, int index) {
        if( value == null && !param.optional ) {
            throw new ProcessUnrecoverableException("[${safeTaskName(task)}] input at index ${index} cannot be null")
        }
        if( value !instanceof List ) {
            throw new ProcessUnrecoverableException("[${safeTaskName(task)}] input at index ${index} expected a tuple but received: ${value} [${value.class.simpleName}]")
        }
        final tupleParams = param.getComponents()
        final tupleValues = value as List
        if( tupleParams.size() != tupleValues.size() ) {
            throw new ProcessUnrecoverableException("[${safeTaskName(task)}] input at index ${index} expected a tuple with ${tupleParams.size()} elements but received ${tupleValues.size()} -- offending value: $tupleValues")
        }
        for( int i = 0; i < tupleParams.size(); i++ ) {
            assignTaskInput(task, tupleParams[i], tupleValues[i], index)
        }
    }

    @CompileStatic
    private void assignTaskInput(TaskRun task, ProcessInput param, Object value, int index) {
        if( value == null && !param.optional ) {
            throw new ProcessUnrecoverableException("[${safeTaskName(task)}] input at index ${index} cannot be null -- append `?` to the type annotation to mark it as nullable")
        }
        if( value != null ) {
            final expectedType = param.type
            final actualType = value.getClass()
            if( expectedType != null && !isAssignableFrom(expectedType, actualType) )
                log.warn "[${safeTaskName(task)}] invalid argument type at index ${index} -- expected a ${Types.getName(expectedType)} but got a ${Types.getName(actualType)}"
        }
        task.context.put(param.getName(), value)
        task.setInput(param, value)
    }

    private static boolean isAssignableFrom(Class targetType, Class sourceType) {
        // treat all record types as compatible
        // record types are validated at compile-time
        if( Record.class.isAssignableFrom(targetType) && Record.class.isAssignableFrom(sourceType) )
            return true
        return targetType.isAssignableFrom(sourceType)
    }

    /**
     * Check if exists a *storeDir* for the specified task. When if exists
     * and contains the expected result files, the process execution is skipped.
     *
     * @param task The task for which check the stored output
     * @return {@code true} when the folder exists and it contains the expected outputs,
     *      {@code false} otherwise
     */
    final boolean checkStoredOutput( TaskRun task ) {
        if( !task.config.storeDir ) {
            log.trace "[${safeTaskName(task)}] storeDir not set -- return false"
            return false
        }

        // -- when store path is set, only output params of type 'file' can be specified
        if( isInvalidStoreDir(task) ) {
            checkWarn "[${safeTaskName(task)}] storeDir can only be used with `val` and `path` outputs"
            return false
        }

        if( !task.config.getStoreDir().exists() ) {
            log.trace "[${safeTaskName(task)}] storeDir does not exist > ${task.config.storeDir} -- return false"
            // no folder -> no cached result
            return false
        }


        try {
            // -- expose task exit status to make accessible as output value
            task.config.exitStatus = TaskConfig.EXIT_ZERO
            // -- check if all output resources are available
            collectOutputs(task)
            log.info "[skipping] Stored process > ${safeTaskName(task)}"
            // set the exit code in to the task object
            task.exitStatus = TaskConfig.EXIT_ZERO
            task.cached = true
            session.notifyTaskCached(new StoredTaskHandler(task))

            // -- now bind the results
            onComplete.accept(task)
            return true
        }
        catch( MissingFileException | MissingValueException e ) {
            log.trace "[${safeTaskName(task)}] Missed storeDir > ${e.getMessage()} -- folder: ${task.config.storeDir}"
            task.exitStatus = Integer.MAX_VALUE
            task.workDir = null
            return false
        }
    }

    protected boolean isInvalidStoreDir(TaskRun task) {
        final ctx = task.context

        if( config instanceof ProcessConfigV1 ) {
            return task.getOutputs().keySet().any { param ->
                if( param instanceof ValueOutParam )
                    return false
                if( param instanceof FileOutParam )
                    return false
                return true
            }
        }

        if( config instanceof ProcessConfigV2 ) {
            return config.getOutputs().getFiles().isEmpty()
        }

        return false
    }

    /**
     * Try to check if exists a previously executed process result in the a cached folder. If it exists
     * use the that result and skip the process execution, otherwise the task is sumitted for execution.
     *
     * @param task
     *      The {@code TaskRun} instance to be executed
     * @param hash
     *      The unique {@code HashCode} for the given task inputs
     * @param script
     *      The script to be run (only when it's a merge task)
     * @return
     *      {@code false} when a cached result has been found and the execution has skipped,
     *      or {@code true} if the task has been submitted for execution
     *
     */
    @CompileStatic
    final protected void checkCachedOrLaunchTask( TaskRun task, HashCode hash, boolean shouldTryCache ) {

        int tries = task.failCount +1
        while( true ) {
            hash = HashBuilder.defaultHasher().putBytes(hash.asBytes()).putInt(tries).hash()

            Path resumeDir = null
            boolean exists = false
            try {
                final entry = session.cache.getTaskEntry(hash, processor)
                resumeDir = entry ? FileHelper.asPath(entry.trace.getWorkDir()) : null
                if( resumeDir )
                    exists = resumeDir.exists()

                log.trace "[${safeTaskName(task)}] Cacheable folder=${resumeDir?.toUriString()} -- exists=$exists; try=$tries; shouldTryCache=$shouldTryCache; entry=$entry"
                final cached = shouldTryCache && exists && entry.trace.isCompleted() && checkCachedOutput(task.clone(), resumeDir, hash, entry)
                if( cached )
                    break
            }
            catch (Throwable t) {
                log.warn1("[${safeTaskName(task)}] Unable to resume cached task -- See log file for details", causedBy: t)
            }

            if( exists ) {
                tries++
                continue
            }

            final lock = lockManager.acquire(hash)
            final workDir = task.getWorkDirFor(hash)
            try {
                if( resumeDir != workDir )
                    exists = workDir.exists()
                if( exists ) {
                    tries++
                    continue
                }
                else if( !workDir.mkdirs() )
                    throw new IOException("Unable to create directory=$workDir -- check file system permissions")
            }
            finally {
                lock.release()
            }

            // submit task for execution
            submitTask( task, hash, workDir )
            break
        }

    }

    /**
     * Check whenever the outputs for the specified task already exist
     *
     * @param task The task instance
     * @param folder The folder where the outputs are stored (eventually)
     * @return {@code true} when all outputs are available, {@code false} otherwise
     */
    final boolean checkCachedOutput(TaskRun task, Path folder, HashCode hash, TaskEntry entry) {

        // check if exists the task exit code file
        def exitCode = null
        def exitFile = folder.resolve(TaskRun.CMD_EXIT)
        if( task.type == ScriptType.SCRIPTLET ) {
            def str
            try {
                str = exitFile.text?.trim()
            }
            catch( IOException e ) {
                log.trace "[${safeTaskName(task)}] Exit file can't be read > $exitFile -- return false -- Cause: ${e.message}"
                return false
            }

            exitCode = str.isInteger() ? str.toInteger() : null
            if( !task.isSuccess(exitCode) ) {
                log.trace "[${safeTaskName(task)}] Exit code is not valid > $str -- return false"
                return false
            }
        }

        /*
         * verify cached context map
         */
        if( !entry ) {
            log.trace "[${safeTaskName(task)}] Missing cache entry -- return false"
            return false
        }

        if( task.hasCacheableValues() && !entry.context ) {
            log.trace "[${safeTaskName(task)}] Missing cache context -- return false"
            return false
        }

        /*
         * verify stdout file
         */
        final stdoutFile = folder.resolve( TaskRun.CMD_OUTFILE )

        if( entry.context != null ) {
            task.context = entry.context
            task.config.context = entry.context
            task.code?.delegate = entry.context
        }

        try {
            // -- set task properties in order to resolve task outputs
            task.workDir = folder
            task.stdout = stdoutFile
            task.config.exitStatus = exitCode
            // -- check if all output resources are available
            collectOutputs(task)

            // set the exit code in to the task object
            task.cached = true
            task.hash = hash
            if( exitCode != null ) {
                task.exitStatus = exitCode
            }

            log.info "[${task.hashLog}] Cached process > ${task.name}"
            // -- notify cached event
            if( entry )
                session.notifyTaskCached(new CachedTaskHandler(task,entry.trace))

            // -- now bind the results
            onComplete.accept(task)
            return true
        }
        catch( MissingFileException | MissingValueException e ) {
            log.trace "[${safeTaskName(task)}] Missed cache > ${e.getMessage()} -- folder: $folder"
            task.exitStatus = Integer.MAX_VALUE
            task.workDir = null
            return false
        }
    }

    /**
     * Execute the specified task shell script
     *
     * @param script The script string to be execute, e.g. a BASH script
     * @return {@code TaskDef}
     */
    final protected void submitTask( TaskRun task, HashCode hash, Path folder ) {
        log.trace "[${safeTaskName(task)}] actual run folder: ${folder}"

        // set name, hash, and working directory
        task.hash = hash
        task.workDir = folder
        task.config.workDir = folder
        task.config.hash = hash.toString()
        task.config.name = task.getName()

        // when no collector is define OR it's a task retry, then submit directly for execution
        final arrayCollector = processor.getArrayCollector()
        if( !arrayCollector || task.config.getAttempt() > 1 )
            executor.submit(task)
        // add the task to the collection of running tasks
        else
            arrayCollector.collect(task)
    }

    protected boolean checkWhenGuard(TaskRun task) {

        try {
            def pass = task.config.getWhenGuard()
            if( pass ) {
                return true
            }

            log.trace "Task ${safeTaskName(task)} is not executed because `when` condition is not verified"
            onComplete.accept(task)
            return false
        }
        catch ( FailedGuardException error ) {
            handleException(error, task)
            return false
        }
    }

    /**
     * Finalize the task execution, checking the exit status
     * and binding output values accordingly
     *
     * @param task The {@code TaskRun} instance to finalize
     */
    @PackageScope
    final finalizeTask( TaskHandler handler) {
        def task = handler.task
        log.trace "finalizing process > ${safeTaskName(task)} -- $task"

        def fault = null
        try {
            // -- verify task exit status
            if( task.error )
                throw new ProcessFailedException("Process `${safeTaskName(task)}` failed", task.error)

            if( task.type == ScriptType.SCRIPTLET ) {
                if( task.exitStatus == Integer.MAX_VALUE )
                    throw new ProcessFailedException("Process `${safeTaskName(task)}` terminated for an unknown reason -- Likely it has been terminated by the external system")

                if ( !task.isSuccess() )
                    throw new ProcessFailedException("Process `${safeTaskName(task)}` terminated with an error exit status (${task.exitStatus})")
            }

            // -- expose task exit status to make accessible as output value
            task.config.exitStatus = task.exitStatus
            // -- if it's OK collect results and finalize
            collectOutputs(task)
        }
        catch ( Throwable error ) {
            fault = resumeOrDie(task, error, handler.getTraceRecord())
            log.trace "Task fault (3): $fault"
        }

        // -- finalize the task
        if( fault != ErrorStrategy.RETRY )
            onComplete.accept(task)

        return fault
    }

    /**
     * Once the task has completed this method is invoked to collected all the task results
     *
     * @param task
     */
    @CompileStatic
    protected void collectOutputs( TaskRun task ) {
        if( config instanceof ProcessConfigV2 )
            collectOutputsV2( task )
        else if( config instanceof ProcessConfigV1 )
            collectOutputsV1( task )
    }

    @CompileStatic
    protected void collectOutputsV2(TaskRun task) {
        final declaredOutputs = configV2().getOutputs()
        final resolver = new TaskOutputResolver(declaredOutputs.getFiles(), task)

        for( final param : declaredOutputs.getParams() ) {
            final value = resolver.resolve(param.getLazyValue())
            task.setOutput(param, value)
        }

        for( final topic : declaredOutputs.getTopics() ) {
            final value = resolver.resolve(topic.getLazyValue())
            topic.getChannel().bind(value)
        }

        task.canBind = true
    }

    @CompileStatic
    protected void collectOutputsV1( TaskRun task ) {
        log.trace "<$name> collecting output: ${task.outputs}"

        final resolver = new TaskOutputResolverV1(task)

        for( OutParam param : task.outputs.keySet() )
            resolver.resolve(param)

        // mark ready for output binding
        task.canBind = true
    }

    /**
     * Handles an error raised during the processor execution
     *
     * @param error The exception raised during the task execution
     * @param task The {@code TaskDef} instance which raised the exception
     * @return {@code true} to terminate the processor execution,
     *         {@code false} ignore the error and continue to process other pending tasks
     */
    final protected boolean handleException( Throwable error, TaskRun task = null ) {
        log.trace "Handling error: $error -- task: $task"
        def fault = resumeOrDie(task, error)
        log.trace "Task fault (2): $fault"

        if (fault instanceof TaskFault) {
            session.fault(fault)
            // when a `TaskFault` is returned a `TERMINATE` is implicit, thus return `true`
            return true
        }

        return fault == TERMINATE || fault == FINISH
    }

    /**
     * @param task The {@code TaskRun} instance that raised an error
     * @param error The error object
     * @return
     *      Either a value of value of {@link ErrorStrategy} representing the error strategy chosen
     *      or an instance of {@TaskFault} representing the cause of the error (that implicitly means
     *      a {@link ErrorStrategy#TERMINATE})
     */
    @PackageScope
    final synchronized resumeOrDie( TaskRun task, Throwable error, TraceRecord traceRecord = null) {
        log.debug "Handling unexpected condition for\n  task: name=${safeTaskName(task)}; work-dir=${task?.workDirStr}\n  error [${error?.class?.name}]: ${error?.getMessage()?:error}"

        ErrorStrategy errorStrategy = TERMINATE
        final List<String> message = []
        try {
            // -- do not recoverable error, just re-throw it
            if( error instanceof Error ) throw error

            // -- retry without increasing the error counts
            if( task && (error.cause instanceof ProcessRetryableException || error.cause instanceof CloudSpotTerminationException) ) {
                if( error.cause instanceof ProcessRetryableException )
                    log.info "[$task.hashLog] NOTE: ${error.message} -- Execution is retried"
                else
                    log.info "[$task.hashLog] NOTE: ${error.message} -- Cause: ${error.cause.message} -- Execution is retried"
                task.failCount+=1
                final taskCopy = task.makeCopy()
                session.getExecService().submit {
                    try {
                        taskCopy.runType = TaskProcessor.RunType.RETRY
                        checkCachedOrLaunchTask( taskCopy, taskCopy.hash, false )
                    }
                    catch( Throwable e ) {
                        log.error("Unable to re-submit task `${taskCopy.name}`", e)
                        session.abort(e)
                    }
                }
                task.failed = true
                task.errorAction = RETRY
                return RETRY
            }

            final submitTimeout = error.cause instanceof ProcessSubmitTimeoutException
            final submitErrMsg = submitTimeout ? error.cause.message : null
            final int submitRetries = submitTimeout ? ++task.submitRetries : 0
            final int taskErrCount = !submitTimeout && task ? ++task.failCount : 0
            final int procErrCount = !submitTimeout ? ++errorCount : errorCount

            // -- when is a task level error and the user has chosen to ignore error,
            //    just report and error message and DO NOT stop the execution
            if( task && error instanceof ProcessException ) {
                // expose current task exit status
                task.config.exitStatus = task.exitStatus
                task.config.errorCount = procErrCount
                task.config.retryCount = taskErrCount
                //Add trace of the previous execution in the task context for next execution
                if ( traceRecord )
                    task.config.previousTrace = traceRecord
                task.config.previousException = error

                errorStrategy = checkErrorStrategy(task, error, taskErrCount, procErrCount, submitRetries)
                if( errorStrategy.soft ) {
                    def msg = "[$task.hashLog] NOTE: ${submitTimeout ? submitErrMsg : error.message}"
                    if( errorStrategy == IGNORE )
                        msg += " -- Error is ignored"
                    else if( errorStrategy == RETRY )
                        msg += " -- Execution is retried (${submitTimeout ? submitRetries : taskErrCount})"
                    log.info msg
                    task.failed = true
                    task.errorAction = errorStrategy
                    return errorStrategy
                }
            }

            // -- mark the task as failed
            if( task ) {
                task.failed = true
                task.errorAction = errorStrategy
            }

            // -- make sure the error is showed only the very first time across all processes
            if( errorShown.getAndSet(true) || session.aborted ) {
                log.trace "Task errorShown=${errorShown.get()}; aborted=${session.aborted}"
                return errorStrategy
            }

            def dumpStackTrace = log.isTraceEnabled()
            def errorFormatter = new TaskErrorFormatter()
            message << "Error executing process > '${safeTaskName(task)}'"
            switch( error ) {
                case ProcessException:
                    errorFormatter.formatTaskError( message, error, task )
                    break

                case ProcessEvalException:
                    errorFormatter.formatCommandError( message, error, task )
                    break

                case FailedGuardException:
                    errorFormatter.formatGuardError( message, error as FailedGuardException, task )
                    break;

                default:
                    message << errorFormatter.formatErrorCause(error)
                    dumpStackTrace = true
            }

            if( dumpStackTrace )
                log.error(message.join('\n'), error)
            else
                log.error(message.join('\n'))
        }
        catch( Throwable e ) {
            // no recoverable error
            log.error("Execution aborted due to an unexpected error", e )
        }

        return new TaskFault(error: error, task: task, report: message.join('\n'))
    }

    protected ErrorStrategy checkErrorStrategy( TaskRun task, ProcessException error, final int taskErrCount, final int procErrCount, final submitRetries ) {

        final action = task.config.getErrorStrategy()

        // retry is not allowed when the script cannot be compiled or similar errors
        if( error instanceof ProcessUnrecoverableException || error.cause instanceof ProcessUnrecoverableException ) {
            return !action.soft ? action : TERMINATE
        }

        // IGNORE strategy -- just continue
        if( action == IGNORE ) {
            return IGNORE
        }

        // RETRY strategy -- check that process do not exceed 'maxError' and the task do not exceed 'maxRetries'
        if( action == RETRY ) {
            final int maxErrors = task.config.getMaxErrors()
            final int maxRetries = task.config.getMaxRetries()

            if( (procErrCount < maxErrors || maxErrors == -1) && taskErrCount <= maxRetries && submitRetries <= maxRetries ) {
                final taskCopy = task.makeCopy()
                session.getExecService().submit({
                    try {
                        taskCopy.config.attempt = taskErrCount+1
                        taskCopy.config.submitAttempt = submitRetries+1
                        taskCopy.runType = TaskProcessor.RunType.RETRY
                        taskCopy.resolve(taskBody)
                        checkCachedOrLaunchTask( taskCopy, taskCopy.hash, false )
                    }
                    catch( Throwable e ) {
                        log.error("Unable to re-submit task `${safeTaskName(taskCopy)}`", e)
                        session.abort(e)
                    }
                } as Runnable)
                return RETRY
            }

            return TERMINATE
        }

        return action
    }
}
