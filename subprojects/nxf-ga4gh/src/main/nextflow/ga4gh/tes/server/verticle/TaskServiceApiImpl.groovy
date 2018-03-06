package nextflow.ga4gh.tes.server.verticle

import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import io.vertx.core.AsyncResult
import io.vertx.core.Handler
import nextflow.Const
import nextflow.Session
import nextflow.executor.Executor
import nextflow.executor.LocalExecutor
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.ga4gh.tes.server.model.TesCancelTaskResponse
import nextflow.ga4gh.tes.server.model.TesCreateTaskResponse
import nextflow.ga4gh.tes.server.model.TesExecutor
import nextflow.ga4gh.tes.server.model.TesInput
import nextflow.ga4gh.tes.server.model.TesListTasksResponse
import nextflow.ga4gh.tes.server.model.TesOutput
import nextflow.ga4gh.tes.server.model.TesResources
import nextflow.ga4gh.tes.server.model.TesServiceInfo
import nextflow.ga4gh.tes.server.model.TesState
import nextflow.ga4gh.tes.server.model.TesTask
import nextflow.ga4gh.tes.server.model.TesTaskLog
import nextflow.processor.TaskContext
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.CacheHelper
/**
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
@CompileStatic
class TaskServiceApiImpl implements TaskServiceApi {

    final protected AtomicInteger indexCount = new AtomicInteger()

    /**
     * Lock the work dir creation process
     */
    private static Lock lockWorkDirCreation = new ReentrantLock()

    private Session session

    private Executor executor

    TaskServiceApiImpl() {
        init()
    }

    TaskServiceApiImpl(Session session) {
        init(session)
    }

    protected init(Session session=null) {
        this.session = session ?: new Session()
        this.session.runName = 'tes-worker'
        this.session.init(null)
        this.session.start()

        this.executor = new LocalExecutor(session: session, name: 'local')
        this.executor.init()

    }

    static class Response<T> implements AsyncResult<T> {

        Throwable error

        Closure<T> action

        Response(Closure<T> action) {
            this.action = action
        }

        @Override
        T result() {
            try {
                return action.call()
            }
            catch( Throwable e ) {
                error = e
                e.printStackTrace()
                return null
            }
        }

        @Override
        Throwable cause() {
            return error
        }

        @Override
        boolean succeeded() {
            return error == null
        }

        @Override
        boolean failed() {
            return error != null
        }
    }

    protected <V> Response<V> response( Closure<V> action )  {
        new Response<V>(action)
    }

    @Override
    void cancelTask(String id, Handler<AsyncResult<TesCancelTaskResponse>> handler) {

    }

    protected Path createPath(HashCode hash) {

        int tries = 1
        while( true ) {
            hash = CacheHelper.defaultHasher().newHasher().putBytes(hash.asBytes()).putInt(tries).hash()

            boolean exists=false
            final folder = FileHelper.getWorkFolder(session.workDir, hash)
            lockWorkDirCreation.lock()
            try {
                exists = folder.exists()
                if( !exists && !folder.mkdirs() )
                    throw new IOException("Unable to create folder=$folder -- check file system permission")
            }
            finally {
                lockWorkDirCreation.unlock()
            }

            if( exists ) {
                tries++
                continue
            }

            return folder
        }

    }

    protected HashCode getTaskHash(TaskRun task) {
        List keys = [ session.uniqueId, task.name, task.source, task.container ] as List<Object>

        // add all the input name-value pairs to the key generator
        task.inputs.each {
            keys.add( it.key.name )
            keys.add( it.value )
        }

        CacheHelper.hasher(keys).hash()
    }

    protected TesCreateTaskResponse submitTask(TesTask req) {

        if( req.executors?.size()==0 )
            throw new TaskServiceApiException(400, "Malformed TES req -- Missing executor definition")

        def proc = new TesTaskProcessorAdaptor(session,executor)
        def config = proc.config.createTaskConfig()
        def task = new TaskRun(
                id: TaskId.next(),
                index: indexCount.incrementAndGet(),
                processor: proc,
                type: ScriptType.SCRIPTLET,
                config: config,
                context: new TaskContext(proc)
        )

        int count=0
        for( TesInput input : req.inputs ) {
            final url = FileHelper.asPath( input.getUrl() )
            final path = Paths.get( input.getPath() )
            task.setInput( new TesFileInputParamAdaptor("in-file-${++count}"), new FileHolder(url).withName(path) )
        }

        TesExecutor exec = req.executors.get(0)
        if( !exec.image )
            throw new TaskServiceApiException(400, "Malformed TES request -- Missing container image")
        if( !exec.command?.size() )
            throw new TaskServiceApiException(400, "Malformed TES request -- Missing task command")

        config.container = exec.image
        task.name = req.name ?: "tes-task-${indexCount.get()}"
        task.script = exec.command.get(0)
        task.source = exec.command.get(0)
        task.hash = getTaskHash(task)
        task.workDir = createPath(task.hash)

        // -- submit the task execution
        session.dispatcher.submit( task, false )

        // -- create the response object
        return new TesCreateTaskResponse(task.hash.toString())
    }

    @Override
    void createTask(TesTask request, Handler<AsyncResult<TesCreateTaskResponse>> handler) {


        handler.handle( response { submitTask(request) })
    }

    @Override
    void getServiceInfo(Handler<AsyncResult<TesServiceInfo>> handler) {

        handler.handle(response {
            def info = new TesServiceInfo()
            info.name = "Nextflow GA4GH TES back-end v${Const.APP_VER}"
            return info
        })
    }

    @Override
    void getTask(String id, String view, Handler<AsyncResult<TesTask>> handler) {

    }

    @Override
    void listTasks(String namePrefix, Long pageSize, String pageToken, String view, Handler<AsyncResult<TesListTasksResponse>> handler) {

        AsyncResult<TesListTasksResponse> resp = response {
            TesTask t1 = new TesTask(
                    id: "t1",
                    state: TesState.CANCELED,
                    name: "test 100",
                    description: "test 1",
                    inputs: [new TesInput(name: "input 1")],
                    outputs: [new  TesOutput(name: "output 1")],
                    resources: new TesResources(cpuCores: 1),
                    executors: [new TesExecutor(image: "busybox")],
                    volumes: ['/my/volume'],
                    tags: [t1: 'my tag 1'],
                    logs: [new TesTaskLog()],
                    creationTime: "now"
            );

            TesTask t2 = new TesTask(
                    id: "t2",
                    state: TesState.COMPLETE,
                    name: "test 2",
                    description: "test 2",
                    inputs: [new TesInput(name: "input 2")],
                    outputs: [new  TesOutput(name: "output 2")],
                    resources: new TesResources(ramGb: 2),
                    executors: [new TesExecutor(image: "debian")],
                    volumes: ['/my/volume2'],
                    tags: [t1: 'my tag 2'],
                    logs: [new TesTaskLog()],
                    creationTime: "yesterday"
            );

            new TesListTasksResponse([t1, t2], null)
        }
        handler.handle(resp)
    }
}
