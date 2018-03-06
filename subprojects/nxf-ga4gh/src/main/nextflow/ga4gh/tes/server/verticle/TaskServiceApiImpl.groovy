package nextflow.ga4gh.tes.server.verticle;

import io.vertx.core.AsyncResult;
import io.vertx.core.Handler
import nextflow.ga4gh.tes.server.model.TesInput
import nextflow.ga4gh.tes.server.model.TesOutput
import nextflow.ga4gh.tes.server.model.TesResources;
import nextflow.ga4gh.tes.server.model.TesCancelTaskResponse;
import nextflow.ga4gh.tes.server.model.TesCreateTaskResponse
import nextflow.ga4gh.tes.server.model.TesExecutor;
import nextflow.ga4gh.tes.server.model.TesListTasksResponse;
import nextflow.ga4gh.tes.server.model.TesServiceInfo
import nextflow.ga4gh.tes.server.model.TesState;
import nextflow.ga4gh.tes.server.model.TesTask
import nextflow.ga4gh.tes.server.model.TesTaskLog;

/**
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
class TaskServiceApiImpl implements TaskServiceApi {
    @Override
    void cancelTask(String id, Handler<AsyncResult<TesCancelTaskResponse>> handler) {

    }

    @Override
    void createTask(TesTask body, Handler<AsyncResult<TesCreateTaskResponse>> handler) {

    }

    @Override
    void getServiceInfo(Handler<AsyncResult<TesServiceInfo>> handler) {

    }

    @Override
    void getTask(String id, String view, Handler<AsyncResult<TesTask>> handler) {

    }

    @Override
    void listTasks(String namePrefix, Long pageSize, String pageToken, String view, Handler<AsyncResult<TesListTasksResponse>> handler) {
        TesTask t1 = new TesTask(
                id: "t1",
                state: TesState.CANCELED,
                name: "test 1",
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
        TesListTasksResponse res = new TesListTasksResponse([t1, t2], null)
        handler.handle(new AsyncResult<TesListTasksResponse>() {
            @Override
            TesListTasksResponse result() {
                return res
            }

            @Override
            Throwable cause() {
                return null
            }

            @Override
            boolean succeeded() {
                return true
            }

            @Override
            boolean failed() {
                return false
            }
        })
    }
}
