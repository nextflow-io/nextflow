package nextflow.ga4gh.tes.server.verticle;

import nextflow.ga4gh.tes.server.MainApiException;
import nextflow.ga4gh.tes.server.model.TesCancelTaskResponse;
import nextflow.ga4gh.tes.server.model.TesCreateTaskResponse;
import nextflow.ga4gh.tes.server.model.TesListTasksResponse;
import nextflow.ga4gh.tes.server.model.TesServiceInfo;
import nextflow.ga4gh.tes.server.model.TesTask;

import io.vertx.core.AsyncResult;
import io.vertx.core.Handler;

import java.util.List;
import java.util.Map;

public interface TaskServiceApi  {
    //CancelTask
    void cancelTask(String id, Handler<AsyncResult<TesCancelTaskResponse>> handler);
    
    //CreateTask
    void createTask(TesTask body, Handler<AsyncResult<TesCreateTaskResponse>> handler);
    
    //GetServiceInfo
    void getServiceInfo(Handler<AsyncResult<TesServiceInfo>> handler);
    
    //GetTask
    void getTask(String id, String view, Handler<AsyncResult<TesTask>> handler);
    
    //ListTasks
    void listTasks(String namePrefix, Long pageSize, String pageToken, String view, Handler<AsyncResult<TesListTasksResponse>> handler);
    
}
