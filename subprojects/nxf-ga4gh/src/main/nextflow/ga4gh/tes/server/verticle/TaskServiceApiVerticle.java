package nextflow.ga4gh.tes.server.verticle;

import io.vertx.core.AbstractVerticle;
import io.vertx.core.eventbus.Message;
import io.vertx.core.json.Json;
import io.vertx.core.json.JsonObject;
import io.vertx.core.logging.Logger;
import io.vertx.core.logging.LoggerFactory;
import nextflow.ga4gh.tes.server.MainApiException;
import nextflow.ga4gh.tes.server.model.TesTask;

public class TaskServiceApiVerticle extends AbstractVerticle {
    final static Logger LOGGER = LoggerFactory.getLogger(TaskServiceApiVerticle.class); 
    
    final static String CANCELTASK_SERVICE_ID = "CancelTask";
    final static String CREATETASK_SERVICE_ID = "CreateTask";
    final static String GETSERVICEINFO_SERVICE_ID = "GetServiceInfo";
    final static String GETTASK_SERVICE_ID = "GetTask";
    final static String LISTTASKS_SERVICE_ID = "ListTasks";
    
    TaskServiceApi service = new TaskServiceApiImpl();

    @Override
    public void start() throws Exception {
        
        //Consumer for CancelTask
        vertx.eventBus().<JsonObject> consumer(CANCELTASK_SERVICE_ID).handler(message -> {
            try {
                String id = message.body().getString("id");
                service.cancelTask(id, result -> {
                    if (result.succeeded())
                        message.reply(new JsonObject(Json.encode(result.result())).encodePrettily());
                    else {
                        Throwable cause = result.cause();
                        manageError(message, cause, "CancelTask");
                    }
                });
            } catch (Exception e) {
                logUnexpectedError("CancelTask", e);
                message.fail(MainApiException.INTERNAL_SERVER_ERROR.getStatusCode(), MainApiException.INTERNAL_SERVER_ERROR.getStatusMessage());
            }
        });
        
        //Consumer for CreateTask
        vertx.eventBus().<JsonObject> consumer(CREATETASK_SERVICE_ID).handler(message -> {
            try {
                TesTask body = Json.mapper.readValue(message.body().getJsonObject("body").encode(), TesTask.class);
                service.createTask(body, result -> {
                    if (result.succeeded())
                        message.reply(new JsonObject(Json.encode(result.result())).encodePrettily());
                    else {
                        Throwable cause = result.cause();
                        manageError(message, cause, "CreateTask");
                    }
                });
            } catch (Exception e) {
                logUnexpectedError("CreateTask", e);
                message.fail(MainApiException.INTERNAL_SERVER_ERROR.getStatusCode(), MainApiException.INTERNAL_SERVER_ERROR.getStatusMessage());
            }
        });
        
        //Consumer for GetServiceInfo
        vertx.eventBus().<JsonObject> consumer(GETSERVICEINFO_SERVICE_ID).handler(message -> {
            try {
                service.getServiceInfo(result -> {
                    if (result.succeeded())
                        message.reply(new JsonObject(Json.encode(result.result())).encodePrettily());
                    else {
                        Throwable cause = result.cause();
                        manageError(message, cause, "GetServiceInfo");
                    }
                });
            } catch (Exception e) {
                logUnexpectedError("GetServiceInfo", e);
                message.fail(MainApiException.INTERNAL_SERVER_ERROR.getStatusCode(), MainApiException.INTERNAL_SERVER_ERROR.getStatusMessage());
            }
        });
        
        //Consumer for GetTask
        vertx.eventBus().<JsonObject> consumer(GETTASK_SERVICE_ID).handler(message -> {
            try {
                String id = message.body().getString("id");
                String view = message.body().getString("view");
                service.getTask(id, view, result -> {
                    if (result.succeeded())
                        message.reply(new JsonObject(Json.encode(result.result())).encodePrettily());
                    else {
                        Throwable cause = result.cause();
                        manageError(message, cause, "GetTask");
                    }
                });
            } catch (Exception e) {
                logUnexpectedError("GetTask", e);
                message.fail(MainApiException.INTERNAL_SERVER_ERROR.getStatusCode(), MainApiException.INTERNAL_SERVER_ERROR.getStatusMessage());
            }
        });
        
        //Consumer for ListTasks
        vertx.eventBus().<JsonObject> consumer(LISTTASKS_SERVICE_ID).handler(message -> {
            try {
                String namePrefix = message.body().getString("name_prefix");
                Long pageSize = Json.mapper.readValue(message.body().getString("page_size"), Long.class);
                String pageToken = message.body().getString("page_token");
                String view = message.body().getString("view");
                service.listTasks(namePrefix, pageSize, pageToken, view, result -> {
                    if (result.succeeded())
                        message.reply(new JsonObject(Json.encode(result.result())).encodePrettily());
                    else {
                        Throwable cause = result.cause();
                        manageError(message, cause, "ListTasks");
                    }
                });
            } catch (Exception e) {
                logUnexpectedError("ListTasks", e);
                message.fail(MainApiException.INTERNAL_SERVER_ERROR.getStatusCode(), MainApiException.INTERNAL_SERVER_ERROR.getStatusMessage());
            }
        });
        
    }
    
    private void manageError(Message<JsonObject> message, Throwable cause, String serviceName) {
        int code = MainApiException.INTERNAL_SERVER_ERROR.getStatusCode();
        String statusMessage = MainApiException.INTERNAL_SERVER_ERROR.getStatusMessage();
        if (cause instanceof MainApiException) {
            code = ((MainApiException)cause).getStatusCode();
            statusMessage = ((MainApiException)cause).getStatusMessage();
        } else {
            logUnexpectedError(serviceName, cause); 
        }
            
        message.fail(code, statusMessage);
    }
    
    private void logUnexpectedError(String serviceName, Throwable cause) {
        LOGGER.error("Unexpected error in "+ serviceName, cause);
    }
}
