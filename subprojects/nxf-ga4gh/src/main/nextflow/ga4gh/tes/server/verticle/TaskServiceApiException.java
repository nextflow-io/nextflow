package nextflow.ga4gh.tes.server.verticle;

import nextflow.ga4gh.tes.server.MainApiException;
import nextflow.ga4gh.tes.server.model.TesCancelTaskResponse;
import nextflow.ga4gh.tes.server.model.TesCreateTaskResponse;
import nextflow.ga4gh.tes.server.model.TesListTasksResponse;
import nextflow.ga4gh.tes.server.model.TesServiceInfo;
import nextflow.ga4gh.tes.server.model.TesTask;

public final class TaskServiceApiException extends MainApiException {
    public TaskServiceApiException(int statusCode, String statusMessage) {
        super(statusCode, statusMessage);
    }
}