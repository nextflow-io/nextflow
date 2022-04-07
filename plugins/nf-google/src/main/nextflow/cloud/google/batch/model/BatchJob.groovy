/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.model


import groovy.transform.CompileStatic
import nextflow.cloud.google.batch.json.JsonHelper
/**
 * Model a Google Batch task submit request
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BatchJob {

    List<TaskGroup> taskGroups = []
    LogsPolicy logsPolicy = LogsPolicy.CLOUD_LOGGING

    BatchJob addTaskGroup(TaskGroup it) {
        taskGroups.add(it)
        return this
    }

    String toJson(boolean pretty=false) {
        JsonHelper.toJson(this,pretty)
    }

    static BatchJob create(Map req ) {
        assert req.imageUri
        assert req.command
        assert req.command instanceof List

        final container = new TaskContainer()
                .withImageUri(req.imageUri as String)
                .withEntrypoint(req.entrypoint as String ?: '/bin/bash')
                .withCommands(req.command as List)

        final taskSpec = new TaskSpec().addRunnable( new TaskRunnable().withContainer(container) )
        final taskGroup = new TaskGroup().withTaskSpec(taskSpec)

        if( req.volumes instanceof Map ) {
            final volumes = (Map)req.volumes
            final contVols = []
            final opts = ["-o rw,allow_other", "-implicit-dirs"]
            for( String remote : volumes.keySet() ) {
                final target = volumes.get(remote)
                contVols << "$target:$target:rw"
                taskSpec.withVolumes( [new TaskVolume(gcs: [remotePath: remote], mountPath: "/mnt/$remote".toString(), mountOptions: opts)] )
            }

            container.withVolumes( contVols )
        }

        return new BatchJob().addTaskGroup(taskGroup)
    }

}
