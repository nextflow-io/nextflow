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
 *
 */

package nextflow.data.cid

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.data.cid.model.TaskOutput
import nextflow.data.config.DataConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.params.FileOutParam
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CidObserver implements TraceObserver {

    private CidStore store

    @Override
    void onFlowCreate(Session session) {
        store = new DefaultCidStore()
        store.open(DataConfig.create(session))
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        storeTaskInfo(handler.task)
    }

    void storeTaskInfo(TaskRun task) {
        // store the task run entry
        storeTaskRun(task)
        // store all task outputs files
        final outputs = task.getOutputsByType(FileOutParam)
        for( Map.Entry entry : outputs ) {
            final value = entry.value
            if( value instanceof Path ) {
                storeTaskOutput(task, (Path)value)
            }
            else if( value instanceof Collection<Path> ) {
                for( Path it : value )
                    storeTaskOutput(task, (Path)it)
            }
        }
    }

    protected void storeTaskRun(TaskRun task) {
        final value = new nextflow.data.cid.model.TaskRun(
            task.id.value,
            task.getName(),
            task.hash.toString() )
        // store in the underlying persistence
        final key = "${value.hash}/.task"
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    protected void storeTaskOutput(TaskRun task, Path path) {
        final attrs = readAttributes(path)
        final rel = task.workDir.relativize(path).toString()
        final key = "${task.hash}/${rel}"
        final value = new TaskOutput(
            "cid://$key",
            attrs.size(),
            attrs.creationTime().toMillis(),
            attrs.lastModifiedTime().toMillis() )
        // store in the underlying persistence
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    protected BasicFileAttributes readAttributes(Path path) {
        Files.readAttributes(path, BasicFileAttributes)
    }
}
