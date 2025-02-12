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

import com.google.common.hash.HashCode
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowRun
import nextflow.file.FileHelper

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.data.cid.model.DataType
import nextflow.data.cid.model.Output
import nextflow.data.config.DataConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.params.FileOutParam
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CidObserver implements TraceObserver {

    private CidStore store
    private Session session

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        store = new DefaultCidStore()
        store.open(DataConfig.create(session))
    }

    void onFlowBegin() {
        storeWorkflowRun()
    }

    protected void storeWorkflowRun() {
        final workflow = new Workflow(
            DataType.Workflow,
            session.workflowMetadata.scriptFile.toString(),
            session.workflowMetadata.scriptId.toString(),
            session.workflowMetadata.repository,
            session.workflowMetadata.commitId
        )
        final value = new WorkflowRun(
            DataType.WorkflowRun,
            workflow,
            session.uniqueId.toString(),
            session.runName,
            session.params
        )
        final content = JsonOutput.prettyPrint(JsonOutput.toJson(value))
        store.save("${session.executionHash}/.data.json", content)
    }
    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        storeTaskInfo(handler.task)
    }

    protected void storeTaskInfo(TaskRun task) {
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
            DataType.Task,
            task.id.value,
            task.getName(),
            task.hash.toString(),
            task.inputFilesMap ? convertToReferences(task.inputFilesMap): null
            )
        // store in the underlying persistence
        final key = "${value.hash}/.data.json"
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    protected void storeTaskOutput(TaskRun task, Path path) {
        final attrs = readAttributes(path)
        final rel = task.workDir.relativize(path).toString()
        final cid = "${task.hash}/${rel}"
        final key = "${cid}/.data.json"
        final hash = CacheHelper.hasher(path).hash().toString()
        final value = new Output(
            DataType.Output,
            path.toString(),
            hash,
            "cid://$task.hash",
            attrs.size(),
            attrs.creationTime().toMillis(),
            attrs.lastModifiedTime().toMillis() )
        // store in the underlying persistence
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    protected BasicFileAttributes readAttributes(Path path) {
        Files.readAttributes(path, BasicFileAttributes)
    }

    @Override
    void onFilePublish(Path destination, Path source){
        final hash = CacheHelper.hasher(destination).hash().toString()
        final rel = session.outputDir.relativize(destination).toString()
        final key = "${rel}/.data.json"
        final sourceReference = getSourceReference(source)
        final attrs = readAttributes(destination)
        final value = new Output(
            DataType.Output,
            destination.toString(),
            hash,
            sourceReference,
            attrs.size(),
            attrs.creationTime().toMillis(),
            attrs.lastModifiedTime().toMillis() )
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    String getSourceReference(Path source){
        final hash = FileHelper.getTaskHashFromPath(source, session.workDir)
        if (hash) {
            final target = FileHelper.getWorkFolder(session.workDir, hash).relativize(source).toString()
            return "cid://$hash/$target"
        }
        return null
    }

    @Override
    void onFilePublish(Path destination){
        final hash = CacheHelper.hasher(destination).hash().toString()
        final rel = session.outputDir.relativize(destination).toString()
        final attrs = readAttributes(destination)
        final value = new Output(
            DataType.Output,
            destination.toString(),
            hash,
            session.executionHash,
            attrs.size(),
            attrs.creationTime().toMillis(),
            attrs.lastModifiedTime().toMillis() )
        store.save(rel, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    protected Map convertToReferences(Map<String, Path> inputs) {
        Map<String, String> references = new HashMap<String, String>()
        inputs.each { name, path ->
            final ref = getSourceReference(path)
            references.put(name, ref ? ref : path.toString())}
        return references
    }
}
