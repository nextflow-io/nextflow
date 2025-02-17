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

import groovy.util.logging.Slf4j
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowRun
import nextflow.file.FileHelper
import nextflow.script.ScriptMeta
import nextflow.util.PathNormalizer

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.data.cid.model.DataType
import nextflow.data.cid.model.Output
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
@Slf4j
@CompileStatic
class CidObserver implements TraceObserver {
    public static final String METADATA_FILE = '.data.json'
    public static final String CID_PROT = 'cid://'
    private CidStore store
    private Session session

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.store = session.cidStore
    }

    void onFlowBegin() {
        storeWorkflowRun()
    }

    protected void storeWorkflowRun() {
        final normalizer = new PathNormalizer(session.workflowMetadata)
        final mainScript = normalizer.normalizePath(session.workflowMetadata.scriptFile.normalize())
        final workflow = new Workflow(
            DataType.Workflow,
            mainScript,
            ScriptMeta.allScriptNames().values().collect {normalizer.normalizePath(it.normalize())},
            session.workflowMetadata.repository,
            session.workflowMetadata.commitId
        )
        final value = new WorkflowRun(
            DataType.WorkflowRun,
            workflow,
            session.uniqueId.toString(),
            session.runName,
            getNormalizedParams(session.params, normalizer)
        )
        final content = JsonOutput.prettyPrint(JsonOutput.toJson(value))
        store.save("${session.executionHash}/$METADATA_FILE", content)
    }

    private static Map getNormalizedParams(Map<String, Object> params, PathNormalizer normalizer){
        final normalizedParams = new HashMap<String,Object>()
        params.each{String key, Object value ->
            log.debug("Managing parameter $key , class ${value.class}")
            if (value instanceof Path)
                normalizedParams.put(key,normalizer.normalizePath(value as Path))
            else if (value instanceof String || value instanceof GString)
                normalizedParams.put(key,normalizer.normalizePath(value.toString()))
            else
                normalizedParams.put(key, value)
        }
        return normalizedParams
    }


    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        storeTaskInfo(handler.task)
    }

    protected void storeTaskInfo(TaskRun task) {
        final pathNormalizer = new PathNormalizer(session.workflowMetadata)
        // store the task run entry
        storeTaskRun(task, pathNormalizer)
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

    protected void storeTaskRun(TaskRun task, PathNormalizer normalizer) {
        final value = new nextflow.data.cid.model.TaskRun(
            DataType.TaskRun,
            session.uniqueId.toString(),
            task.getName(),
            session.stubRun ? task.stubSource: task.source,
            task.inputFilesMap ? convertToReferences(task.inputFilesMap, normalizer): null,
            task.isContainerEnabled() ? task.getContainerFingerprint(): null,
            normalizer.normalizePath(task.getCondaEnv()),
            normalizer.normalizePath(task.getSpackEnv()),
            task.config?.getArchitecture()?.toString(),
            task.processor.getTaskGlobalVars(task),
            task.processor.getTaskBinEntries(task.source).collect { Path p -> normalizer.normalizePath(p.normalize()) }
            )
        // store in the underlying persistence
        final key = "${task.hash}/$METADATA_FILE"
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    protected void storeTaskOutput(TaskRun task, Path path) {
        final attrs = readAttributes(path)
        final rel = task.workDir.relativize(path).toString()
        final cid = "${task.hash}/${rel}"
        final key = "${cid}/$METADATA_FILE"
        final hash = CacheHelper.hasher(path).hash().toString()
        final value = new Output(
            DataType.TaskOutput,
            path.toString(),
            hash,
            "$CID_PROT$task.hash",
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
        final key = "$session.executionHash/${rel}/$METADATA_FILE"
        final sourceReference = getSourceReference(source)
        final attrs = readAttributes(destination)
        final value = new Output(
            DataType.WorkflowOutput,
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
            return "$CID_PROT$hash/$target"
        }
        return null
    }

    @Override
    void onFilePublish(Path destination){
        final hash = CacheHelper.hasher(destination).hash().toString()
        final rel = session.outputDir.relativize(destination).toString()
        final key = "$session.executionHash/${rel}/$METADATA_FILE"
        final attrs = readAttributes(destination)
        final value = new Output(
            DataType.WorkflowOutput,
            destination.toString(),
            hash,
            session.executionHash,
            attrs.size(),
            attrs.creationTime().toMillis(),
            attrs.lastModifiedTime().toMillis() )
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
    }

    protected List<String> convertToReferences(Map<String, Path> inputs, PathNormalizer normalizer) {
        List<String> references = new LinkedList<String>()
        inputs.each { name, path ->
            final ref = getSourceReference(path)
            references.add(ref ? ref : normalizer.normalizePath(path))}
        return references
    }
}
