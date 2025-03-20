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
import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.WorkflowResults
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowRun
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.script.ScriptMeta
import nextflow.script.params.DefaultInParam
import nextflow.script.params.FileInParam
import nextflow.script.params.InParam
import nextflow.util.PathNormalizer
import nextflow.util.TestOnly

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

import static nextflow.data.cid.fs.CidPath.CID_PROT

/**
 * Observer to write the generated workflow metadata in a CID store.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CidObserver implements TraceObserver {

    private String executionHash
    private CidStore store
    private Session session
    private WorkflowResults workflowResults
    private Map<String,String> outputsStoreDirCid = new HashMap<String,String>(10)

    CidObserver(Session session, CidStore store){
        this.session = session
        this.store = store
    }

    @Override
    void onFlowCreate(Session session) {
        this.store.getHistoryLog().write(session.runName, session.uniqueId, '-', '-')
    }

    @TestOnly
    String getExecutionHash(){ executionHash }

    @Override
    void onFlowBegin() {
        this.executionHash = storeWorkflowRun()
        workflowResults = new WorkflowResults(
            DataType.WorkflowResults,
            System.currentTimeMillis(),
            new HashMap<String, Object>(),
            new ArrayList<String>(),
        )
        this.store.getHistoryLog().updateRunCid(session.uniqueId, "${CID_PROT}${this.executionHash}")
    }

    @Override
    void onFlowComplete(){
        if (this.workflowResults){
            final content = JsonOutput.prettyPrint(JsonOutput.toJson(workflowResults))
            workflowResults.creationTime = System.currentTimeMillis()
            this.store.save("${this.executionHash}/outputs", content)
            this.store.getHistoryLog().updateResultsCid(session.uniqueId, "${CID_PROT}${executionHash}/outputs")
        }
    }

    protected String storeWorkflowRun() {
        final normalizer = new PathNormalizer(session.workflowMetadata)
        final mainScript = new DataPath(
            normalizer.normalizePath(session.workflowMetadata.scriptFile.normalize()),
            new Checksum(session.workflowMetadata.scriptId, "nextflow", CacheHelper.HashMode.DEFAULT().toString().toLowerCase())
        )
        List<DataPath> otherScripts  = new LinkedList<>()
        for (Path p: ScriptMeta.allScriptNames().values()) {
            if (p && p != session.workflowMetadata.scriptFile) {
                otherScripts.add(
                    new DataPath(
                        normalizer.normalizePath(p.normalize()),
                        new Checksum(
                            CacheHelper.hasher(p.text).hash().toString(),
                            "nextflow",
                            CacheHelper.HashMode.DEFAULT().toString().toLowerCase()
                        )
                    )
                )
            }
        }
        final workflow = new Workflow(
            DataType.Workflow,
            mainScript,
            otherScripts,
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
        final executionHash = CacheHelper.hasher(content).hash().toString()
        store.save(executionHash, content)
        return executionHash
    }

    private static List<Parameter> getNormalizedParams(Map<String, Object> params, PathNormalizer normalizer){
        final normalizedParams = new LinkedList<Parameter>()
        params.each{String key, Object value ->
            if( value instanceof Path )
                normalizedParams.add( new Parameter( Path.class.simpleName, key, normalizer.normalizePath( value as Path ) ) )
            else if ( value instanceof CharSequence )
                normalizedParams.add( new Parameter( String.class.simpleName, key, normalizer.normalizePath( value.toString() ) ) )
            else
                normalizedParams.add( new Parameter( value.class.simpleName, key, value) )
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
        outputs.forEach { FileOutParam key, Object value -> manageFileOutParams(value, task)}

    }

    private void manageFileOutParams( Object value, TaskRun task) {
        if (value instanceof Path) {
            storeTaskOutput(task, (Path) value)
        } else if (value instanceof Collection<Path>) {
            for (Path it : value) {
                storeTaskOutput(task, (Path) it)
            }
        }
    }

    protected String storeTaskRun(TaskRun task, PathNormalizer normalizer) {
        final codeChecksum = new Checksum(CacheHelper.hasher(session.stubRun ? task.stubSource: task.source).hash().toString(),
            "nextflow", CacheHelper.HashMode.DEFAULT().toString().toLowerCase())
        final value = new nextflow.data.cid.model.TaskRun(
            DataType.TaskRun,
            session.uniqueId.toString(),
            task.getName(),
            codeChecksum,
            task.inputs ? manageInputs(task.inputs, normalizer): null,
            task.isContainerEnabled() ? task.getContainerFingerprint(): null,
            normalizer.normalizePath(task.getCondaEnv()),
            normalizer.normalizePath(task.getSpackEnv()),
            task.config?.getArchitecture()?.toString(),
            task.processor.getTaskGlobalVars(task),
            task.processor.getTaskBinEntries(task.source).collect { Path p -> new DataPath(
                normalizer.normalizePath(p.normalize()),
                new Checksum(CacheHelper.hasher(p).hash().toString(), "nextflow",
                    CacheHelper.HashMode.DEFAULT().toString().toLowerCase()) )
            }
        )

        // store in the underlying persistence
        final key = task.hash.toString()
        store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
        return key
    }

    protected void storeTaskOutput(TaskRun task, Path path) {
        try {
            final attrs = readAttributes(path)
            final rel = getTaskRelative(task, path)
            final cid = "${task.hash}/${rel}"
            final key = cid.toString()
            final checksum = new Checksum( CacheHelper.hasher(path).hash().toString(),
                "nextflow", CacheHelper.HashMode.DEFAULT().toString().toLowerCase() )
            final value = new Output(
                DataType.TaskOutput,
                path.toUriString(),
                checksum,
                "$CID_PROT$task.hash",
                attrs.size(),
                attrs.creationTime().toMillis(),
                attrs.lastModifiedTime().toMillis())
            store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
        } catch (Throwable e) {
            log.warn("Exception storing CID output $path for task ${task.name}. ${e.getLocalizedMessage()}")
        }
    }

    protected String getTaskRelative(TaskRun task, Path path){
        if (path.isAbsolute()) {
            final rel = getTaskRelative0(task, path)
            if (rel) return rel
            throw new Exception("Cannot asses the relative path for output $path of ${task.name}")
        } else {
            //Check if contains workdir or storeDir
            final rel = getTaskRelative0(task, path.toAbsolutePath())
            if (rel) return rel
            if (path.normalize().getName(0).toString() == "..")
                throw new Exception("Cannot asses the relative path for output $path of ${task.name}" )
            return path.normalize().toString()
        }

    }

    private String getTaskRelative0(TaskRun task, Path path){
        final workDirAbsolute = task.workDir.toAbsolutePath()
            if (path.startsWith(workDirAbsolute)) {
                return workDirAbsolute.relativize(path).toString()
            }
            //If task output is not in the workDir check if output is stored in the task's storeDir
            final storeDir = task.getConfig().getStoreDir().toAbsolutePath()
            if( storeDir && path.startsWith(storeDir)) {
                final rel = storeDir.relativize(path)
                //If output stored in storeDir, keep the path in case it is used as workflow output
                this.outputsStoreDirCid.put(path.toString(), "$CID_PROT${task.hash}/$rel".toString())
                return rel
            }
    }

    protected BasicFileAttributes readAttributes(Path path) {
        Files.readAttributes(path, BasicFileAttributes)
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        storePublishedFile(destination, source)
    }

    protected void storePublishedFile(Path destination, Path source = null, Map annotations = null){
        try {
            final checksum = new Checksum(
                CacheHelper.hasher(destination).hash().toString(),
                "nextflow",
                CacheHelper.HashMode.DEFAULT().toString().toLowerCase()
            )
            final rel = getWorkflowRelative(destination)
            final key = "$executionHash/outputs/${rel}"

            final sourceReference = source ? getSourceReference(source) : "${CID_PROT}${executionHash}".toString()
            final attrs = readAttributes(destination)
            final value = new Output(
                DataType.WorkflowOutput,
                destination.toUriString(),
                checksum,
                sourceReference,
                attrs.size(),
                attrs.creationTime().toMillis(),
                attrs.lastModifiedTime().toMillis(),
                annotations)
            store.save(key, JsonOutput.prettyPrint(JsonOutput.toJson(value)))
            workflowResults.publishedFiles.add("${CID_PROT}${key}".toString())
        } catch (Throwable e) {
            log.warn("Exception storing published file $destination for workflow ${executionHash}.", e)
        }
    }



    String getSourceReference(Path source){
        final hash = FileHelper.getTaskHashFromPath(source, session.workDir)
        if (hash) {
            final target = FileHelper.getWorkFolder(session.workDir, hash).relativize(source).toString()
            return "$CID_PROT$hash/$target"
        } else {
            final storeDirReference = outputsStoreDirCid.get(source.toString())
            if (storeDirReference)
                return "$CID_PROT$storeDirReference"
        }
        return null
    }

    @Override
    void onFilePublish(Path destination){
        storePublishedFile (destination)
    }

    @Override
    void onWorkflowPublish(String name, Object value){
        workflowResults.outputs.put(name,convertPathsToCidReferences(value))
    }

    private Object convertPathsToCidReferences(Object value){
        if( value instanceof Path ) {
            final rel = getWorkflowRelative(value)
            return rel ? "${CID_PROT}${executionHash}/outputs/${rel}".toString() : value
        }

        if( value instanceof Collection ) {
            return value.collect { el -> convertPathsToCidReferences(el) }
        }

        if( value instanceof Map ) {
            return value
                .findAll { k, v -> v != null }
                .collectEntries { k, v -> Map.entry(k, convertPathsToCidReferences(v)) }
        }
        return value
    }

    @Override
    void onFilePublish(Path destination, Path source, Map annotations){
        storePublishedFile( destination, source, annotations)
    }

    protected String getWorkflowRelative(Path path){
        final outputDirAbs = session.outputDir.toAbsolutePath()
        if (path.isAbsolute()) {
            if (path.startsWith(outputDirAbs)) {
                return outputDirAbs.relativize(path).toString()
            } else {
                throw new Exception("Cannot asses the relative path for workflow output $path")
            }
        } else {
            final pathAbs = path.toAbsolutePath()
            if (pathAbs.startsWith(outputDirAbs)) {
                return outputDirAbs.relativize(pathAbs).toString()
            }
            if (path.normalize().getName(0).toString() == "..")
                throw new Exception("Cannot asses the relative path for workflow output $path")
            return path.normalize().toString()
        }

    }

    protected List<Parameter> manageInputs(Map<InParam, Object> inputs, PathNormalizer normalizer) {
        List<Parameter> managedInputs = new LinkedList<Parameter>()
        inputs.forEach{ param, value ->
            final type = param.class.simpleName
            final name = param.name
            if( param instanceof FileInParam )
                managedInputs.add( new Parameter( type, name, manageFileInParam( (List<FileHolder>)value , normalizer) ) )
            else if( !(param instanceof DefaultInParam) )
                managedInputs.add( new Parameter( type, name, value) )
        }
        return managedInputs
    }

    private List<Object> manageFileInParam(List<FileHolder> files, PathNormalizer normalizer){
        final paths = new LinkedList<Object>();
        for( FileHolder it : files ) {
            final ref = getSourceReference(it.storePath)
            paths.add(ref ? new DataPath(ref) : new DataPath(
                normalizer.normalizePath(it.storePath),
                new Checksum(CacheHelper.hasher(it.storePath).hash().toString(), "nextflow",
                    CacheHelper.HashMode.DEFAULT().toString().toLowerCase()))
            )
        }
        return paths
    }
}
