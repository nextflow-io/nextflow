/*
 * Copyright 2013-2025, Seqera Labs
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

import static nextflow.data.cid.fs.CidPath.*

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.time.Instant

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataOutput
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.TaskOutputs
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowOutputs
import nextflow.data.cid.model.WorkflowRun
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.ScriptMeta
import nextflow.script.params.DefaultInParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InParam
import nextflow.script.params.OutParam
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
import nextflow.util.PathNormalizer
import nextflow.util.TestOnly
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
    private WorkflowOutputs workflowResults
    private Map<String,String> outputsStoreDirCid = new HashMap<String,String>(10)

    CidObserver(Session session, CidStore store){
        this.session = session
        this.store = store
    }

    @Override
    void onFlowCreate(Session session) {
        this.store.getHistoryLog().write(session.runName, session.uniqueId, '-')
    }

    @TestOnly
    String getExecutionHash(){ executionHash }

    @Override
    void onFlowBegin() {
        executionHash = storeWorkflowRun()
        final executionUri = asUriString(executionHash)
        workflowResults = new WorkflowOutputs(
            Instant.now(),
            executionUri,
            new HashMap<String, Object>()
        )
        this.store.getHistoryLog().updateRunCid(session.uniqueId, executionUri)
    }

    @Override
    void onFlowComplete(){
        if (this.workflowResults){
            workflowResults.createdAt = Instant.now()
            final key = executionHash + SEPARATOR + 'outputs'
            this.store.save(key, workflowResults)
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
            mainScript,
            otherScripts,
            session.workflowMetadata.repository,
            session.workflowMetadata.commitId
        )
        final value = new WorkflowRun(
            workflow,
            session.uniqueId.toString(),
            session.runName,
            getNormalizedParams(session.params, normalizer)
        )
        final executionHash = CacheHelper.hasher(value).hash().toString()
        store.save(executionHash, value)
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
        // store all task results
        storeTaskResults(task, pathNormalizer)
    }

    protected String storeTaskResults(TaskRun task, PathNormalizer normalizer){
        final outputParams = getNormalizedTaskOutputs(task, normalizer)
        final value = new TaskOutputs( asUriString(task.hash.toString()), asUriString(executionHash), Instant.now(), outputParams )
        final key = CacheHelper.hasher(value).hash().toString()
        store.save(key,value)
        return key
    }

    private List<Parameter> getNormalizedTaskOutputs( TaskRun task, PathNormalizer normalizer){
        final outputs = task.getOutputs()
        final outputParams = new LinkedList<Parameter>()
        outputs.forEach { OutParam key, Object value ->
            if (key instanceof FileOutParam) {
                outputParams.add( new Parameter( key.class.simpleName, key.name, manageFileOutParams(value, task) ) )
            }
            else {
                if( value instanceof Path )
                    outputParams.add( new Parameter( key.class.simpleName, key.name, normalizer.normalizePath( value as Path ) ) )
                else if ( value instanceof CharSequence )
                    outputParams.add( new Parameter( key.class.simpleName, key.name, normalizer.normalizePath( value.toString() ) ) )
                else
                    outputParams.add( new Parameter( key.class.simpleName, key.name, value) )
            }
        }
        return outputParams
    }

    private Object manageFileOutParams(Object value, TaskRun task) {
        if (value instanceof Path) {
            return asUriString(storeTaskOutput(task, (Path) value))
        }
        if (value instanceof Collection<Path>) {
            final files = new LinkedList<String>()
            for (Path it : value) {
                files.add( asUriString(storeTaskOutput(task, (Path)it)) )
            }
            return files
        }
        // unexpected task output
        final msg = value!=null
            ? "Unexepected output [${value.getClass().getName()}] '${value}' for task '${task.name}'"
            : "Unexpected output null for task '${task.name}'"
        throw new IllegalArgumentException(msg)
    }

    protected String storeTaskRun(TaskRun task, PathNormalizer normalizer) {
        final codeChecksum = new Checksum(CacheHelper.hasher(session.stubRun ? task.stubSource: task.source).hash().toString(),
            "nextflow", CacheHelper.HashMode.DEFAULT().toString().toLowerCase())
        final scriptChecksum = new Checksum(CacheHelper.hasher(task.script).hash().toString(),
            "nextflow", CacheHelper.HashMode.DEFAULT().toString().toLowerCase())
        final value = new nextflow.data.cid.model.TaskRun(
            session.uniqueId.toString(),
            task.getName(),
            codeChecksum,
            scriptChecksum,
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
            },
            asUriString(executionHash)
        )

        // store in the underlying persistence
        final key = task.hash.toString()
        store.save(key, value)
        return key
    }

    protected String storeTaskOutput(TaskRun task, Path path) {
        try {
            final attrs = readAttributes(path)
            final key = getTaskOutputKey(task, path)
            final checksum = new Checksum( CacheHelper.hasher(path).hash().toString(),
                "nextflow", CacheHelper.HashMode.DEFAULT().toString().toLowerCase() )
            final value = new DataOutput(
                path.toUriString(),
                checksum,
                asUriString(task.hash.toString()),
                asUriString(task.hash.toString()),
                attrs.size(),
                CidUtils.toDate(attrs?.creationTime()),
                CidUtils.toDate(attrs?.lastModifiedTime()))
            store.save(key, value)
            return key
        } catch (Throwable e) {
            log.warn("Unexpected error storing CID output '${path.toUriString()}' for task '${task.name}'", e)
            return path.toUriString()
        }
    }

    protected String getTaskOutputKey(TaskRun task, Path path) {
        final rel = getTaskRelative(task, path)
        return task.hash.toString() + SEPARATOR + 'outputs' + SEPARATOR + rel
    }

    protected String getWorkflowOutputKey(Path destination) {
        final rel = getWorkflowRelative(destination)
        return executionHash + SEPARATOR + 'outputs' + SEPARATOR + rel
    }

    protected String getTaskRelative(TaskRun task, Path path){
        if (path.isAbsolute()) {
            final rel = getTaskRelative0(task, path)
            if (rel)
                return rel
            throw new IllegalArgumentException("Cannot access the relative path for output '${path.toUriString()}' and task '${task.name}'")
        }
        //Check if contains workdir or storeDir
        final rel = getTaskRelative0(task, path.toAbsolutePath())
        if (rel) return rel
        if (path.normalize().getName(0).toString() == "..")
            throw new IllegalArgumentException("Cannot access the relative path for output '${path.toUriString()}' and task '${task.name}'" )
        return path.normalize().toString()
    }

    private String getTaskRelative0(TaskRun task, Path path){
        final workDirAbsolute = task.workDir.toAbsolutePath()
        if (path.startsWith(workDirAbsolute)) {
            return workDirAbsolute.relativize(path).toString()
        }
        //If task output is not in the workDir check if output is stored in the task's storeDir
        final storeDir = task.getConfig().getStoreDir().toAbsolutePath()
        if( storeDir && path.startsWith(storeDir) ) {
            final rel = storeDir.relativize(path)
            //If output stored in storeDir, keep the path in case it is used as workflow output
            this.outputsStoreDirCid.put(path.toString(), asUriString(task.hash.toString(),rel.toString()))
            return rel
        }
        return null
    }

    protected BasicFileAttributes readAttributes(Path path) {
        return Files.readAttributes(path, BasicFileAttributes)
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
            final key = getWorkflowOutputKey(destination)
            final sourceReference = source ? getSourceReference(source) : asUriString(executionHash)
            final attrs = readAttributes(destination)
            final value = new DataOutput(
                destination.toUriString(),
                checksum,
                sourceReference,
                asUriString(executionHash),
                attrs.size(),
                CidUtils.toDate(attrs?.creationTime()),
                CidUtils.toDate(attrs?.lastModifiedTime()),
                annotations)
            store.save(key, value)
        } catch (Throwable e) {
            log.warn("Unexpected error storing published file '${destination.toUriString()}' for workflow '${executionHash}'", e)
        }
    }

    String getSourceReference(Path source){
        final hash = FileHelper.getTaskHashFromPath(source, session.workDir)
        if (hash) {
            final target = FileHelper.getWorkFolder(session.workDir, hash).relativize(source).toString()
            return asUriString(hash.toString(), 'outputs', target)
        }
        final storeDirReference = outputsStoreDirCid.get(source.toString())
        return storeDirReference ? asUriString(storeDirReference) : null
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
            try {
                final key = getWorkflowOutputKey(value)
                return asUriString(key)
            } catch (Throwable e){
                //Workflow output key not found
                return value
            }
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
            }
            throw new IllegalArgumentException("Cannot access relative path for workflow output '${path.toUriString()}'")
        }
        final pathAbs = path.toAbsolutePath()
        if (pathAbs.startsWith(outputDirAbs)) {
            return outputDirAbs.relativize(pathAbs).toString()
        }
        if (path.normalize().getName(0).toString() == "..")
            throw new IllegalArgumentException("Cannot access relative path for workflow output '${path.toUriString()}'")
        return path.normalize().toString()
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
