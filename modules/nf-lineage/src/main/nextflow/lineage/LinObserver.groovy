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
 */

package nextflow.lineage

import nextflow.util.SecretHelper

import java.time.OffsetDateTime

import static nextflow.lineage.fs.LinPath.*

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.lineage.model.Annotation
import nextflow.lineage.model.Checksum
import nextflow.lineage.model.DataOutput
import nextflow.lineage.model.DataPath
import nextflow.lineage.model.Parameter
import nextflow.lineage.model.TaskOutputs
import nextflow.lineage.model.Workflow
import nextflow.lineage.model.WorkflowOutputs
import nextflow.lineage.model.WorkflowRun
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.ScriptMeta

import nextflow.script.params.BaseParam
import nextflow.script.params.CmdEvalParam
import nextflow.script.params.DefaultInParam
import nextflow.script.params.EachInParam
import nextflow.script.params.EnvInParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InParam
import nextflow.script.params.OutParam
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.ValueOutParam
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
import nextflow.util.PathNormalizer
import nextflow.util.TestOnly

/**
 * Observer to write the generated workflow metadata in a lineage store.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LinObserver implements TraceObserver {
    private static Map<Class<? extends BaseParam>, String> TaskParamToValue = [
        (StdOutParam)  : "stdout",
        (StdInParam)   : "stdin",
        (FileInParam)  : "path",
        (FileOutParam) : "path",
        (ValueInParam) : "val",
        (ValueOutParam): "val",
        (EnvInParam)   : "env",
        (EnvOutParam)  : "env",
        (CmdEvalParam) : "eval",
        (EachInParam)  : "each"
    ]
    private String executionHash
    private LinStore store
    private Session session
    private WorkflowOutputs workflowOutputs
    private Map<String,String> outputsStoreDirLid = new HashMap<String,String>(10)
    private PathNormalizer normalizer

    LinObserver(Session session, LinStore store){
        this.session = session
        this.store = store
    }

    @Override
    void onFlowCreate(Session session) {
        this.store.getHistoryLog().write(session.runName, session.uniqueId, '-')
    }

    @TestOnly
    String getExecutionHash(){ executionHash }

    @TestOnly
    String setExecutionHash(String hash){ this.executionHash = hash }

    @TestOnly
    String setNormalizer(PathNormalizer normalizer){  this.normalizer = normalizer }

    @Override
    void onFlowBegin() {
        normalizer = new PathNormalizer(session.workflowMetadata)
        executionHash = storeWorkflowRun(normalizer)
        final executionUri = asUriString(executionHash)
        workflowOutputs = new WorkflowOutputs(
            OffsetDateTime.now(),
            executionUri,
            new LinkedList<Parameter>()
        )
        this.store.getHistoryLog().updateRunLid(session.uniqueId, executionUri)
    }

    @Override
    void onFlowComplete(){
        if (this.workflowOutputs){
            workflowOutputs.createdAt = OffsetDateTime.now()
            final key = executionHash + '#outputs'
            this.store.save(key, workflowOutputs)
        }
    }

    protected Collection<Path> allScriptFiles() {
        return ScriptMeta.allScriptNames().values()
    }

    protected List<DataPath> collectScriptDataPaths(PathNormalizer normalizer) {
        final allScripts = allScriptFiles()
        final result = new ArrayList<DataPath>(allScripts.size()+1)
        // the main script
        result.add( new DataPath(
            normalizer.normalizePath(session.workflowMetadata.scriptFile.normalize()),
            Checksum.of(session.workflowMetadata.scriptId, "nextflow", CacheHelper.HashMode.DEFAULT())
        ) )

        // all other scripts
        for (Path it: allScripts) {
            if( it==null || it == session.workflowMetadata.scriptFile )
                continue
            final dataPath = new DataPath(normalizer.normalizePath(it.normalize()), Checksum.ofNextflow(it.text))
            result.add(dataPath)
        }
        return result
    }

    protected String storeWorkflowRun(PathNormalizer normalizer) {
        // create the workflow object holding script files and repo tracking info
        final workflow = new Workflow(
            collectScriptDataPaths(normalizer),
            session.workflowMetadata.repository,
            session.workflowMetadata.commitId
        )
        // create the workflow run main object
        final value = new WorkflowRun(
            workflow,
            session.uniqueId.toString(),
            session.runName,
            getNormalizedParams(session.params, normalizer),
            SecretHelper.hideSecrets(session.config.deepClone()) as Map
        )
        final executionHash = CacheHelper.hasher(value).hash().toString()
        store.save(executionHash, value)
        return executionHash
    }

    protected static List<Parameter> getNormalizedParams(Map<String, Object> params, PathNormalizer normalizer){
        final normalizedParams = new LinkedList<Parameter>()
        params.each{ String key, Object value ->
            normalizedParams.add( new Parameter( getParameterType(value), key, normalizeValue(value, normalizer) ) )
        }
        return normalizedParams
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        storeTaskInfo(handler.task)
    }

    protected void storeTaskInfo(TaskRun task) {
        // store the task run entry
        storeTaskRun(task, normalizer)
        // store all task results
        storeTaskResults(task, normalizer)
    }

    protected String storeTaskResults(TaskRun task, PathNormalizer normalizer){
        final outputParams = getNormalizedTaskOutputs(task, normalizer)
        final value = new TaskOutputs( asUriString(task.hash.toString()), asUriString(executionHash), OffsetDateTime.now(), outputParams )
        final key = task.hash.toString() + '#outputs'
        store.save(key,value)
        return key
    }

    private List<Parameter> getNormalizedTaskOutputs( TaskRun task, PathNormalizer normalizer){
        final outputs = task.getOutputs()
        final outputParams = new LinkedList<Parameter>()
        outputs.forEach { OutParam key, Object value ->
            manageTaskOutputParameter(key, outputParams, value, task, normalizer)
        }
        return outputParams
    }

    private void manageTaskOutputParameter(OutParam key, LinkedList<Parameter> outputParams, value, TaskRun task, PathNormalizer normalizer) {
        if (key instanceof FileOutParam) {
            outputParams.add(new Parameter(getParameterType(key), key.name, manageFileOutParam(value, task)))
        } else {
            outputParams.add(new Parameter(getParameterType(key), key.name, normalizeValue(value, normalizer)))
        }
    }

    private static Object normalizeValue(Object value, PathNormalizer normalizer) {
        if (value instanceof Path)
            return normalizer.normalizePath(value as Path)
        else if (value instanceof CharSequence)
            return normalizer.normalizePath(value.toString())
        else
            return value
    }

    private Object manageFileOutParam(Object value, TaskRun task) {
        if (value == null) {
            throw new IllegalArgumentException("Unexpected output null for task '${task.name}'")
        }
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
        throw new IllegalArgumentException("Unexpected output [${value.getClass().getName()}] '${value}' for task '${task.name}'")
    }

    protected String storeTaskRun(TaskRun task, PathNormalizer normalizer) {
        final codeChecksum = Checksum.ofNextflow(session.stubRun ? task.stubSource: task.source)
        final scriptChecksum = Checksum.ofNextflow(task.script)
        final value = new nextflow.lineage.model.TaskRun(
            session.uniqueId.toString(),
            task.getName(),
            codeChecksum,
            scriptChecksum,
            task.inputs ? manageTaskInputParameters(task.inputs, normalizer): null,
            task.isContainerEnabled() ? task.getContainerFingerprint(): null,
            normalizer.normalizePath(task.getCondaEnv()),
            normalizer.normalizePath(task.getSpackEnv()),
            task.config?.getArchitecture()?.toString(),
            task.processor.getTaskGlobalVars(task),
            task.processor.getTaskBinEntries(task.source).collect { Path p -> new DataPath(
                normalizer.normalizePath(p.normalize()),
                Checksum.ofNextflow(p) )
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
            final checksum = Checksum.ofNextflow(path)
            final value = new DataOutput(
                path.toUriString(),
                checksum,
                asUriString(task.hash.toString()),
                asUriString(executionHash),
                asUriString(task.hash.toString()),
                attrs.size(),
                LinUtils.toDate(attrs?.creationTime()),
                LinUtils.toDate(attrs?.lastModifiedTime()))
            store.save(key, value)
            return key
        } catch (Throwable e) {
            log.warn("Unexpected error storing lineage output '${path.toUriString()}' for task '${task.name}'", e)
            return path.toUriString()
        }
    }

    protected String getTaskOutputKey(TaskRun task, Path path) {
        final rel = getTaskRelative(task, path)
        return task.hash.toString() + SEPARATOR + rel
    }

    protected String getWorkflowOutputKey(Path destination) {
        final rel = getWorkflowRelative(destination)
        return executionHash + SEPARATOR + rel
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
            this.outputsStoreDirLid.put(path.toString(), asUriString(task.hash.toString(),rel.toString()))
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
            final checksum = Checksum.ofNextflow(destination)
            final key = getWorkflowOutputKey(destination)
            final sourceReference = source ? getSourceReference(source) : asUriString(executionHash)
            final attrs = readAttributes(destination)
            final value = new DataOutput(
                destination.toUriString(),
                checksum,
                sourceReference,
                asUriString(executionHash),
                null,
                attrs.size(),
                LinUtils.toDate(attrs?.creationTime()),
                LinUtils.toDate(attrs?.lastModifiedTime()),
                convertAnnotations(annotations))
            store.save(key, value)
        } catch (Throwable e) {
            log.warn("Unexpected error storing published file '${destination.toUriString()}' for workflow '${executionHash}'", e)
        }
    }

    private static List<Annotation> convertAnnotations(Map annotations){
        if( !annotations )
            return null
        final converted = new LinkedList<Annotation>()
        annotations.forEach { Object key, Object value -> converted.add(new Annotation(key.toString(), value)) }
        return converted
    }
    String getSourceReference(Path source){
        final hash = FileHelper.getTaskHashFromPath(source, session.workDir)
        if (hash) {
            final target = FileHelper.getWorkFolder(session.workDir, hash).relativize(source).toString()
            return asUriString(hash.toString(), target)
        }
        final storeDirReference = outputsStoreDirLid.get(source.toString())
        return storeDirReference ? asUriString(storeDirReference) : null
    }

    @Override
    void onFilePublish(Path destination){
        storePublishedFile (destination)
    }

    @Override
    void onWorkflowPublish(String name, Object value){
        workflowOutputs.outputs.add(new Parameter(getParameterType(value), name, convertPathsToLidReferences(value)))
    }

    protected static String getParameterType(Object param) {
        if( param instanceof BaseParam )
            return TaskParamToValue.get(param.class)
        // return generic types
        if( param instanceof Path )
            return Path.simpleName
        if (param instanceof CharSequence)
            return String.simpleName
        if( param instanceof Collection )
            return Collection.simpleName
        if( param instanceof Map)
            return Map.simpleName
        return param.class.simpleName
    }

    private Object convertPathsToLidReferences(Object value){
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
            return value.collect { el -> convertPathsToLidReferences(el) }
        }

        if( value instanceof Map ) {
            return value
                .findAll { k, v -> v != null }
                .collectEntries { k, v -> Map.entry(k, convertPathsToLidReferences(v)) }
        }
        return value
    }

    @Override
    void onFilePublish(Path destination, Path source, Map annotations){
        storePublishedFile( destination, source, annotations)
    }
    /**
     * Relativizes a path from the workflow's output dir.
     *
     * @param path Path to relativize
     * @return Path String with the relative path
     * @throws IllegalArgumentException
     */
    protected String getWorkflowRelative(Path path) throws IllegalArgumentException{
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

    protected List<Parameter> manageTaskInputParameters(Map<InParam, Object> inputs, PathNormalizer normalizer) {
        List<Parameter> managedInputs = new LinkedList<Parameter>()
        inputs.forEach{ param, value ->
            if( param instanceof FileInParam )
                managedInputs.add( new Parameter( getParameterType(param), param.name, manageFileInParam( (List<FileHolder>)value , normalizer) ) )
            else if( !(param instanceof DefaultInParam) )
                managedInputs.add( new Parameter( getParameterType(param), param.name, value) )
        }
        return managedInputs
    }

    private List<Object> manageFileInParam(List<FileHolder> files, PathNormalizer normalizer){
        final paths = new LinkedList<Object>();
        for( FileHolder it : files ) {
            final ref = getSourceReference(it.storePath)
            paths.add(ref ? ref : new DataPath(
                normalizer.normalizePath(it.storePath),
                Checksum.ofNextflow(it.storePath))
            )
        }
        return paths
    }
}
