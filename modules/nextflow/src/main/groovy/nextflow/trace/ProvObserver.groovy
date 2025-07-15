/*
 * Copyright 2022, Seqera Labs
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

package nextflow.trace

import java.nio.file.Path
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.params.FileOutParam
import nextflow.script.WorkflowMetadata
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord
import nextflow.util.StringUtils

/**
 * Plugin observer of workflow events
 *
 * @author Bruno Grande <bruno.grande@sagebase.org>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ProvObserver implements TraceObserver {

    private Session session

    private Set<TaskRun> tasks = []

    private Map<Path,Path> workflowOutputs = [:]

    private Lock lock = new ReentrantLock()

    @Override
    void onFlowCreate(Session session) {
        this.session = session
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        // skip failed tasks
        final task = handler.task
        if( !task.isSuccess() )
            return

        lock.withLock {
            tasks << task
        }
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {
        lock.withLock {
            tasks << handler.task
        }
    }

    @Override
    void onFilePublish(Path target, Path source) {
        lock.withLock {
            workflowOutputs[source] = target
        }
    }

    @Override
    void onFlowComplete() {
        if( !session.isSuccess() )
            return

        try {
            new TaskDagRenderer().render(session, tasks, workflowOutputs)
        }
        catch( Exception e ) {
            log.warn "Exception while rendering provenance report: ${e}"
            throw e
        }
    }

}


/**
 * Renderer for the task graph format.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskDagRenderer {

    private Path path

    @Delegate
    private PathNormalizer normalizer

    TaskDagRenderer() {
        path = ('dag.html' as Path).complete()
    }

    void render(Session session, Set<TaskRun> tasks, Map<Path,Path> workflowOutputs) {
        // get workflow metadata
        final metadata = session.workflowMetadata
        this.normalizer = new PathNormalizer(metadata)

        // get workflow inputs
        final taskLookup = ProvHelper.getTaskLookup(tasks)
        final workflowInputs = ProvHelper.getWorkflowInputs(tasks, taskLookup)

        // construct task graph
        final dag = new Dag(getVertices(tasks), taskLookup)

        path.text = renderHtml(dag)
    }

    private Map<TaskRun,Vertex> getVertices(Set<TaskRun> tasks) {
        Map<TaskRun,Vertex> result = [:]
        for( final task : tasks ) {
            final inputs = task.getInputFilesMap()
            final outputs = ProvHelper.getTaskOutputs(task)

            result[task] = new Vertex(result.size(), task.name, inputs, outputs)
        }

        return result
    }

    /**
     * Render the task graph as a Mermaid diagram embedded
     * in an HTML document.
     *
     * @param dag
     */
    private String renderHtml(Dag dag) {
        // load html template
        final writer = new StringWriter()
        final res = nextflow.dag.DagRenderer.class.getResourceAsStream('mermaid.dag.template.html')
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char)
        }
        final template = writer.toString()

        // render html document
        final mmd = renderDiagram(dag)
        return template.replace('REPLACE_WITH_NETWORK_DATA', mmd)
    }

    /**
     * Render the task graph as a Mermaid diagram.
     *
     * @param dag
     */
    private String renderDiagram(Dag dag) {
        // construct task tree
        final taskTree = getTaskTree(dag.vertices)

        // render diagram
        List<String> lines = []
        lines << "flowchart TD"

        // render workflow inputs
        Map<Path,String> inputs = [:]

        lines << "    subgraph \" \""

        dag.vertices.each { task, vertex ->
            vertex.inputs.each { name, path ->
                if( dag.getSourceVertex(path) || path in inputs )
                    return

                inputs[path] = "in${inputs.size()}".toString()
                lines << "    ${inputs[path]}[\"${path.name}\"]".toString()

                // add hyperlink if path is remote URL
                final sourcePath = normalizePath(path)
                if( isRemotePath(sourcePath) )
                    lines << "    click ${inputs[path]} href \"${sourcePath}\" _blank".toString()
            }
        }

        lines << "    end"

        // render tasks
        renderTaskTree(lines, null, taskTree)

        // render task inputs
        Set<Path> taskOutputs = []

        dag.vertices.each { task, vertex ->
            // render inputs from upstream tasks
            task.upstreamTasks.each { id ->
                final upstreamTask = dag.vertices.keySet().find { t -> t.id == id }
                final pred = dag.vertices[upstreamTask]
                if( pred != null ) {
                    final files = (pred.outputs + vertex.inputs.values()) as Set<Path>
                    if( !files )
                        lines << "    ${pred.name} --> ${vertex.name}".toString()
                    for( final file : files )
                        lines << "    ${pred.name} -->|${file.name}| ${vertex.name}".toString()
                    taskOutputs.addAll(files)
                }
            }

            // render inputs from workflow inputs
            vertex.inputs.each { name, path ->
                final pred = dag.getSourceVertex(path)
                if( !pred )
                    lines << "    ${inputs[path]} --> ${vertex.name}".toString()
            }
        }

        // render task outputs
        Map<Path,String> outputs = [:]

        dag.vertices.each { task, vertex ->
            vertex.outputs.each { path ->
                if( path in taskOutputs || path in outputs )
                    return

                outputs[path] = "out${outputs.size()}".toString()
                lines << "    ${vertex.name} --> ${outputs[path]}".toString()
            }
        }

        // render workflow outputs
        lines << "    subgraph \" \""

        outputs.each { path, label ->
            lines << "    ${label}[${path.name}]".toString()
        }

        lines << "    end"

        return lines.join('\n')
    }

    /**
     * Construct a task tree with a subgraph for each subworkflow.
     *
     * @param vertices
     */
    private Map getTaskTree(Map<TaskRun,Vertex> vertices) {
        final taskTree = [:]

        for( final entry : vertices ) {
            final task = entry.key
            final vertex = entry.value

            // infer subgraph keys from fully qualified process name
            final result = getSubgraphKeys(task.processor.name)
            final keys = (List)result[0]

            // update vertex label
            final hash = task.hash.toString()
            vertex.label = task.name.replace(task.processor.name, (String)result[1])

            // navigate to given subgraph
            def subgraph = taskTree
            for( final key : keys ) {
                if( key !in subgraph )
                    subgraph[key] = [:]
                subgraph = subgraph[key]
            }

            // add vertex to tree
            subgraph[vertex.name] = vertex
        }

        return taskTree
    }

    /**
     * Get the subgraph keys from a fully qualified process name.
     *
     * @param name
     */
    private List getSubgraphKeys(String name) {
        final tokens = name.tokenize(':')
        return [
            tokens[0..<tokens.size()-1],
            tokens.last()
        ]
    }

    private boolean isRemotePath(String path) {
        if( !path ) return false
        final result = StringUtils.getUrlProtocol(path)
        return result != null && result != 'file'
    }

    /**
     * Render a tree of tasks and subgraphs.
     *
     * @param lines
     * @param name
     * @param taskTree
     */
    private void renderTaskTree(List<String> lines, String name, Map<String,Object> taskTree) {
        if( name )
            lines << "    subgraph ${name}".toString()

        taskTree.each { key, value ->
            if( value instanceof Map )
                renderTaskTree(lines, key, value)
            else if( value instanceof Vertex )
                lines << "    ${renderTask(value)}".toString()
        }

        if( name )
            lines << "    end"
    }

    private String renderTask(Vertex vertex) {
        return "${vertex.name}([\"${vertex.label}\"])"
    }

    @TupleConstructor
    static class Dag {
        Map<TaskRun,Vertex> vertices
        Map<Path,TaskRun> taskLookup

        Vertex getSourceVertex(Path path) {
            vertices[taskLookup[path]]
        }
    }

    @TupleConstructor
    static class Vertex {
        int index
        String label
        Map<String,Path> inputs
        Set<Path> outputs

        String getName() { "t${index}" }

        @Override
        String toString() { label }
    }

}


/**
 * Helper methods for provenance reports.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProvHelper {

    /**
     * Get the set of output files for a task.
     *
     * @param task
     */
    static Set<Path> getTaskOutputs(TaskRun task) {
        return task
            .getOutputsByType(FileOutParam)
            .values()
            .flatten() as Set<Path>
    }

    /**
     * Get a mapping of output file to the task that produced it.
     *
     * @param tasks
     */
    static Map<Path,TaskRun> getTaskLookup(Set<TaskRun> tasks) {
        Map<Path,TaskRun> result = [:]

        for( def task : tasks )
            for( def output : getTaskOutputs(task) )
                result[output] = task

        return result
    }

    /**
     * Get the list of workflow inputs. A workflow input is an input file
     * to a task that was not produced by another task.
     *
     * @param tasks
     * @param taskLookup
     */
    static Set<Path> getWorkflowInputs(Set<TaskRun> tasks, Map<Path,TaskRun> taskLookup) {
        Set<Path> result = []

        tasks.each { task ->
            task.getInputFilesMap().each { name, path ->
                if( taskLookup[path] )
                    return

                result << path
            }
        }

        return result
    }

}


/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class PathNormalizer {

    private URL repository

    private String commitId

    private String launchDir

    private String projectDir

    private String workDir

    PathNormalizer(WorkflowMetadata metadata) {
        repository = metadata.repository ? new URL(metadata.repository) : null
        commitId = metadata.commitId
        projectDir = metadata.projectDir.toUriString()
        launchDir = metadata.launchDir.toUriString()
        workDir = metadata.workDir.toUriString()
    }

    /**
     * Normalize paths so that local absolute paths become
     * relative paths, and local paths derived from remote URLs
     * become the URLs.
     *
     * @param path
     */
    String normalizePath(Path path) {
        normalizePath(path.toUriString())
    }

    String normalizePath(String path) {
        // replace work directory with relative path
        if( path.startsWith(workDir) )
            return path.replace(workDir, 'work')

        // replace project directory with source URL (if applicable)
        if( repository && path.startsWith(projectDir) )
            return getProjectSourceUrl(path)

        // replace launch directory with relative path
        if( path.startsWith(launchDir) )
            return path.replace(launchDir + '/', '')

        return path
    }

    /**
     * Get the source URL for a project asset.
     *
     * @param path
     */
    private String getProjectSourceUrl(String path) {
        switch( repository.host ) {
        case 'bitbucket.org':
            return path.replace(projectDir, "${repository}/src/${commitId}")
        case 'github.com':
            return path.replace(projectDir, "${repository}/tree/${commitId}")
        case 'gitlab.com':
            return path.replace(projectDir, "${repository}/-/tree/${commitId}")
        default:
            return path
        }
    }

}
