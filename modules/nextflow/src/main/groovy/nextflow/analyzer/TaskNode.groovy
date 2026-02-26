package nextflow.analyzer

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.dag.DAG
import nextflow.processor.TaskProcessor

/**
 * TaskNode represents a process in the workflow DAG for analysis purposes.
 * It holds references to the DAG vertex, task processor, resource requirements, and dependencies.
 * This class will be used as the node representation in the DependencyGraph for Slurm task grouping analysis.
 * 
 * @author Beatriz Cepa
 */
@Slf4j
@CompileStatic
class TaskNode {

    private DAG.Vertex vertex // Reference to the DAG vertex
    private TaskProcessor processor // Reference to the task processor
    private String name // Name of the process from Vertex.label
    private long id // Vertex ID from Vertex.id
    private Map<String, Object> resources // Process resource requirements
    private List<Long> upstreamDependencies // List of vertex ids of upstream process dependencies
    private List<Long> downstreamDependencies // List of vertex ids of downstream process dependencies
    private boolean isParallelizable // Flag indicating if the process is parallelizable based on resources


    TaskNode(DAG.Vertex vertex, TaskProcessor processor, Map<String, Object> resources) {
        this.vertex = vertex
        this.processor = processor
        this.name = vertex.label
        this.id = vertex.id
        this.resources = resources
        this.isParallelizable = false // Default value 
        this.upstreamDependencies = new ArrayList<>()
        this.downstreamDependencies = new ArrayList<>()
    }

    DAG.Vertex getVertex() {
        return vertex
    }

    TaskProcessor getProcessor() {
        return processor
    }

    String getName() {
        return name
    }

    long getId() {
        return id
    }

    Map<String, Object> getResources() {
        return resources
    }

    List<Long> getUpstreamDependencies() {
        return upstreamDependencies
    }

    List<Long> getDownstreamDependencies() {
        return downstreamDependencies
    }

    boolean getIsParallelizable() {
        return isParallelizable
    }

    void setUpstreamDependencies(List<Long> dependencies) {
        this.upstreamDependencies = dependencies
    }

    void setDownstreamDependencies(List<Long> dependencies) {
        this.downstreamDependencies = dependencies
    }

    void setIsParallelizable(boolean parallelizable) {
        this.isParallelizable = parallelizable
    }

    void addUpstreamDependency(long processId) {
        if (upstreamDependencies == null) {
            upstreamDependencies = new ArrayList<>()
        }
        upstreamDependencies.add(processId)
    }

    void addDownstreamDependency(long processId) {
        if (downstreamDependencies == null) {
            downstreamDependencies = new ArrayList<>()
        }
        downstreamDependencies.add(processId)
    }

}