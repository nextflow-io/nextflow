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

    enum ProcessingState {
    UNVISITED,    // Not yet processed
    VISITING,     // Currently being processed 
    VISITED       // Fully processed
    }

    private ProcessingState state // State for cycle detection during level assignment
    private DAG.Vertex vertex // Reference to the DAG vertex
    private TaskProcessor processor // Reference to the task processor
    private String name // Name of the process from Vertex.label
    private long id // Vertex ID from Vertex.id
    private Map<String, Object> resources // Process resource requirements
    private List<Long> upstreamDependencies // List of vertex ids of upstream process dependencies
    private List<Long> downstreamDependencies // List of vertex ids of downstream process dependencies
    private boolean isParallelizable // Flag indicating if the process is parallelizable based on resources
    private int level // Dependency level in the DAG 

    TaskNode(DAG.Vertex vertex, TaskProcessor processor, Map<String, Object> resources) {
        this.vertex = vertex
        this.processor = processor
        this.name = vertex.label
        this.id = vertex.id
        this.resources = resources
        this.isParallelizable = false // Default value 
        this.upstreamDependencies = new ArrayList<>()
        this.downstreamDependencies = new ArrayList<>()
        this.level = 0
        this.state = ProcessingState.UNVISITED
    }

    DAG.Vertex getVertex() {
        return this.vertex
    }

    TaskProcessor getProcessor() {
        return this.processor
    }

    String getName() {
        return this.name
    }

    long getId() {
        return this.id
    }

    Map<String, Object> getResources() {
        return this.resources
    }

    List<Long> getUpstreamDependencies() {
        return this.upstreamDependencies
    }

    List<Long> getDownstreamDependencies() {
        return this.downstreamDependencies
    }

    boolean getIsParallelizable() {
        return this.isParallelizable
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

    int getLevel() {
        return this.level
    }

    void setLevel(int level) {
        this.level = level
    }

    ProcessingState getState() {
        return this.state
    }

    void setState(ProcessingState state) {
        this.state = state
    }

    /**
     * Add an upstream dependency by vertex ID
     * @param vertexId the vertex ID of the upstream dependency
     */
    void addUpstreamDependency(long vertexId) {
        if (this.upstreamDependencies == null) {
            this.upstreamDependencies = new ArrayList<>()
        }
        this.upstreamDependencies.add(vertexId)
    }

    /**
     * Add a downstream dependency by vertex ID
     * @param vertexId the vertex ID of the downstream dependency
     */
    void addDownstreamDependency(long vertexId) {
        if (this.downstreamDependencies == null) {
            this.downstreamDependencies = new ArrayList<>()
        }
        this.downstreamDependencies.add(vertexId)
    }
}