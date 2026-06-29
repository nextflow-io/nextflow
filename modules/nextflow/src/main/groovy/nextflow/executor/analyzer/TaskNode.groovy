package nextflow.executor.analyzer

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.dag.DAG
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.util.Duration
import nextflow.util.MemoryUnit

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
    VISITED       // Fully processed
    }

    private ProcessingState state // State for cycle detection during level assignment
    private DAG.Vertex vertex // Reference to the DAG vertex
    private TaskProcessor processor // Reference to the task processor
    private String name // Name of the process from Vertex.label
    private long id // Vertex ID from Vertex.id
    private int cpus
    private MemoryUnit memory
    private Duration time
    private String queue
    private List<Long> upstreamDependencies // List of vertex ids of upstream process dependencies
    private List<Long> downstreamDependencies // List of vertex ids of downstream process dependencies
    private int level // Dependency level in the DAG 

    TaskNode(DAG.Vertex vertex, TaskProcessor processor) {
        this.vertex = vertex
        this.processor = processor
        this.name = vertex.label
        this.id = vertex.id
        final TaskConfig config = processor.config.createTaskConfig()
        try {
            this.cpus = config.getCpus()
            this.memory = config.getMemory()
            this.time = config.getTime()
            this.queue = config.queue as String
        }
        catch( IllegalStateException e ) {
            log.warn "[SLURM TASK GROUPING] Process '${vertex.label}' uses dynamic resource directives — static grouping will use defaults: ${e.message}"
            if( !this.cpus ) this.cpus = 1
        }
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

    int getCpus() { return this.cpus }

    MemoryUnit getMemory() { return this.memory }

    Duration getTime() { return this.time }

    String getQueue() { return this.queue }

    List<Long> getUpstreamDependencies() {
        return this.upstreamDependencies
    }

    List<Long> getDownstreamDependencies() {
        return this.downstreamDependencies
    }

    void setUpstreamDependencies(List<Long> dependencies) {
        this.upstreamDependencies = dependencies
    }

    void setDownstreamDependencies(List<Long> dependencies) {
        this.downstreamDependencies = dependencies
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
        this.upstreamDependencies.add(vertexId)
    }

    /**
     * Add a downstream dependency by vertex ID
     * @param vertexId the vertex ID of the downstream dependency
     */
    void addDownstreamDependency(long vertexId) {
        this.downstreamDependencies.add(vertexId)
    }
}