package nextflow.analyzer

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.DAG
import nextflow.processor.TaskProcessor

/**
 * Analyzer for Slurm task grouping.
 * 
 * It analyzes the workflow DAG to identify processes that can be grouped together
 * based on static resource requirements and dependency relationships.
 * 
 * @author Beatriz Cepa
 */
@Slf4j
@CompileStatic
class SlurmTaskGroupAnalyzer {

    private Session session
    private DAG dag

    SlurmTaskGroupAnalyzer(Session session) {
        this.session = session
        this.dag = session.dag  // Access the DAG from the session for analysis
    }

    /**
     * Analyze the workflow DAG for task grouping opportunities.
     * 
     * This method:
     * 1. Accesses the DAG to iterate through processes (vertices)
     * 2. Builds dependency graph using DAG edges
     * 3. Identifies processes with static resource requirements
     */
    void analyze() {
        if( !isTaskGroupingEnabled() ) {
            log.debug "[SLURM TASK GROUPING] Task grouping is disabled"
            return
        }
        log.debug "[SLURM TASK GROUPING] Task grouping enabled. Starting Slurm Analyzer"
        try {
            if( dag == null ) {
                log.debug "[SLURM TASK GROUPING] No DAG available for analysis"
                return
            }
            // Collect process vertices from the DAG
            final List<DAG.Vertex> processVertices = collectProcessVertices()
            log.debug "[SLURM TASK GROUPING] Found ${processVertices.size()} processes to analyze"
            
            // Build dependency graph
            DependencyGraph dependencyGraph = buildDependencyGraph(processVertices)
            log.debug "[SLURM TASK GROUPING] Dependency graph built with ${dependencyGraph.getNodes().size()} nodes"
            
            // Identify parallelizable tasks
            List<TaskNode> parallelTasks = identifyParallelTasks(dependencyGraph)
            log.debug "[SLURM TASK GROUPING] Identified ${parallelTasks.size()} parallelizable tasks for grouping"

        }
        catch( Exception e ) {
            log.error "[SLURM TASK GROUPING] Error during analysis: ${e.message}", e
        }
    }

    /**
      * Check if task grouping is enabled via config
      * @return true if enabled, false otherwise
      */
    private boolean isTaskGroupingEnabled() {
         return session.config?.executor?.slurm?.taskGrouping == true
    }

    /**
     * Collect all process vertices from the DAG
     * @return list of process vertices
     */
    private List<DAG.Vertex> collectProcessVertices() {
        final List<DAG.Vertex> result = []
        for (DAG.Vertex vertex : dag.vertices) {
            if (vertex.type == DAG.Type.PROCESS && vertex.process != null) {
                result << vertex
            }
        }
        return result
    }

    /**
     * Build a dependency graph from the DAG for grouping analysis
     * @param processVertices the list of process vertices to include in the graph
     */
    private DependencyGraph buildDependencyGraph(List<DAG.Vertex> processVertices) {
        DependencyGraph dependencyGraph = new DependencyGraph(session)

        // Add all process nodes to the graph
        for( DAG.Vertex vertex : processVertices ) {
            TaskNode node = analyzeVertex(vertex)
            dependencyGraph.addNode(node)
        }
        
        // Build dependencies between nodes based on DAG edges
        dependencyGraph.buildDependencies()
        
        return dependencyGraph
    } 
    
    /**
      * Analyze a single vertex for grouping potential
      * @param vertex the DAG vertex to analyze
      */
    private TaskNode analyzeVertex(DAG.Vertex vertex) {
        final TaskProcessor processor = vertex.process
        log.debug "[SLURM TASK GROUPING] Analyzing vertex: ${vertex.label}"

        // Extract static resource requirements from processor config
        final config = processor.config

        final Integer cpus = config.get('cpus') as Integer
        final def memory = config.get('memory')
        final def time = config.get('time')
        final def disk = config.get('disk')
        final String queue = config.get('queue') as String

        log.debug "[SLURM TASK GROUPING] Vertex ${vertex.label} resources: cpus=${cpus}, memory=${memory}, time=${time}"

        Map<String, Object> resources = [
            cpus: cpus, 
            memory: memory, 
            time: time, 
            disk: disk,
            queue: queue
        ]

        // Create and return TaskNode with vertex reference
        TaskNode node = new TaskNode(vertex, processor, resources)
        return node
    } 

    /**
     * Identify parallelizable tasks based on the dependency graph
     * @param dependencyGraph the graph representing process dependencies
     * @return list of parallelizable TaskNodes
     */
    private List<TaskNode> identifyParallelTasks(DependencyGraph dependencyGraph) {
        List<TaskNode> parallelTasks = new ArrayList<>()
        for (TaskNode node : dependencyGraph.getNodes().values()) {
            if (node.getUpstreamDependencies().isEmpty()) {
                node.setIsParallelizable(true)
                parallelTasks.add(node)
            }
        }
        return parallelTasks
    }

}