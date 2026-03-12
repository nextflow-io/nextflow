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
     * Main method to perform analysis for Slurm task grouping.
     * It checks if task grouping is enabled, collects process vertices, builds a dependency graph,
     * and identifies parallelizable tasks based on their levels in the graph.
     */
    void analyze() {
        if( !isTaskGroupingEnabled() ) {
            log.debug "[SLURM TASK GROUPING] Task grouping is disabled"
            return
        }
        log.debug "[SLURM TASK GROUPING] Task grouping enabled. Starting Slurm Analyzer"
        try {
            if( this.dag == null ) {
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
            Map<Integer, List<Long>> parallelTasks = identifyParallelTasks(dependencyGraph)
            log.debug "[SLURM TASK GROUPING] Tasks grouped by levels: ${parallelTasks.collect { k, v -> "Level $k: ${v.size()} tasks" }.join(', ')}"
        

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
        return session.config.navigate('executor.slurm.taskGrouping', false) as boolean
    }

    /**
     * Collect all process vertices from the DAG
     * @return list of process vertices
     */
    private List<DAG.Vertex> collectProcessVertices() {
        final List<DAG.Vertex> result = []
        for (DAG.Vertex vertex : this.dag.vertices) {
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
        DependencyGraph dependencyGraph = new DependencyGraph(this.session)

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
        final TaskNode node = new TaskNode(vertex, processor)
        log.debug "[SLURM TASK GROUPING] Analyzing vertex: ${vertex.label} — cpus=${node.getCpus()}, memory=${node.getMemory()}, time=${node.getTime()}"
        return node
    } 

    /**
     * Identify parallelizable tasks based on dependency graph levels.
     * Tasks at the same level have no dependency relationship with 
     * each other and can potentially be executed in parallel.
     * Uses iterative topological sort (Kahn's algorithm) for O(V+E) performance.
     *
     * @param dependencyGraph the graph representing process dependencies
     * @return a map of level to list of TaskNode IDs that can be grouped together
     */
    Map<Integer, List<Long>> identifyParallelTasks(DependencyGraph dependencyGraph) {
        // Use iterative topological sort for better performance
        dependencyGraph.assignLevelsIteratively()
        
        final Map<Integer, List<Long>> groups = dependencyGraph.groupByLevels()
        
        // Log summary of grouping results
        for (Map.Entry<Integer, List<Long>> entry : groups.entrySet()) {
            log.debug "[SLURM TASK GROUPING] Level ${entry.getKey()}: ${entry.getValue().size()} task(s)"
        }
        
        return groups
    }

}