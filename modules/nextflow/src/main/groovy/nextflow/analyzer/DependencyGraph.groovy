package nextflow.analyzer

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.DAG

/**
 * DependencyGraph represents the relationships between processes in the workflow DAG for analysis purposes.
 * It builds a graph structure where each node is a TaskNode representing a process, and edges represent dependencies.
 * This graph will be used to identify parallelizable tasks for Slurm task grouping.
 * 
 * @author Beatriz Cepa
 */
@Slf4j
@CompileStatic
class DependencyGraph {

    private Map<Long, TaskNode> nodes // Map of vertex ID to TaskNode representing each process in the DAG
    private DAG dag // Reference to the workflow DAG for analysis
    private Session session 

    DependencyGraph(Session session) {
        this.nodes = new HashMap<>()
        this.session = session
        this.dag = session.dag
    }

    Map<Long, TaskNode> getNodes() {
        return this.nodes
    }

    DAG getDag() {
        return this.dag
    }

    /**
     * Add a TaskNode to the graph
     * @param node the TaskNode to add
     */
    void addNode(TaskNode node) {
        this.nodes.put(node.getId(), node)
    }

    /**
     * Build dependencies between nodes based on DAG edges.
     * Only considers edges between PROCESS vertices.
     */
    void buildDependencies() {
        for (DAG.Edge edge : this.dag.edges) {
            // Only process edges between actual process vertices
            if (edge.from?.type == DAG.Type.PROCESS && edge.to?.type == DAG.Type.PROCESS) {
                final TaskNode fromNode = this.nodes.get(edge.from.id)
                final TaskNode toNode = this.nodes.get(edge.to.id)
                
                if (fromNode != null && toNode != null) {
                    toNode.addUpstreamDependency(edge.from.id)
                    fromNode.addDownstreamDependency(edge.to.id)
                }
            }
        }
    }

    /**
     * Assign levels to each node in the graph based on dependencies using iterative topological sort.
     * This method uses Kahn's algorithm for better performance O(V+E) and built-in cycle detection.
     * Level 0 means no dependencies, level 1 means depends on level 0, etc.
     */
    void assignLevelsIteratively() {
        // Queue for nodes with no unprocessed upstream dependencies (incoming degree of 0)
        final Queue<TaskNode> queue = new LinkedList<>()
        
        // Track incoming degree (number of unprocessed upstream dependencies)
        final Map<Long, Integer> inDegree = new HashMap<>()
        
        // Initialize incoming degrees and find initial nodes with no dependencies
        for (TaskNode node : this.nodes.values()) {
            final int degree = node.getUpstreamDependencies().size()
            inDegree.put(node.getId(), degree)
            
            // Enqueue nodes with no dependencies and set their level to 0
            if (degree == 0) {
                node.setLevel(0)
                node.setState(TaskNode.ProcessingState.VISITED)
                queue.offer(node)
                log.debug "[SLURM TASK GROUPING] Node ${node.getName()} has no dependencies, assigned level 0"
            }
        }
        
        // Process nodes level by level 
        int processedCount = 0
        while (!queue.isEmpty()) {
            final TaskNode current = queue.poll()
            processedCount++
            
            // Process all downstream dependencies
            for (Long downstreamId : current.getDownstreamDependencies()) {
                final TaskNode downstream = this.nodes.get(downstreamId)
                if (downstream != null) {
                    // Update incoming degree
                    final int newInDegree = inDegree.get(downstreamId) - 1
                    inDegree.put(downstreamId, newInDegree)
                    
                    // Calculate level based on maximum level of upstream dependencies
                    int maxUpstreamLevel = 0
                    for (Long upstreamId : downstream.getUpstreamDependencies()) {
                        final TaskNode upstream = this.nodes.get(upstreamId)
                        if (upstream != null && upstream.getState() == TaskNode.ProcessingState.VISITED) {
                            maxUpstreamLevel = Math.max(maxUpstreamLevel, upstream.getLevel())
                        }
                    }
                    downstream.setLevel(maxUpstreamLevel + 1)
                    
                    // If all dependencies are processed, add node to queue
                    if (newInDegree == 0) {
                        downstream.setState(TaskNode.ProcessingState.VISITED)
                        queue.offer(downstream)
                        log.debug "[SLURM TASK GROUPING] Node ${downstream.getName()} assigned level ${downstream.getLevel()}"
                    }
                }
            }
        }
        
        // Verify all nodes were processed (cycle detection)
        if (processedCount < this.nodes.size()) {
            log.error "[SLURM TASK GROUPING] Cycle detected: only ${processedCount} of ${this.nodes.size()} nodes processed"
            for (Map.Entry<Long, Integer> entry : inDegree.entrySet()) {
                if (entry.getValue() > 0) {
                    final TaskNode node = this.nodes.get(entry.getKey())
                    log.error "[SLURM TASK GROUPING] Node ${node?.getName()} (id=${entry.getKey()}) has ${entry.getValue()} unprocessed dependencies"
                }
            }
            throw new IllegalStateException("Cycle detected in dependency graph: ${this.nodes.size() - processedCount} nodes unreachable")
        }
        
        log.debug "[SLURM TASK GROUPING] Successfully assigned levels to ${processedCount} nodes"
    }
    
    /**
     * Group tasks by their assigned levels for parallel execution.
     * @return a map of level to list of TaskNode IDs at that level
     */
    Map<Integer, List<Long>> groupByLevels() {
        Map<Integer, List<Long>> tasksByLevels = new HashMap<>()
        for (TaskNode node : this.nodes.values()) {
            int level = node.getLevel()
            tasksByLevels.computeIfAbsent(level, k -> new ArrayList<>()).add(node.getId())
        }
        return tasksByLevels
    }
}