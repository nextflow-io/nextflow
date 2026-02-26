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
        return nodes
    }

    DAG getDag() {
        return dag
    }

    /**
     * Add a TaskNode to the graph
     * @param node the TaskNode to add
     */
    void addNode(TaskNode node) {
        nodes.put(node.getId(), node)
    }

    /**
     * Build dependencies between nodes based on DAG edges.
     * Only considers edges between PROCESS vertices.
     */
    void buildDependencies() {
        for (DAG.Edge edge : dag.edges) {
            // Only process edges between actual process vertices
            if (edge.from?.type == DAG.Type.PROCESS && edge.to?.type == DAG.Type.PROCESS) {
                final TaskNode fromNode = nodes.get(edge.from.id)
                final TaskNode toNode = nodes.get(edge.to.id)
                
                if (fromNode != null && toNode != null) {
                    toNode.addUpstreamDependency(edge.from.id)
                    fromNode.addDownstreamDependency(edge.to.id)
                }
            }
        }
    }
}