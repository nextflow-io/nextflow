/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.script.control;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Resolve all fully-qualified process names invoked
 * (directly or indirectly) by an entry workflow.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ProcessNameResolver {

    private Map<WorkflowNode, Map<String, MethodNode>> callSites;

    public ProcessNameResolver(Map<WorkflowNode, Map<String, MethodNode>> callSites) {
        this.callSites = callSites;
    }

    public Set<String> resolve(SourceUnit main) {
        var result = new HashSet<String>();
        var queue = new LinkedList<CallScope>();

        if( main.getAST() instanceof ScriptNode sn && sn.getEntry() != null )
            queue.add(new CallScope("", sn.getEntry()));

        while( !queue.isEmpty() ) {
            var scope = queue.remove();
            var calls = callSites.get(scope.node());
            calls.forEach((name, mn) -> {
                if( mn instanceof WorkflowNode wn ) {
                    var workflowName = fullyQualifiedName(scope.name(), name);
                    queue.add(new CallScope(workflowName, wn));
                }
                else if( mn instanceof ProcessNode pn ) {
                    var processName = fullyQualifiedName(scope.name(), name);
                    result.add(processName);
                }
            });
        }

        return result;
    }

    private String fullyQualifiedName(String scope, String name) {
        return scope.isEmpty()
            ? name
            : scope + ":" + name;
    }

    private static record CallScope(
        String name,
        WorkflowNode node
    ) {}
}
