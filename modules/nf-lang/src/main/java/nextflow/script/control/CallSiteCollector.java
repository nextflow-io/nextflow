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

import java.util.IdentityHashMap;
import java.util.HashMap;
import java.util.Map;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Collect call sites for each workflow in a script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class CallSiteCollector {

    public Map<WorkflowNode, Map<String, MethodNode>> apply(SourceUnit source) {
        var callSites = new IdentityHashMap<WorkflowNode, Map<String, MethodNode>>();
        if( source.getAST() instanceof ScriptNode sn ) {
            for ( var wn : sn.getWorkflows() )
                callSites.put(wn, new WorkflowVisitor().apply(wn));
        }
        return callSites;
    }

    public class WorkflowVisitor extends CodeVisitorSupport {

        private Map<String, MethodNode> calls;

        public Map<String, MethodNode> apply(WorkflowNode node) {
            calls = new HashMap<>();
            visit(node.main);
            visit(node.emits);
            visit(node.publishers);
            return calls;
        }

        @Override
        public void visitMethodCallExpression(MethodCallExpression node) {
            visit(node.getObjectExpression());
            visit(node.getArguments());

            if( node.isImplicitThis() ) {
                var mn = (MethodNode) node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
                if( mn instanceof WorkflowNode || mn instanceof ProcessNode )
                    calls.put(node.getMethodAsString(), mn);
            }
        }
    }
}
