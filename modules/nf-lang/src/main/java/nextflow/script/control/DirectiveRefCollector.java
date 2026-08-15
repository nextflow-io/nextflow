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

import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import nextflow.script.dsl.ProcessDsl;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.Statement;

/**
 * Collect the `task` directives referenced by an AST node, e.g. `memory` for an
 * expression interpolating `"-Xmx${task.memory.toGiga()}g"`.
 *
 * A rendered command carries the directive values that were *requested*, so an executor
 * that allows them to be adjusted at schedule time needs to know which ones the command
 * depends on. Collecting the references from the AST does not require any value to be
 * evaluated, therefore it also covers the directives that are resolved independently of
 * the task script e.g. `beforeScript`, `ext`.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class DirectiveRefCollector extends CodeVisitorSupport {

    /**
     * The process directive names, used to tell apart the `task` properties that are
     * directives (e.g. `task.memory`) from the ones that are not (e.g. `task.attempt`).
     */
    private static final Set<String> DIRECTIVE_NAMES = Arrays.stream(ProcessDsl.DirectiveDsl.class.getMethods())
        .map(Method::getName)
        .collect(Collectors.toSet());

    private final Set<String> names = new HashSet<>();

    public static Set<String> collect(Statement... nodes) {
        var collector = new DirectiveRefCollector();
        for( var node : nodes ) {
            if( node != null )
                node.visit(collector);
        }
        return collector.names;
    }

    public static Set<String> collect(Expression node) {
        var collector = new DirectiveRefCollector();
        if( node != null )
            node.visit(collector);
        return collector.names;
    }

    @Override
    public void visitPropertyExpression(PropertyExpression node) {
        var name = directiveRef(node);
        if( name != null )
            names.add(name);
        super.visitPropertyExpression(node);
    }

    /**
     * Return the directive name when the given expression is a `task.<name>` property
     * access, e.g. `memory` for the `task.memory` in `task.memory.toGiga()`, or `ext`
     * for the `task.ext` in `task.ext.args`.
     *
     * Note the `task` object also exposes task properties that are not directives e.g.
     * `task.attempt`, which are not reported.
     */
    private static String directiveRef(PropertyExpression node) {
        if( !(node.getObjectExpression() instanceof VariableExpression ve) )
            return null;
        if( !"task".equals(ve.getName()) )
            return null;
        var name = node.getPropertyAsString();
        return DIRECTIVE_NAMES.contains(name) ? name : null;
    }

}
