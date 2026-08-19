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
package nextflow.config.control;

import java.util.ArrayList;
import java.util.stream.Collectors;

import nextflow.config.ast.ConfigApplyNode;
import nextflow.config.ast.ConfigApplyBlockNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigVisitorSupport;
import nextflow.script.control.DirectiveRefCollector;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform a Nextflow config AST into a Groovy AST.
 *
 * @see nextflow.config.parser.v2.ConfigDsl
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ConfigToGroovyVisitor extends ConfigVisitorSupport {

    private SourceUnit sourceUnit;

    private ConfigNode moduleNode;

    /**
     * Whether to attach the referenced `task` directive names to the dynamic directive
     * values -- see {@link #transformDirectiveRefs}.
     *
     * It is disabled for the rendering path (`nextflow config`, `kuberun`), which replaces
     * every closure with its source text and therefore has no use for the names. Keeping the
     * two apart is what allows the wrapper to be introduced without teaching
     * `ClosureToStringVisitor` to unwrap it again: the paths are mutually exclusive, so the
     * value a renderer sees is always the plain closure the user wrote.
     */
    private boolean collectDirectiveRefs;

    public ConfigToGroovyVisitor(SourceUnit sourceUnit) {
        this(sourceUnit, true);
    }

    public ConfigToGroovyVisitor(SourceUnit sourceUnit, boolean collectDirectiveRefs) {
        this.sourceUnit = sourceUnit;
        this.moduleNode = (ConfigNode) sourceUnit.getAST();
        this.collectDirectiveRefs = collectDirectiveRefs;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        if( moduleNode == null )
            return;
        super.visit(moduleNode);
        if( moduleNode.isEmpty() )
            moduleNode.addStatement(ReturnStatement.RETURN_NULL_OR_VOID);
    }

    @Override
    public void visitConfigApplyBlock(ConfigApplyBlockNode node) {
        moduleNode.addStatement(transformConfigApplyBlock(node));
    }

    protected Statement transformConfigApplyBlock(ConfigApplyBlockNode node) {
        var statements = new ArrayList<Statement>();
        for( var call : node.statements )
            statements.add(stmt(call));
        var code = block(new VariableScope(), statements);
        return stmt(callThisX("block", args(constX(node.name), closureX(null, code))));
    }

    @Override
    public void visitConfigAssign(ConfigAssignNode node) {
        moduleNode.addStatement(transformConfigAssign(node));
    }

    protected Statement transformConfigAssign(ConfigAssignNode node) {
        var names = listX(
            node.names.stream()
                .map(name -> (Expression) constX(name))
                .collect(Collectors.toList())
        );
        return stmt(callThisX("assign", args(names, transformDirectiveRefs(node.value))));
    }

    /**
     * Attach to a dynamic directive value the names of the `task` directives that it
     * references, so that they are available at runtime *without* evaluating the value.
     *
     * Why the names cannot be recovered at runtime instead: a config value is compiled to
     * a closure, and a closure carries no source text, so by the time the value reaches the
     * executor there is nothing left to inspect. The names are therefore read off the AST
     * here -- where the source still exists -- and carried by the value itself.
     *
     * Why only a closure is wrapped: `task` is not defined in the config scope, so an
     * interpolated string referencing it (`ext.args = "-Xmx${task.memory}g"`) is rejected by
     * the parser with "`task` is not defined". A closure body is not resolved against the
     * config scope, which makes `ext.args = { "-Xmx${task.memory}g" }` the *only* way a
     * config value can reference a directive -- hence the only shape worth wrapping.
     *
     * @see nextflow.script.DirectiveRefsClosure
     * @see nextflow.script.control.DirectiveRefCollector
     * @param node
     */
    private Expression transformDirectiveRefs(Expression node) {
        if( !collectDirectiveRefs )
            return node;

        // a directive value may hold its closures nested in a list and/or map literal at any
        // depth, e.g. `publishDir = [[path: {..}], ..]` -- recurse so they are reached too.
        // Note this mirrors ClosureToStringVisitor#replaceClosures, which walks the very same
        // shapes to swap closures for their source text; the two must stay in sync
        if( node instanceof ListExpression le ) {
            var items = le.getExpressions();
            for( int i = 0; i < items.size(); i++ )
                items.set(i, transformDirectiveRefs(items.get(i)));
            return le;
        }
        if( node instanceof MapExpression me ) {
            for( var entry : me.getMapEntryExpressions() )
                entry.setValueExpression(transformDirectiveRefs(entry.getValueExpression()));
            return me;
        }

        // anything that is not a closure cannot reference `task` -- see the javadoc above
        if( !(node instanceof ClosureExpression) )
            return node;

        // leave the value untouched when it references no directive, so that the wrapper is
        // only ever introduced where it is actually needed
        var refs = DirectiveRefCollector.collect(node);
        if( refs.isEmpty() )
            return node;

        // sort the names to keep the generated bytecode reproducible, since the collector
        // returns an unordered set
        var names = refs.stream()
            .sorted()
            .map(name -> (Expression) constX(name))
            .collect(Collectors.toList());
        // referenced by name: this module cannot depend on the runtime class, which lives in
        // the `nextflow` module -- the same indirection used for `BodyDef` on the script side
        return createX("nextflow.script.DirectiveRefsClosure", node, listX(names));
    }

    @Override
    public void visitConfigBlock(ConfigBlockNode node) {
        moduleNode.addStatement(transformConfigBlock(node));
    }

    protected Statement transformConfigBlock(ConfigBlockNode node) {
        var statements = new ArrayList<Statement>();
        for( var stmt : node.statements ) {
            if( stmt instanceof ConfigAssignNode can )
                statements.add(transformConfigAssign(can));
            else if( stmt instanceof ConfigBlockNode cbn )
                statements.add(transformConfigBlock(cbn));
            else if( stmt instanceof ConfigIncludeNode cin )
                statements.add(transformConfigInclude(cin));
        }
        var code = block(new VariableScope(), statements);
        var kind = node.kind != null ? node.kind : "block";
        return stmt(callThisX(kind, args(constX(node.name), closureX(null, code))));
    }

    @Override
    public void visitConfigInclude(ConfigIncludeNode node) {
        moduleNode.addStatement(transformConfigInclude(node));
    }

    protected Statement transformConfigInclude(ConfigIncludeNode node) {
        return stmt(callThisX("includeConfig", args(node.source)));
    }

}
