/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.formatter;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeEntryNode;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.OutputBlockNode;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ParamNodeV1;
import nextflow.script.ast.ParamBlockNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Format a script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptFormattingVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private FormattingOptions options;

    private Formatter fmt;

    private int maxIncludeWidth = 0;

    private int maxParamWidth = 0;

    public ScriptFormattingVisitor(SourceUnit sourceUnit, FormattingOptions options) {
        this.sourceUnit = sourceUnit;
        this.options = options;
        this.fmt = new Formatter(options);
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( !(moduleNode instanceof ScriptNode) )
            return;
        var scriptNode = (ScriptNode) moduleNode;
        if( scriptNode.getShebang() != null )
            fmt.append(scriptNode.getShebang());

        // format code snippet if applicable
        var entry = scriptNode.getEntry();
        if( entry != null && entry.isCodeSnippet() ) {
            fmt.visit(entry.main);
            return;
        }

        // format script declarations
        var declarations = scriptNode.getDeclarations();

        // -- declarations are canonically sorted by default
        // -- revert to original order if sorting is disabled
        if( !options.sortDeclarations() ) {
            declarations.sort(Comparator.comparing(node -> node.getLineNumber()));
        }

        // -- prepare alignment widths if needed
        if( options.harshilAlignment() ) {
            maxIncludeWidth = scriptNode.getIncludes().stream()
                .flatMap(in -> in.entries.stream())
                .map(this::getIncludeWidth)
                .max(Integer::compare).orElse(0);

            maxParamWidth = scriptNode.getParamsV1().stream()
                .map(ScriptFormattingVisitor::parameterWidth)
                .max(Integer::compare).orElse(0);
        }

        for( var decl : declarations ) {
            if( decl instanceof ClassNode cn && cn.isEnum() )
                visitEnum(cn);
            else if( decl instanceof FeatureFlagNode ffn )
                visitFeatureFlag(ffn);
            else if( decl instanceof FunctionNode fn )
                visitFunction(fn);
            else if( decl instanceof IncludeNode in )
                visitInclude(in);
            else if( decl instanceof OutputBlockNode obn )
                visitOutputs(obn);
            else if( decl instanceof ParamBlockNode pbn )
                visitParams(pbn);
            else if( decl instanceof ParamNodeV1 pn )
                visitParamV1(pn);
            else if( decl instanceof ProcessNode pn )
                visitProcess(pn);
            else if( decl instanceof WorkflowNode wn )
                visitWorkflow(wn);
        }
    }

    public String toString() {
        return fmt.toString();
    }

    // script declarations

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
        fmt.appendLeadingComments(node);
        fmt.append(node.name);
        fmt.append(" = ");
        fmt.visit(node.value);
        fmt.appendNewLine();
    }

    @Override
    public void visitInclude(IncludeNode node) {
        var wrap = node.getLineNumber() < node.getLastLineNumber();
        fmt.appendLeadingComments(node);
        fmt.append("include {");
        if( wrap )
            fmt.incIndent();
        for( int i = 0; i < node.entries.size(); i++ ) {
            if( wrap ) {
                fmt.appendNewLine();
                fmt.appendIndent();
            }
            else {
                fmt.append(' ');
            }
            var entry = node.entries.get(i);
            fmt.append(entry.name);
            if( entry.alias != null ) {
                fmt.append(" as ");
                fmt.append(entry.alias);
            }
            if( !wrap && node.entries.size() == 1 && options.harshilAlignment() ) {
                var padding = maxIncludeWidth - getIncludeWidth(entry);
                fmt.append(" ".repeat(padding));
            }
            if( i + 1 < node.entries.size() )
                fmt.append(" ;");
        }
        if( wrap ) {
            fmt.appendNewLine();
            fmt.decIndent();
        }
        else {
            fmt.append(' ');
        }
        fmt.append("} from ");
        fmt.visit(node.source);
        fmt.appendNewLine();
    }

    protected int getIncludeWidth(IncludeEntryNode entry) {
        return entry.alias != null
            ? entry.name.length() + 4 + entry.alias.length()
            : entry.name.length();
    }

    @Override
    public void visitParams(ParamBlockNode node) {
        var alignmentWidth = options.harshilAlignment()
            ? maxParameterWidth(node.declarations)
            : 0;

        fmt.appendLeadingComments(node);
        fmt.append("params {\n");
        fmt.incIndent();
        for( var param : node.declarations ) {
            fmt.appendLeadingComments(param);
            fmt.appendIndent();
            fmt.append(param.getName());
            if( fmt.hasType(param) ) {
                if( alignmentWidth > 0 ) {
                    var padding = alignmentWidth - param.getName().length() + 1;
                    fmt.append(" ".repeat(padding));
                }
                fmt.append(": ");
                fmt.visitTypeAnnotation(param.getType());
            }
            if( param.hasInitialExpression() ) {
                fmt.append(" = ");
                fmt.visit(param.getInitialExpression());
            }
            fmt.appendNewLine();
        }
        fmt.decIndent();
        fmt.append("}\n");
    }

    private static int maxParameterWidth(Parameter[] parameters) {
        return Arrays.stream(parameters)
            .map(param -> param.getName().length())
            .max(Integer::compare).orElse(0);
    }

    @Override
    public void visitParamV1(ParamNodeV1 node) {
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        fmt.visit(node.target);
        if( maxParamWidth > 0 ) {
            var padding = maxParamWidth - parameterWidth(node);
            fmt.append(" ".repeat(padding));
        }
        fmt.append(" = ");
        fmt.visit(node.value);
        fmt.appendNewLine();
    }

    private static int parameterWidth(ParamNodeV1 node) {
        var target = (PropertyExpression) node.target;
        var name = target.getPropertyAsString();
        return name != null ? name.length() : 0;
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        fmt.appendLeadingComments(node);
        fmt.append("workflow");
        if( !node.isEntry() ) {
            fmt.append(' ');
            fmt.append(node.getName());
        }
        fmt.append(" {\n");
        fmt.incIndent();
        if( node.takes instanceof BlockStatement ) {
            fmt.appendIndent();
            fmt.append("take:\n");
            visitWorkflowTakes(asBlockStatements(node.takes));
        }
        if( node.main instanceof BlockStatement ) {
            if( node.takes instanceof BlockStatement ) {
                fmt.appendNewLine();
            }
            if( node.takes instanceof BlockStatement || node.emits instanceof BlockStatement || node.publishers instanceof BlockStatement ) {
                fmt.appendIndent();
                fmt.append("main:\n");
            }
            fmt.visit(node.main);
        }
        if( node.emits instanceof BlockStatement ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("emit:\n");
            visitWorkflowEmits(asBlockStatements(node.emits));
        }
        if( node.publishers instanceof BlockStatement ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("publish:\n");
            fmt.visit(node.publishers);
        }
        if( node.onComplete instanceof BlockStatement ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("onComplete:\n");
            fmt.visit(node.onComplete);
        }
        if( node.onError instanceof BlockStatement ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("onError:\n");
            fmt.visit(node.onError);
        }
        fmt.decIndent();
        fmt.append("}\n");
    }

    protected void visitWorkflowTakes(List<Statement> takes) {
        var alignmentWidth = options.harshilAlignment()
            ? maxParameterWidth(takes)
            : 0;

        for( var stmt : takes ) {
            var ve = asVarX(stmt);
            fmt.appendIndent();
            fmt.visit(ve);
            if( fmt.hasTrailingComment(stmt) ) {
                if( alignmentWidth > 0 ) {
                    var padding = alignmentWidth - ve.getName().length();
                    fmt.append(" ".repeat(padding));
                }
                fmt.appendTrailingComment(stmt);
            }
            fmt.appendNewLine();
        }
    }

    protected void visitWorkflowEmits(List<Statement> emits) {
        var alignmentWidth = options.harshilAlignment()
            ? maxParameterWidth(emits)
            : 0;

        for( var stmt : emits ) {
            var stmtX = (ExpressionStatement)stmt;
            var emit = stmtX.getExpression();
            if( emit instanceof AssignmentExpression assign ) {
                var ve = (VariableExpression)assign.getLeftExpression();
                fmt.appendIndent();
                fmt.visit(ve);
                if( alignmentWidth > 0 ) {
                    var padding = alignmentWidth - ve.getName().length();
                    fmt.append(" ".repeat(padding));
                }
                fmt.append(" = ");
                fmt.visit(assign.getRightExpression());
                fmt.appendTrailingComment(stmt);
                fmt.appendNewLine();
            }
            else if( emit instanceof VariableExpression ve ) {
                fmt.appendIndent();
                fmt.visit(ve);
                if( fmt.hasTrailingComment(stmt) ) {
                    if( alignmentWidth > 0 ) {
                        var padding = alignmentWidth - ve.getName().length();
                        fmt.append(" ".repeat(padding));
                    }
                    fmt.appendTrailingComment(stmt);
                }
                fmt.appendNewLine();
            }
            else {
                fmt.visit(stmt);
            }
        }
    }

    private static int maxParameterWidth(List<Statement> statements) {
        if( statements.size() == 1 )
            return 0;

        return statements.stream()
            .map((stmt) -> {
                var stmtX = (ExpressionStatement)stmt;
                var emit = stmtX.getExpression();
                if( emit instanceof VariableExpression ve ) {
                    return ve.getName().length();
                }
                if( emit instanceof AssignmentExpression assign ) {
                    var target = (VariableExpression)assign.getLeftExpression();
                    return target.getName().length();
                }
                return 0;
            })
            .max(Integer::compare).orElse(0);
    }

    @Override
    public void visitProcess(ProcessNode node) {
        fmt.appendLeadingComments(node);
        fmt.append("process ");
        fmt.append(node.getName());
        fmt.append(" {\n");
        fmt.incIndent();
        if( node.directives instanceof BlockStatement ) {
            visitDirectives(node.directives);
            fmt.appendNewLine();
        }
        if( node.inputs instanceof BlockStatement ) {
            fmt.appendIndent();
            fmt.append("input:\n");
            visitDirectives(node.inputs);
            fmt.appendNewLine();
        }
        if( !options.maheshForm() && node.outputs instanceof BlockStatement ) {
            visitProcessOutputs(node.outputs);
            fmt.appendNewLine();
        }
        if( !(node.when instanceof EmptyExpression) ) {
            fmt.appendIndent();
            fmt.append("when:\n");
            fmt.appendIndent();
            fmt.visit(node.when);
            fmt.append("\n\n");
        }
        fmt.appendIndent();
        fmt.append(node.type);
        fmt.append(":\n");
        fmt.visit(node.exec);
        if( !(node.stub instanceof EmptyStatement) ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("stub:\n");
            fmt.visit(node.stub);
        }
        if( options.maheshForm() && node.outputs instanceof BlockStatement ) {
            fmt.appendNewLine();
            visitProcessOutputs(node.outputs);
        }
        fmt.decIndent();
        fmt.append("}\n");
    }

    private void visitProcessOutputs(Statement outputs) {
        fmt.appendIndent();
        fmt.append("output:\n");
        visitDirectives(outputs);
    }

    @Override
    public void visitFunction(FunctionNode node) {
        fmt.appendLeadingComments(node);
        fmt.append("def ");
        if( Formatter.isLegacyType(node.getReturnType()) ) {
            fmt.visitTypeAnnotation(node.getReturnType());
            fmt.append(' ');
        }
        fmt.append(node.getName());
        fmt.append('(');
        fmt.visitParameters(node.getParameters());
        fmt.append(") {\n");
        fmt.incIndent();
        fmt.visit(node.getCode());
        fmt.decIndent();
        fmt.append("}\n");
    }

    @Override
    public void visitEnum(ClassNode node) {
        fmt.appendLeadingComments(node);
        fmt.append("enum ");
        fmt.append(node.getName());
        fmt.append(" {\n");
        fmt.incIndent();
        for( var fn : node.getFields() ) {
            fmt.appendIndent();
            fmt.append(fn.getName());
            fmt.append(',');
            fmt.appendNewLine();
        }
        fmt.decIndent();
        fmt.append("}\n");
    }

    @Override
    public void visitOutputs(OutputBlockNode node) {
        fmt.appendLeadingComments(node);
        fmt.append("output {\n");
        fmt.incIndent();
        super.visitOutputs(node);
        fmt.decIndent();
        fmt.append("}\n");
    }

    @Override
    public void visitOutput(OutputNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        fmt.append(node.name);
        fmt.append(" {\n");
        fmt.incIndent();
        visitOutputBody((BlockStatement) node.body);
        fmt.decIndent();
        fmt.appendIndent();
        fmt.append("}\n");
    }

    protected void visitOutputBody(BlockStatement block) {
        asBlockStatements(block).forEach((stmt) -> {
            var call = asMethodCallX(stmt);
            if( call == null )
                return;

            // treat as index definition
            var name = call.getMethodAsString();
            if( "index".equals(name) ) {
                var code = asDslBlock(call, 1);
                if( code != null ) {
                    fmt.appendLeadingComments(stmt);
                    fmt.appendIndent();
                    fmt.append(name);
                    fmt.append(" {\n");
                    fmt.incIndent();
                    visitDirectives(code);
                    fmt.decIndent();
                    fmt.appendIndent();
                    fmt.append("}\n");
                    return;
                }
            }

            // treat as regular directive
            fmt.appendLeadingComments(stmt);
            fmt.visitDirective(call);
        });
    }

    protected void visitDirectives(Statement statement) {
        asBlockStatements(statement).forEach((stmt) -> {
            var call = asMethodCallX(stmt);
            if( call == null )
                return;
            fmt.appendLeadingComments(stmt);
            fmt.visitDirective(call);
        });
    }

}
