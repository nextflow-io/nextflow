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

import java.util.Comparator;
import java.util.List;

import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeEntryNode;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.OutputBlockNode;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ParamNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.parser.CommentWriter;

import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import groovy.lang.Tuple2;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Format a script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptFormattingVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private FormattingOptions options;

    private CommentWriter commentWriter;

    private Formatter fmt;

    private int maxIncludeWidth = 0;

    private int maxParamWidth = 0;

    public ScriptFormattingVisitor(SourceUnit sourceUnit, FormattingOptions options) {
        this.sourceUnit = sourceUnit;
        this.options = options;
        this.commentWriter = getCommentWriter();
        this.fmt = new Formatter(options, this.commentWriter);
    }

    public CommentWriter getCommentWriter() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ScriptNode ) {
            var scriptNode = (ScriptNode) moduleNode;
            return (CommentWriter)scriptNode.getNodeMetaData("commentWriter");
        }
        return null;
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

            maxParamWidth = scriptNode.getParams().stream()
                .map(this::getParamWidth)
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
            else if( decl instanceof ParamNode pn )
                visitParam(pn);
            else if( decl instanceof ProcessNode pn )
                visitProcess(pn);
            else if( decl instanceof WorkflowNode wn )
                visitWorkflow(wn);
        }
        fmt.appendBlanks(1);
        System.out.print(commentWriter.toNFComments());
    }

    public String toString() {
        return fmt.toString();
    }

    // script declarations

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
        fmt.appendLeadingComments(node);
        fmt.append(node.name);
        fmt.appendSpace();
        fmt.append("=");
        fmt.appendSpace();
        fmt.visit(node.value);
        fmt.appendNewLine();
    }

    @Override
    public void visitInclude(IncludeNode node) {
        var wrap = node.getLineNumber() < node.getLastLineNumber();
        fmt.appendLeadingComments(node);
        fmt.append("include");
        fmt.appendSpace();
        fmt.append('{');
        if( wrap )
            fmt.incIndent();
        for( int i = 0; i < node.entries.size(); i++ ) {
            if( wrap ) {
                fmt.appendNewLine();
                fmt.appendIndent();
            }
            else {
                fmt.appendSpace();
            }
            var entry = node.entries.get(i);
            fmt.append(entry.name);
            if( entry.alias != null ) {
                fmt.appendSpace();
                fmt.append("as");
                fmt.appendSpace();
                fmt.append(entry.alias);
            }
            if( !wrap && node.entries.size() == 1 && options.harshilAlignment() ) {
                var padding = maxIncludeWidth - getIncludeWidth(entry);
                fmt.appendPadding(padding);
            }
            if( i + 1 < node.entries.size() )
                fmt.appendSpace();
                fmt.append(";");
        }
        if( wrap ) {
            fmt.appendNewLine();
            fmt.decIndent();
        }
        else {
            fmt.appendSpace();
        }
        fmt.append('}');
        fmt.appendSpace();
        fmt.append("from");
        fmt.appendSpace();
        fmt.visit(node.source);
        fmt.appendNewLine();
    }

    protected int getIncludeWidth(IncludeEntryNode entry) {
        return entry.alias != null
            ? entry.name.length() + 4 + entry.alias.length()
            : entry.name.length();
    }

    @Override
    public void visitParam(ParamNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        fmt.visit(node.target);
        if( maxParamWidth > 0 ) {
            var padding = maxParamWidth - getParamWidth(node);
            fmt.appendPadding(padding);
        }
        fmt.appendSpace();
        fmt.append("=");
        fmt.appendSpace();
        fmt.visit(node.value);
        fmt.appendNewLine();
    }

    protected int getParamWidth(ParamNode node) {
        var target = (PropertyExpression) node.target;
        var name = target.getPropertyAsString();
        return name != null ? name.length() : 0;
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        fmt.appendLeadingComments(node);
        var withinComments = commentWriter.getWithinComments(node);
        fmt.appendNewLine();
        fmt.append("workflow");
        fmt.appendSpace();
        fmt.writeComments(withinComments, "KEYWORD", false, true, false, false, true);
        if( !node.isEntry() ) {
            fmt.append(node.getName());
            fmt.appendSpace();
            fmt.writeComments(withinComments, "NAME", false, true, false, false, true);
        }
        fmt.append("{");
        fmt.incIndent();
        if( node.takes instanceof BlockStatement ) {
            appendColonStatement("take", "TAKE", node.takes);
            visitWorkflowTakes(asBlockStatements(node.takes));
            fmt.appendNewLine();
        }
        if( node.main instanceof BlockStatement ) {
            if( node.takes instanceof BlockStatement || node.emits instanceof BlockStatement || node.publishers instanceof BlockStatement ) {
                fmt.appendIndent();
                appendColonStatement("main", "MAIN", node.main);
            }
            fmt.visit(node.main);
        }
        if( node.emits instanceof BlockStatement ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            appendColonStatement("emit", "EMIT", node.takes);
            visitWorkflowEmits(asBlockStatements(node.emits));
        }
        if( node.publishers instanceof BlockStatement ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            appendColonStatement("publish", "PUBLISH", node.publishers);
            fmt.visit(node.publishers);
        }
        fmt.decIndent();
        fmt.appendNewLine();
        fmt.append("}");
        fmt.appendTrailingComments(node);
    }

    protected void visitWorkflowTakes(List<Statement> takes) {
        var alignmentWidth = options.harshilAlignment()
            ? getMaxParameterWidth(takes)
            : 0;

        for( var stmt : takes ) {
            var ve = asVarX(stmt);
            fmt.appendIndent();
            fmt.visit(ve);
            // if( fmt.hasTrailingComment(stmt) ) {
            //     if( alignmentWidth > 0 ) {
            //         var padding = alignmentWidth - ve.getName().length();
            //         fmt.append(" ".repeat(padding));
            //     }
            // }
            fmt.appendTrailingComments(stmt);
            fmt.appendNewLine();
        }
    }

    protected void visitWorkflowEmits(List<Statement> emits) {
        var alignmentWidth = options.harshilAlignment()
            ? getMaxParameterWidth(emits)
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
                    fmt.appendPadding(padding);
                }
                fmt.appendSpace();
                fmt.append("=");
                fmt.appendSpace();
                fmt.visit(assign.getRightExpression());
                fmt.appendTrailingComments(stmt);
                fmt.appendNewLine();
            }
            else if( emit instanceof VariableExpression ve ) {
                fmt.appendIndent();
                fmt.visit(ve);
                // if( fmt.hasTrailingComment(stmt) ) {
                //     if( alignmentWidth > 0 ) {
                //         var padding = alignmentWidth - ve.getName().length();
                //         fmt.append(" ".repeat(padding));
                //     }
                // }
                fmt.appendTrailingComments(stmt);
                fmt.appendNewLine();
            }
            else {
                fmt.visit(stmt);
            }
        }
    }

    protected int getMaxParameterWidth(List<Statement> statements) {
        if( statements.size() == 1 )
            return 0;

        int maxWidth = 0;
        for( var stmt : statements ) {
            var stmtX = (ExpressionStatement)stmt;
            var emit = stmtX.getExpression();
            int width = 0;
            if( emit instanceof VariableExpression ve ) {
                width = ve.getName().length();
            }
            else if( emit instanceof AssignmentExpression assign ) {
                var target = (VariableExpression)assign.getLeftExpression();
                width = target.getName().length();
            }

            if( maxWidth < width )
                maxWidth = width;
        }
        return maxWidth;
    }

    @Override
    public void visitProcess(ProcessNode node) {
        var withinComments = commentWriter.getWithinComments(node);
        var trailingComments = commentWriter.getTrailingComments(node);
        System.err.println(trailingComments);
        fmt.appendLeadingComments(node);
        fmt.appendNewLine();
        fmt.append("process");
        fmt.appendSpace();
        fmt.writeComments(withinComments, "KEYWORD", false, false, true, false, true);
        fmt.append(node.getName());
        fmt.writeComments(withinComments, "NAME", false, true, false, false, true);
        fmt.appendSpace();
        fmt.append("{");
        fmt.writeComments(trailingComments, "LBRACE", false, true, true, false, false);
        fmt.incIndent();

        if( node.directives instanceof BlockStatement ) {
            visitDirectives(node.directives);
            fmt.appendBlanks(1);
        }

        if( node.inputs instanceof BlockStatement ) {
            appendColonStatement("input", "INPUT", node.inputs);
            visitDirectives(node.inputs);
            fmt.appendBlanks(1);
        }

        if( !options.maheshForm() && node.outputs instanceof BlockStatement ) {
            appendColonStatement("output", "OUTPUT", node.outputs);
            visitDirectives(node.outputs);
            fmt.appendBlanks(1);
        }

        if( !(node.when instanceof EmptyExpression) ) {
            appendColonStatement("when", "WHEN", node.when);
            fmt.incIndent();
            fmt.appendNewLine();
            fmt.visit(node.when);
            fmt.appendBlanks(1);
        }

        appendColonStatement(node.type, "EXEC", node.exec);
        fmt.visit(node.exec);
        fmt.appendBlanks(1);

        if( !(node.stub instanceof EmptyStatement) ) {
            appendColonStatement("stub", "STUB", node.stub);
            fmt.visit(node.stub);
            fmt.appendBlanks(1);
        }

        if( options.maheshForm() && node.outputs instanceof BlockStatement ) {
            appendColonStatement("output", "OUTPUT", node.outputs);
            visitDirectives(node.outputs);
            fmt.appendBlanks(1);
        }

        fmt.decIndent();
        fmt.appendNewLine();
        fmt.append("}");
        fmt.writeComments(trailingComments, "RBRACE", false, true, true, false, false);
    }

    private void appendColonStatement(String name, String commentKey, ASTNode stmt) {
        fmt.appendLeadingComments(stmt);
        fmt.appendNewLine();
        var colonWithinComments = commentWriter.getWithinComments(stmt);
        var colonTrailingComments = commentWriter.getTrailingComments(stmt);
        // Write the name of the statement
        fmt.append(name);
        // Write commend inbetween the name and the colon
        fmt.appendWithinComments(colonWithinComments, commentKey);
        fmt.append(":");
        // Write trailing comments after the colon
        fmt.appendTrailingInside(colonTrailingComments, "COLON");
    }

    @Override
    public void visitFunction(FunctionNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendNewLine();
        fmt.append("def");
        fmt.appendSpace();
        if( Formatter.isLegacyType(node.getReturnType()) ) {
            fmt.visitTypeAnnotation(node.getReturnType());
            fmt.appendSpace();
        }
        fmt.append(node.getName());
        fmt.append('(');
        fmt.visitParameters(node.getParameters());
        fmt.append(')');
        fmt.appendSpace();
        fmt.append('{');
        fmt.incIndent();
        fmt.visit(node.getCode());
        fmt.decIndent();
        fmt.appendNewLine();
        fmt.append("}");
    }

    @Override
    public void visitEnum(ClassNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendNewLine();
        fmt.append("enum");
        fmt.appendSpace();
        fmt.append(node.getName());
        fmt.incIndent();
        fmt.appendSpace();
        fmt.append("{");
        for( var fn : node.getFields() ) {
            fmt.appendNewLine();
            fmt.append(fn.getName());
            fmt.append(',');
        }
        fmt.decIndent();
        fmt.appendNewLine();
        fmt.append("}");
    }

    @Override
    public void visitOutputs(OutputBlockNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendNewLine();
        fmt.append("output");
        fmt.appendSpace();
        fmt.append('{');
        fmt.incIndent();
        super.visitOutputs(node);
        fmt.decIndent();
        fmt.appendNewLine();
        fmt.append("}");
    }

    @Override
    public void visitOutput(OutputNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendNewLine();
        fmt.append(node.name);
        fmt.appendSpace();
        fmt.append("{");
        visitOutputBody((BlockStatement) node.body);
        fmt.decIndent();
        fmt.appendNewLine();
        fmt.append("}");
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
                    fmt.appendSpace();
                    fmt.append("{");
                    fmt.incIndent();
                    fmt.appendNewLine();
                    visitDirectives(code);
                    fmt.decIndent();
                    fmt.appendIndent();
                    fmt.append("}");
                    fmt.appendNewLine();
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
