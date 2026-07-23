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
package nextflow.script.formatter;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import nextflow.script.ast.ASTNodeMarker;
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
import nextflow.script.ast.ProcessNodeV1;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.TupleParameter;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Format a script.
 *
 * All comments are preserved: comment metadata is re-derived from the
 * source text with {@link CommentReattacher} and emitted as leading,
 * trailing and dangling comments through {@link Formatter}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptFormattingVisitor extends ScriptVisitorSupport {

    private SourceUnit sourceUnit;

    private FormattingOptions options;

    private String sourceText;

    private Formatter fmt;

    private int maxIncludeWidth = 0;

    private int maxParamWidth = 0;

    public ScriptFormattingVisitor(SourceUnit sourceUnit, FormattingOptions options, String sourceText) {
        this.sourceUnit = sourceUnit;
        this.options = options;
        this.sourceText = sourceText;
        this.fmt = new Formatter(options);
    }

    public ScriptFormattingVisitor(SourceUnit sourceUnit, FormattingOptions options) {
        this(sourceUnit, options, null);
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

        // re-derive comment metadata from the source so that no comment is
        // lost; if the source is not available, fall back to the comment
        // metadata attached by the parser
        if( sourceText == null )
            sourceText = CommentReattacher.readSourceText(sourceUnit);
        var reattached = sourceText != null;
        if( reattached )
            CommentReattacher.apply(scriptNode, sourceText);

        if( scriptNode.getShebang() != null ) {
            fmt.append(scriptNode.getShebang());
            // the line terminator is part of the comment metadata unless
            // the comment metadata was re-derived
            if( reattached )
                fmt.appendNewLine();
        }

        // format code snippet if applicable
        var entry = scriptNode.getEntry();
        if( entry != null && entry.isCodeSnippet() ) {
            fmt.visit(entry.main);
            fmt.appendDanglingComments(scriptNode);
            return;
        }

        // format script declarations
        var declarations = scriptNode.getDeclarations();

        // -- declarations are canonically sorted by default
        // -- revert to original order if sorting is disabled
        if( !options.sortDeclarations() ) {
            declarations.sort(Comparator.comparing(node -> node.getLineNumber()));
        }
        else {
            sortIncludes(declarations);
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

        ASTNode prevDecl = null;
        for( var decl : declarations ) {
            // a declaration that is part of a verbatim region emitted by an
            // earlier declaration emits nothing -- do not emit a blank line
            // for it
            if( fmt.isVerbatimSuppressed(decl) )
                continue;
            // a declaration that is sorted into first position may carry a
            // blank-line prefix from its original location -- strip it so
            // that the output does not start with blank lines
            if( prevDecl == null )
                fmt.stripLeadingBlankPrefix(decl);
            // enforce a blank line above blocks (process, workflow, ...)
            // and between declarations of different kinds; consecutive
            // declarations of the same kind (e.g. includes) may be grouped
            if( needsBlankLineBetween(prevDecl, decl) && !fmt.leadingStartsWithBlankLine(decl) )
                fmt.appendNewLine();
            prevDecl = decl;

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
            else if( decl instanceof RecordNode rn )
                visitRecord(rn);
            else if( decl instanceof WorkflowNode wn )
                visitWorkflow(wn);
        }

        // emit comments at the end of the file
        fmt.appendDanglingComments(scriptNode);
    }

    /**
     * Determine whether a blank line should separate two consecutive
     * declarations: block declarations always get a blank line above, and
     * simple declarations (includes, feature flags, legacy params) get one
     * when the kind changes. The first declaration gets a blank line only
     * after a shebang.
     */
    private boolean needsBlankLineBetween(ASTNode prevDecl, ASTNode decl) {
        if( prevDecl == null ) {
            var scriptNode = (ScriptNode) sourceUnit.getAST();
            return scriptNode.getShebang() != null;
        }
        var prevKind = declKind(prevDecl);
        var kind = declKind(decl);
        return prevKind == DeclKind.BLOCK || kind == DeclKind.BLOCK || prevKind != kind;
    }

    private enum DeclKind {
        FEATURE_FLAG,
        INCLUDE,
        PARAM_V1,
        BLOCK
    }

    private static DeclKind declKind(ASTNode decl) {
        if( decl instanceof FeatureFlagNode )
            return DeclKind.FEATURE_FLAG;
        if( decl instanceof IncludeNode )
            return DeclKind.INCLUDE;
        if( decl instanceof ParamNodeV1 )
            return DeclKind.PARAM_V1;
        return DeclKind.BLOCK;
    }

    /**
     * Sort includes alphabetically by module path within their blank-line
     * separated groups (see nextflow-io/language-server#54). Comments above
     * the first include of a group are treated as the group header and
     * stay at the top of the group; comments above other includes move
     * with their include.
     */
    private void sortIncludes(List<ASTNode> declarations) {
        int i = 0;
        while( i < declarations.size() ) {
            if( !(declarations.get(i) instanceof IncludeNode) ) {
                i++;
                continue;
            }
            // find a contiguous run of includes and split it into groups
            // at blank lines
            int start = i;
            while( i < declarations.size() && declarations.get(i) instanceof IncludeNode )
                i++;
            int groupStart = start;
            for( int j = start + 1; j <= i; j++ ) {
                if( j == i || fmt.leadingStartsWithBlankLine(declarations.get(j)) ) {
                    sortIncludeGroup(declarations, groupStart, j);
                    groupStart = j;
                }
            }
        }
    }

    private void sortIncludeGroup(List<ASTNode> declarations, int start, int end) {
        if( end - start < 2 )
            return;
        var group = new ArrayList<IncludeNode>();
        for( int i = start; i < end; i++ )
            group.add((IncludeNode) declarations.get(i));

        var originalFirst = group.get(0);
        group.sort(
            Comparator.comparing((IncludeNode in) -> in.source.getText(), String.CASE_INSENSITIVE_ORDER)
                .thenComparing(in -> in.source.getText())
        );

        // keep the group header (the leading comments and blank line of the
        // original first include) at the top of the group
        var newFirst = group.get(0);
        if( newFirst != originalFirst ) {
            var marker = ASTNodeMarker.LEADING_COMMENTS;
            var header = (List<String>) originalFirst.getNodeMetaData(marker);
            if( header != null ) {
                originalFirst.removeNodeMetaData(marker);
                var own = (List<String>) newFirst.getNodeMetaData(marker);
                var combined = own != null ? new ArrayList<String>(own) : new ArrayList<String>();
                // lists are in reverse source order -- the header comes first
                combined.addAll(header);
                newFirst.putNodeMetaData(marker, combined);
            }
        }

        for( int i = start; i < end; i++ )
            declarations.set(i, group.get(i - start));
    }

    public String toString() {
        return fmt.toString();
    }

    // script declarations

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append(node.name);
        fmt.append(" = ");
        fmt.visit(node.value);
        fmt.appendTrailingComment(node);
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
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    private int getIncludeWidth(IncludeEntryNode entry) {
        return entry.alias != null
            ? entry.name.length() + 4 + entry.alias.length()
            : entry.name.length();
    }

    @Override
    public void visitParams(ParamBlockNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("params {\n");
        fmt.incIndent();
        for( var param : node.declarations ) {
            if( fmt.appendVerbatim(param) )
                continue;
            fmt.appendLeadingComments(param);
            fmt.emitWrappable(() -> {
                fmt.appendIndent();
                fmt.append(param.getName());
                if( fmt.hasType(param) ) {
                    fmt.append(": ");
                    fmt.visitTypeAnnotation(param.getType());
                }
                if( param.hasInitialExpression() ) {
                    fmt.append(" = ");
                    fmt.visitRootExpression(param.getInitialExpression());
                }
                fmt.appendTrailingComment(param);
                fmt.appendNewLine();
            });
            fmt.appendDanglingAfter(param);
        }
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    @Override
    public void visitParamV1(ParamNodeV1 node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.emitWrappable(() -> {
            fmt.appendIndent();
            fmt.visit(node.target);
            if( maxParamWidth > 0 ) {
                var padding = maxParamWidth - parameterWidth(node);
                fmt.append(" ".repeat(padding));
            }
            fmt.append(" = ");
            fmt.visitRootExpression(node.value);
            fmt.appendTrailingComment(node);
            fmt.appendNewLine();
        });
    }

    private static int parameterWidth(ParamNodeV1 node) {
        var target = (PropertyExpression) node.target;
        var name = target.getPropertyAsString();
        return name != null ? name.length() : 0;
    }

    @Override
    public void visitWorkflow(WorkflowNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("workflow");
        if( !node.isEntry() ) {
            fmt.append(' ');
            fmt.append(node.getName());
        }
        fmt.append(" {\n");
        fmt.incIndent();
        var takes = node.getParameters();
        if( takes.length > 0 ) {
            fmt.appendIndent();
            fmt.append("take:\n");
            visitTypedInputs(takes);
        }
        if( !node.main.isEmpty() ) {
            // NOTE: the main: label is also required when there is an
            // onComplete/onError handler
            if( takes.length > 0 || !node.emits.isEmpty() || !node.publishers.isEmpty() || !node.onComplete.isEmpty() || !node.onError.isEmpty() ) {
                // separate the main: section from the take: section above
                // it, without emitting a blank line at the start of the body
                if( takes.length > 0 )
                    fmt.appendNewLine();
                fmt.appendIndent();
                fmt.append("main:\n");
            }
            fmt.visit(node.main);
        }
        if( !node.emits.isEmpty() ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("emit:\n");
            visitTypedOutputs(asBlockStatements(node.emits));
            fmt.appendDanglingComments(node.emits);
        }
        if( !node.publishers.isEmpty() ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("publish:\n");
            visitWorkflowPublishers(asBlockStatements(node.publishers));
            fmt.appendDanglingComments(node.publishers);
        }
        if( !node.onComplete.isEmpty() ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("onComplete:\n");
            fmt.visit(node.onComplete);
        }
        if( !node.onError.isEmpty() ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("onError:\n");
            fmt.visit(node.onError);
        }
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    private void visitTypedInputs(Parameter[] inputs) {
        for( var input : inputs ) {
            if( fmt.appendVerbatim(input) ) {
                fmt.appendDanglingAfter(input);
                continue;
            }
            fmt.appendLeadingComments(input);
            fmt.appendIndent();
            if( input instanceof TupleParameter tp ) {
                visitStructuredInput(tp);
            }
            else {
                visitTypedInput(input);
            }
            fmt.appendTrailingComment(input);
            fmt.appendNewLine();
            fmt.appendDanglingAfter(input);
        }
    }

    private void visitStructuredInput(TupleParameter tp) {
        var isRecord = "Record".equals(tp.getType().getNameWithoutPackage());
        var wrap = isRecord;

        fmt.append(isRecord ? "record" : "tuple");
        fmt.append('(');
        if( wrap )
            fmt.incIndent();
        for( int i = 0; i < tp.components.length; i++ ) {
            var p = tp.components[i];
            if( wrap ) {
                fmt.appendNewLine();
                fmt.appendIndent();
            }
            fmt.append(p.getName());
            if( fmt.hasType(p) ) {
                fmt.append(": ");
                fmt.visitTypeAnnotation(p.getType());
            }
            if( i < tp.components.length - 1 )
                fmt.append(wrap ? "," : ", ");
        }
        if( wrap ) {
            fmt.appendNewLine();
            fmt.decIndent();
            fmt.appendIndent();
        }
        fmt.append(')');
    }

    private void visitTypedInput(Parameter node) {
        fmt.append(node.getName());
        if( fmt.hasType(node) ) {
            fmt.append(": ");
            fmt.visitTypeAnnotation(node.getType());
        }
    }

    private void visitTypedOutputs(List<Statement> outputs) {
        var alignmentWidth = options.harshilAlignment()
            ? maxParameterWidth(outputs)
            : 0;

        for( var stmt : outputs ) {
            if( fmt.appendVerbatim(stmt) )
                continue;
            var stmtX = (ExpressionStatement)stmt;
            var output = stmtX.getExpression();
            var target =
                output instanceof AssignmentExpression ae ? (VariableExpression)ae.getLeftExpression() :
                output instanceof VariableExpression ve ? ve :
                null;
            var source =
                output instanceof AssignmentExpression ae ? ae.getRightExpression() :
                null;

            if( target != null ) {
                fmt.appendLeadingComments(stmt);
                fmt.emitWrappable(() -> {
                    fmt.appendIndent();
                    visitOutputAssignment(target, source, alignmentWidth);
                    fmt.appendTrailingComment(stmt);
                    fmt.appendNewLine();
                });
            }
            else {
                fmt.visit(stmt);
            }
        }
    }

    private void visitWorkflowPublishers(List<Statement> publishers) {
        var alignmentWidth = options.harshilAlignment()
            ? maxParameterWidth(publishers)
            : 0;

        for( var stmt : publishers ) {
            if( fmt.appendVerbatim(stmt) )
                continue;
            var stmtX = (ExpressionStatement)stmt;
            var emit = (AssignmentExpression)stmtX.getExpression();
            var target = (VariableExpression)emit.getLeftExpression();
            var source = emit.getRightExpression();

            fmt.appendLeadingComments(stmt);
            fmt.emitWrappable(() -> {
                fmt.appendIndent();
                visitOutputAssignment(target, source, alignmentWidth);
                fmt.appendTrailingComment(stmt);
                fmt.appendNewLine();
            });
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

    private void visitOutputAssignment(VariableExpression target, Expression source, int alignmentWidth) {
        fmt.append(target.getText());
        if( (fmt.hasType(target) || source != null) && alignmentWidth > 0 ) {
            var padding = alignmentWidth - target.getName().length();
            fmt.append(" ".repeat(padding));
        }
        if( fmt.hasType(target) ) {
            if( alignmentWidth > 0 )
                fmt.append(' ');
            fmt.append(": ");
            fmt.visitTypeAnnotation(target.getType());
        }
        if( source != null ) {
            fmt.append(" = ");
            fmt.visitRootExpression(source);
        }
    }

    @Override
    public void visitProcessV2(ProcessNodeV2 node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("process ");
        fmt.append(node.getName());
        fmt.append(" {\n");
        fmt.incIndent();
        if( !node.directives.isEmpty() ) {
            visitDirectives(node.directives);
            fmt.appendNewLine();
        }
        var inputs = node.inputs;
        if( inputs.length > 0 ) {
            fmt.appendIndent();
            fmt.append("input:\n");
            visitTypedInputs(inputs);
            fmt.appendNewLine();
        }
        if( !node.stagers.isEmpty() ) {
            fmt.appendIndent();
            fmt.append("stage:\n");
            visitDirectives(node.stagers);
            fmt.appendNewLine();
        }
        if( !options.maheshForm() ) {
            if( !node.outputs.isEmpty() ) {
                visitProcessOutputs(node.outputs);
                fmt.appendNewLine();
            }
            if( !node.topics.isEmpty() ) {
                visitProcessTopics(node.topics);
                fmt.appendNewLine();
            }
        }
        visitProcessWhen(node.when);
        fmt.appendIndent();
        fmt.append(node.type);
        fmt.append(":\n");
        fmt.visit(node.exec);
        if( !node.stub.isEmpty() ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("stub:\n");
            fmt.visit(node.stub);
        }
        if( options.maheshForm() ) {
            if( !node.outputs.isEmpty() ) {
                fmt.appendNewLine();
                visitProcessOutputs(node.outputs);
            }
            if( !node.topics.isEmpty() ) {
                fmt.appendNewLine();
                visitProcessTopics(node.topics);
            }
        }
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    private void visitProcessWhen(Expression when) {
        if( when instanceof EmptyExpression )
            return;
        if( fmt.appendVerbatim(when) ) {
            fmt.appendNewLine();
        }
        else {
            fmt.appendLeadingComments(when);
            fmt.appendIndent();
            fmt.append("when:\n");
            fmt.emitWrappable(() -> {
                // NOTE: not visitRootExpression -- the when-expression
                // receives its own leading comments from the section
                // emission, so a wrapped chain would emit them twice
                fmt.appendIndent();
                fmt.visit(when);
                fmt.appendTrailingComment(when);
                fmt.appendNewLine();
            });
            fmt.appendDanglingAfter(when);
            fmt.appendNewLine();
        }
    }

    private void visitProcessOutputs(Statement outputs) {
        fmt.appendIndent();
        fmt.append("output:\n");
        visitTypedOutputs(asBlockStatements(outputs));
        fmt.appendDanglingComments(outputs);
    }

    private void visitProcessTopics(Statement topics) {
        fmt.appendIndent();
        fmt.append("topic:\n");
        fmt.visit(topics);
    }

    @Override
    public void visitProcessV1(ProcessNodeV1 node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("process ");
        fmt.append(node.getName());
        fmt.append(" {\n");
        fmt.incIndent();
        if( !node.directives.isEmpty() ) {
            visitDirectives(node.directives);
            fmt.appendNewLine();
        }
        if( !node.inputs.isEmpty() ) {
            fmt.appendIndent();
            fmt.append("input:\n");
            visitDirectives(node.inputs);
            fmt.appendNewLine();
        }
        if( !options.maheshForm() && !node.outputs.isEmpty() ) {
            visitProcessOutputsV1(node.outputs);
            fmt.appendNewLine();
        }
        visitProcessWhen(node.when);
        fmt.appendIndent();
        fmt.append(node.type);
        fmt.append(":\n");
        fmt.visit(node.exec);
        if( !node.stub.isEmpty() ) {
            fmt.appendNewLine();
            fmt.appendIndent();
            fmt.append("stub:\n");
            fmt.visit(node.stub);
        }
        if( options.maheshForm() && !node.outputs.isEmpty() ) {
            fmt.appendNewLine();
            visitProcessOutputsV1(node.outputs);
        }
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    private void visitProcessOutputsV1(Statement outputs) {
        fmt.appendIndent();
        fmt.append("output:\n");
        visitDirectives(outputs);
    }

    @Override
    public void visitFunction(FunctionNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("def ");
        fmt.append(node.getName());
        fmt.append('(');
        fmt.visitParameters(node.getParameters());
        fmt.append(')');
        if( fmt.hasType(node.getReturnType()) ) {
            fmt.append(" -> ");
            fmt.visitTypeAnnotation(node.getReturnType());
        }
        fmt.append(" {\n");
        fmt.incIndent();
        fmt.visit(node.getCode());
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    @Override
    public void visitRecord(RecordNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("record ");
        fmt.append(node.getName());
        visitRecordBody(node);
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    private void visitRecordBody(RecordNode node) {
        fmt.append(" {\n");
        fmt.incIndent();
        for( var fn : node.getFields() ) {
            if( fmt.appendVerbatim(fn) )
                continue;
            fmt.appendLeadingComments(fn);
            fmt.appendIndent();
            fmt.append(fn.getName());
            if( fmt.hasType(fn) ) {
                fmt.append(": ");
                fmt.visitTypeAnnotation(fn.getType());
            }
            fmt.appendTrailingComment(fn);
            fmt.appendNewLine();
        }
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.appendIndent();
        fmt.append("}");
    }

    @Override
    public void visitEnum(ClassNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("enum ");
        fmt.append(node.getName());
        fmt.append(" {\n");
        fmt.incIndent();
        for( var fn : node.getFields() ) {
            if( fmt.appendVerbatim(fn) )
                continue;
            fmt.appendLeadingComments(fn);
            fmt.appendIndent();
            fmt.append(fn.getName());
            fmt.append(',');
            fmt.appendTrailingComment(fn);
            fmt.appendNewLine();
        }
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    @Override
    public void visitOutputs(OutputBlockNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.append("output {\n");
        fmt.incIndent();
        super.visitOutputs(node);
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    @Override
    public void visitOutput(OutputNode node) {
        if( fmt.appendVerbatim(node) )
            return;
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        fmt.append(node.getName());
        if( fmt.hasType(node) ) {
            fmt.append(": ");
            fmt.visitTypeAnnotation(node.getType());
        }
        fmt.append(" {\n");
        fmt.incIndent();
        visitOutputBody((BlockStatement) node.body);
        fmt.appendDanglingComments(node);
        fmt.decIndent();
        fmt.appendIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    private void visitOutputBody(BlockStatement block) {
        asBlockStatements(block).forEach((stmt) -> {
            if( fmt.appendVerbatim(stmt) )
                return;
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
            fmt.visitDirective(call, stmt);
        });
        fmt.appendDanglingComments(block);
    }

    private void visitDirectives(Statement statement) {
        asBlockStatements(statement).forEach((stmt) -> {
            if( fmt.appendVerbatim(stmt) )
                return;
            var call = asMethodCallX(stmt);
            if( call == null )
                return;
            fmt.appendLeadingComments(stmt);
            fmt.visitDirective(call, stmt);
        });
        fmt.appendDanglingComments(statement);
    }

}
