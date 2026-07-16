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
package nextflow.config.formatter;

import java.util.regex.Pattern;

import nextflow.config.ast.ConfigApplyNode;
import nextflow.config.ast.ConfigApplyBlockNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigStatement;
import nextflow.config.ast.ConfigVisitorSupport;
import nextflow.script.formatter.CommentReattacher;
import nextflow.script.formatter.Formatter;
import nextflow.script.formatter.FormattingOptions;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Format a config file.
 *
 * All comments are preserved: comment metadata is re-derived from the
 * source text with {@link CommentReattacher} and emitted as leading,
 * trailing and dangling comments through {@link Formatter}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ConfigFormattingVisitor extends ConfigVisitorSupport {

    private SourceUnit sourceUnit;

    private FormattingOptions options;

    private String sourceText;

    private Formatter fmt;

    public ConfigFormattingVisitor(SourceUnit sourceUnit, FormattingOptions options, String sourceText) {
        this.sourceUnit = sourceUnit;
        this.options = options;
        this.sourceText = sourceText;
        this.fmt = new Formatter(options);
    }

    public ConfigFormattingVisitor(SourceUnit sourceUnit, FormattingOptions options) {
        this(sourceUnit, options, null);
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ConfigNode cn ) {
            // re-derive comment metadata from the source so that no comment is lost
            if( sourceText == null ) {
                try {
                    sourceText = org.codehaus.groovy.runtime.IOGroovyMethods.getText(sourceUnit.getSource().getReader());
                }
                catch( java.io.IOException e ) {
                    // source is not available -- fall back to the comment
                    // metadata attached by the parser
                }
            }
            if( sourceText != null )
                CommentReattacher.apply(cn, sourceText);
            // enforce a blank line above config blocks and between
            // statements of different kinds; statements of the same kind
            // (e.g. assignments) may be grouped
            ConfigStatement prevStmt = null;
            for( var stmt : cn.getConfigStatements() ) {
                if( prevStmt != null && needsBlankLineBetween(prevStmt, stmt) && !leadingStartsWithBlankLine(stmt) )
                    fmt.appendNewLine();
                prevStmt = stmt;
                super.visit(stmt);
            }
            // emit comments at the end of the file
            fmt.appendDanglingComments(cn);
        }
    }

    private static boolean needsBlankLineBetween(ConfigStatement prevStmt, ConfigStatement stmt) {
        var prevBlock = isBlockStatement(prevStmt);
        var block = isBlockStatement(stmt);
        return prevBlock || block || prevStmt.getClass() != stmt.getClass();
    }

    private static boolean isBlockStatement(ConfigStatement stmt) {
        return stmt instanceof ConfigBlockNode || stmt instanceof ConfigApplyBlockNode;
    }

    private boolean leadingStartsWithBlankLine(ConfigStatement stmt) {
        var comments = (java.util.List<String>) stmt.getNodeMetaData(nextflow.script.ast.ASTNodeMarker.LEADING_COMMENTS);
        if( comments == null || comments.isEmpty() )
            return false;
        // the list is in reverse source order
        return "\n".equals(comments.get(comments.size() - 1));
    }

    public String toString() {
        return fmt.toString();
    }

    // config statements

    @Override
    public void visitConfigApplyBlock(ConfigApplyBlockNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        fmt.append(node.name);
        fmt.append(" {");
        fmt.appendNewLine();

        fmt.incIndent();
        super.visitConfigApplyBlock(node);
        fmt.appendDanglingComments(node);
        fmt.decIndent();

        fmt.appendIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    @Override
    public void visitConfigApply(ConfigApplyNode node) {
        fmt.appendLeadingComments(node);
        fmt.visitDirective(node, node);
    }

    @Override
    public void visitConfigAssign(ConfigAssignNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        var name = String.join(".", node.names);
        fmt.append(name);
        if( currentAlignmentWidth > 0 ) {
            var padding = currentAlignmentWidth - name.length();
            fmt.append(" ".repeat(padding));
        }
        fmt.append(" = ");
        fmt.visit(node.value);
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    private static final Pattern IDENTIFIER = Pattern.compile("[a-zA-Z_]+[a-zA-Z0-9_]*");

    private int currentAlignmentWidth = 0;

    @Override
    public void visitConfigBlock(ConfigBlockNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        if( node.kind != null ) {
            fmt.append(node.kind);
            fmt.append(": ");
        }
        var name = node.name;
        if( IDENTIFIER.matcher(name).matches() ) {
            fmt.append(name);
        }
        else {
            fmt.append('\'');
            fmt.append(name);
            fmt.append('\'');
        }
        fmt.append(" {");
        fmt.appendNewLine();

        int caw = currentAlignmentWidth;
        if( options.harshilAlignment() ) {
            int maxWidth = 0;
            for( var stmt : node.statements ) {
                if( stmt instanceof ConfigAssignNode can ) {
                    var width = String.join(".", can.names).length();
                    if( maxWidth < width )
                        maxWidth = width;
                }
            }
            currentAlignmentWidth = maxWidth;
        }

        fmt.incIndent();
        super.visitConfigBlock(node);
        fmt.appendDanglingComments(node);
        fmt.decIndent();

        if( options.harshilAlignment() )
            currentAlignmentWidth = caw;

        fmt.appendIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    @Override
    public void visitConfigInclude(ConfigIncludeNode node) {
        fmt.appendLeadingComments(node);
        fmt.appendIndent();
        fmt.append("includeConfig ");
        fmt.visit(node.source);
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

}
