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

import java.util.List;
import java.util.function.Consumer;
import java.util.regex.Pattern;

import nextflow.config.ast.ConfigApplyNode;
import nextflow.config.ast.ConfigApplyBlockNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigStatement;
import nextflow.config.ast.ConfigVisitorSupport;
import nextflow.script.formatter.Comment;
import nextflow.script.formatter.CommentAttacher;
import nextflow.script.formatter.FmtDirectives;
import nextflow.script.formatter.FormattingOptions;
import nextflow.script.formatter.Formatter;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Format a config file.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ConfigFormattingVisitor extends ConfigVisitorSupport {

    private SourceUnit sourceUnit;

    private FormattingOptions options;

    private Formatter fmt;

    public ConfigFormattingVisitor(SourceUnit sourceUnit, FormattingOptions options) {
        this.sourceUnit = sourceUnit;
        this.options = options;
        this.fmt = new Formatter(options);
        var commentAttacher = CommentAttacher.of(sourceUnit.getAST());
        this.fmt.setComments(commentAttacher);
        this.fmt.setFmtDirectives(FmtDirectives.of(sourceUnit.getAST(), () -> FmtDirectives.readSourceText(sourceUnit), commentAttacher));
    }

    /**
     * Get the comments that were not printed. Should always be empty --
     * a non-empty result means formatting would delete source.
     */
    public List<Comment> getMissingComments() {
        return fmt.getComments().getMissingComments();
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ConfigNode cn ) {
            visitStatements(cn.getConfigStatements(), (stmt) -> visit(stmt));
            fmt.appendDanglingComments(cn);
        }
    }

    public String toString() {
        return fmt.toString();
    }

    /**
     * Emit a list of sibling config statements, with the blank line above
     * each one owned by this loop rather than by the statement's own
     * (source-derived) leading-comment logic -- mirrors the top-level
     * declaration loop in {@code ScriptFormattingVisitor}. Used both for the
     * top level of the file and for the body of a config block, since config
     * has no separate "top-level" concept.
     *
     * @param statements
     */
    private <T extends ASTNode> void visitStatements(List<T> statements, Consumer<T> visit) {
        T prev = null;
        for( var stmt : statements ) {
            // a statement suppressed inside a verbatim region emits nothing,
            // so it must not get a blank line above it
            if( !fmt.isSuppressed(stmt) )
                fmt.blankLines(blankLinesBetween(prev, stmt));
            visit.accept(stmt);
            prev = stmt;
        }
    }

    /**
     * The blank lines that belong above a config statement, given its
     * previous sibling (or {@code null} for the first statement of the file
     * or block): a config block is set off by exactly 1 blank line (config is
     * denser than a script's top-level definitions, which get 2); two
     * statements of the same kind (e.g. consecutive assignments) keep their
     * source grouping, capped at 1; of different kinds, exactly 1.
     *
     * @param prev
     * @param decl
     */
    private int blankLinesBetween(ASTNode prev, ASTNode decl) {
        if( prev == null )
            return 0;
        if( isConfigBlock(prev) || isConfigBlock(decl) )
            return 1;
        if( prev.getClass() != decl.getClass() )
            return 1;
        var comments = fmt.getComments();
        return Math.min(1, comments.blankLinesBefore(comments.leadingLine(decl)));
    }

    private static boolean isConfigBlock(ASTNode node) {
        return node instanceof ConfigBlockNode || node instanceof ConfigApplyBlockNode;
    }

    // config statements

    @Override
    public void visitConfigApplyBlock(ConfigApplyBlockNode node) {
        if( fmt.appendVerbatimInner(node) )
            return;
        fmt.appendInnerComments(node);
        fmt.appendIndent();
        fmt.append(node.name);
        fmt.append(" {");
        fmt.appendNewLine();
        fmt.markBlockStart();

        fmt.incIndent();
        visitStatements(node.statements, (stmt) -> visitConfigApply(stmt));
        fmt.appendDanglingComments(node);
        fmt.decIndent();

        fmt.appendIndent();
        fmt.append('}');
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

    @Override
    public void visitConfigApply(ConfigApplyNode node) {
        if( fmt.appendVerbatimInner(node) )
            return;
        fmt.appendInnerComments(node);
        fmt.visitDirective(node);
    }

    @Override
    public void visitConfigAssign(ConfigAssignNode node) {
        if( fmt.appendVerbatimInner(node) )
            return;
        fmt.appendInnerComments(node);
        fmt.emitWrappable(() -> {
            fmt.appendIndent();
            var name = String.join(".", node.names);
            fmt.append(name);
            if( currentAlignmentWidth > 0 ) {
                var padding = currentAlignmentWidth - name.length();
                fmt.append(" ".repeat(padding));
            }
            fmt.append(" = ");
            var sid = fmt.beginStatement(node);
            fmt.visit(node.value);
            fmt.endStatement(sid);
            fmt.appendTrailingComment(node);
            fmt.appendNewLine();
        });
    }

    private static final Pattern IDENTIFIER = Pattern.compile("[a-zA-Z_]+[a-zA-Z0-9_]*");

    private int currentAlignmentWidth = 0;

    @Override
    public void visitConfigBlock(ConfigBlockNode node) {
        if( fmt.appendVerbatimInner(node) )
            return;
        fmt.appendInnerComments(node);
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
        fmt.markBlockStart();

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
        visitStatements(node.statements, (stmt) -> visit(stmt));
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
        if( fmt.appendVerbatimInner(node) )
            return;
        fmt.appendInnerComments(node);
        fmt.appendIndent();
        fmt.append("includeConfig ");
        fmt.visit(node.source);
        fmt.appendTrailingComment(node);
        fmt.appendNewLine();
    }

}
