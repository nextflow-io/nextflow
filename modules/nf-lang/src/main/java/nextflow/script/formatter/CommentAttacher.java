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
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import nextflow.config.ast.ConfigApplyBlockNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigStatement;
import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.OutputBlockNode;
import nextflow.script.ast.ParamBlockNode;
import nextflow.script.ast.ProcessNodeV1;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.ModuleNode;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.IfStatement;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.ast.stmt.TryCatchStatement;
import org.codehaus.groovy.syntax.Types;

import static nextflow.script.ast.ASTUtils.*;
import static nextflow.script.formatter.Comment.*;

/**
 * Attach every comment in a source file to an AST node, so that
 * the formatter can print it back out.
 *
 * The AST is organized into a tree of containers -- constructs that
 * are printed with a set of braces. Each container owns the nodes that
 * can carry a comment directly. A comment is resolved against the
 * innermost container that contains it:
 *
 * - end-of-line comment after a child node -> trailing on that node
 * - inside the extent of a child node -> hoisted to leading on that node
 * - before a child node -> leading on that node
 * - otherwise -> dangling on the container itself
 *
 * The final rule is the reason no comment can be lost: every comment is
 * inside some container, and the module itself is a container.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class CommentAttacher {

    private final Comments comments;

    private final Map<ASTNode,List<Comment>> leading = new IdentityHashMap<>();
    private final Map<ASTNode,List<Comment>> trailing = new IdentityHashMap<>();
    private final Map<ASTNode,List<Comment>> dangling = new IdentityHashMap<>();

    private final Set<Comment> emitted = Collections.newSetFromMap(new IdentityHashMap<>());

    private CommentAttacher(Comments comments) {
        this.comments = comments;
    }

    /**
     * Derive the comment attachments for a module. The attachments are owned
     * by the caller (i.e. a single formatting pass) rather than stored in node
     * metadata, since the same AST may be formatted many times.
     *
     * @param moduleNode
     */
    public static CommentAttacher of(ModuleNode moduleNode) {
        var comments = moduleNode != null
            ? (Comments) moduleNode.getNodeMetaData(ASTNodeMarker.COMMENTS)
            : null;
        if( comments == null )
            comments = Comments.EMPTY;
        var result = new CommentAttacher(comments);
        if( !comments.getComments().isEmpty() )
            result.attach(moduleNode);
        return result;
    }

    public List<Comment> leading(ASTNode node) {
        return get(leading, node);
    }

    public List<Comment> trailing(ASTNode node) {
        return get(trailing, node);
    }

    public List<Comment> dangling(ASTNode node) {
        return get(dangling, node);
    }

    public boolean hasComments(ASTNode node) {
        return !leading(node).isEmpty() || !trailing(node).isEmpty() || !dangling(node).isEmpty();
    }

    public int blankLinesBefore(int line) {
        return comments.blankLinesBefore(line);
    }

    /**
     * Determine whether a node is the first thing on its source line, i.e.
     * whether its column is an indentation rather than a position within a
     * line of code.
     */
    public boolean startsLine(ASTNode node) {
        return hasPosition(node) && comments.startsLine(node.getLineNumber(), node.getColumnNumber());
    }

    public void markEmitted(Comment comment) {
        emitted.add(comment);
    }

    /**
     * Get the comments that were attached but never printed. Should always
     * be empty -- a non-empty result is a formatter bug that would silently
     * delete source.
     */
    public List<Comment> getMissingComments() {
        return comments.getComments().stream()
            .filter(c -> !emitted.contains(c))
            .toList();
    }

    private static List<Comment> get(Map<ASTNode,List<Comment>> map, ASTNode node) {
        if( node == null )
            return Collections.emptyList();
        var result = map.get(node);
        return result != null ? result : Collections.emptyList();
    }

    private void add(Map<ASTNode,List<Comment>> map, ASTNode node, Comment comment) {
        map.computeIfAbsent(node, (k) -> new ArrayList<>()).add(comment);
    }

    /// ATTACHMENT

    private void attach(ModuleNode moduleNode) {
        var root = new Container(moduleNode, Long.MIN_VALUE, Long.MAX_VALUE, false);
        if( moduleNode instanceof ScriptNode sn )
            script(sn, root);
        else if( moduleNode instanceof ConfigNode cn )
            configStatements(cn.getConfigStatements(), root);
        for( var comment : comments.getComments() )
            attach(comment, root);
    }

    private void attach(Comment comment, Container container) {
        for( var nested : container.nested ) {
            if( nested.contains(comment) ) {
                attach(comment, nested);
                return;
            }
        }

        Slot preceding = null;
        Slot following = null;
        Slot enclosing = null;
        for( var child : container.children ) {
            if( child.end <= comment.start() ) {
                if( preceding == null || child.end >= preceding.end )
                    preceding = child;
            }
            else if( child.start >= comment.end() ) {
                if( following == null || child.start < following.start )
                    following = child;
            }
            else if( child.start <= comment.start() && comment.end() <= child.end ) {
                // innermost enclosing child
                if( enclosing == null || child.start >= enclosing.start )
                    enclosing = child;
            }
        }

        if( !comment.ownLine && preceding != null && preceding.node.getLastLineNumber() == comment.line ) {
            add(trailing, preceding.node, comment);
        }
        else if( !comment.ownLine && enclosing != null && enclosing.trailing != null
                && enclosing.trailing.getLastLineNumber() == comment.line ) {
            add(trailing, enclosing.trailing, comment);
        }
        else if( enclosing != null ) {
            add(leading, enclosing.node, comment);
        }
        else if( following != null ) {
            add(leading, following.node, comment);
        }
        else if( !comment.ownLine && preceding != null ) {
            add(trailing, preceding.node, comment);
        }
        else if( container.hoist ) {
            add(leading, container.node, comment);
        }
        else {
            add(dangling, container.node, comment);
        }
    }

    /**
     * A child of a container: a node together with the source range in which
     * a comment belongs to it. The range is usually the node's own extent,
     * but for a method-chain link it is the gap between the receiver and the
     * method name, where the formatter can break the line.
     */
    private record Slot(ASTNode node, long start, long end, ASTNode trailing) {

        Slot(ASTNode node, long start, long end) {
            this(node, start, end, null);
        }
    }

    /**
     * A construct that can own comments which belong to no child node.
     */
    private static class Container {
        final ASTNode node;
        final long start;
        final long end;
        /** emit leftover comments as leading on the node instead of dangling */
        final boolean hoist;
        final List<Slot> children = new ArrayList<>();
        final List<Container> nested = new ArrayList<>();

        Container(ASTNode node, long start, long end, boolean hoist) {
            this.node = node;
            this.start = start;
            this.end = end;
            this.hoist = hoist;
        }

        boolean contains(Comment comment) {
            return start <= comment.start() && comment.end() <= end;
        }
    }

    private Container container(ASTNode node, Container parent) {
        return container(node, start(node), end(node), parent);
    }

    private Container container(ASTNode node, long start, long end, Container parent) {
        return container(node, start, end, parent, false);
    }

    private Container container(ASTNode node, long start, long end, Container parent, boolean hoist) {
        var result = new Container(node, start, end, hoist);
        parent.nested.add(result);
        return result;
    }

    private void addChild(Container container, ASTNode node) {
        if( hasPosition(node) )
            container.children.add(new Slot(node, start(node), end(node)));
    }

    private void addChild(Container container, ASTNode node, long start, long end, ASTNode trailing) {
        if( start < end )
            container.children.add(new Slot(node, start, end, trailing));
    }

    /// SCRIPT DECLARATIONS

    private void script(ScriptNode node, Container root) {
        for( var decl : node.getDeclarations() ) {
            // a code snippet is an implicit entry workflow with no position
            if( decl instanceof WorkflowNode wn && wn.isCodeSnippet() ) {
                statements(wn.main, root);
                continue;
            }
            if( !hasPosition(decl) )
                continue;
            addChild(root, decl);
            scriptDeclaration(decl, root);
        }
    }

    private void scriptDeclaration(ASTNode decl, Container parent) {
        if( decl instanceof ClassNode cn && cn.isEnum() ) {
            var container = container(cn, parent);
            for( var field : cn.getFields() )
                addChild(container, field);
        }
        else if( decl instanceof FunctionNode fn ) {
            statements(fn.getCode(), container(fn, parent));
        }
        else if( decl instanceof OutputBlockNode obn ) {
            var container = container(obn, parent);
            for( var on : obn.declarations ) {
                addChild(container, on);
                statements(on.body, container(on, container));
            }
        }
        else if( decl instanceof ParamBlockNode pbn ) {
            var container = container(pbn, parent);
            for( var param : pbn.declarations )
                addChild(container, param);
        }
        else if( decl instanceof ProcessNodeV1 pn ) {
            var container = container(pn, parent);
            statements(pn.directives, container);
            statements(pn.inputs, container);
            statements(pn.outputs, container);
            statements(pn.exec, container);
            statements(pn.stub, container);
        }
        else if( decl instanceof ProcessNodeV2 pn ) {
            var container = container(pn, parent);
            statements(pn.directives, container);
            for( var input : pn.inputs )
                addChild(container, input);
            statements(pn.stagers, container);
            statements(pn.outputs, container);
            statements(pn.topics, container);
            statements(pn.exec, container);
            statements(pn.stub, container);
        }
        else if( decl instanceof RecordNode rn ) {
            var container = container(rn, parent);
            for( var field : rn.getFields() )
                addChild(container, field);
        }
        else if( decl instanceof WorkflowNode wn ) {
            var container = container(wn, parent);
            for( var take : wn.getParameters() )
                addChild(container, take);
            statements(wn.main, container);
            statements(wn.emits, container);
            statements(wn.publishers, container);
            statements(wn.onComplete, container);
            statements(wn.onError, container);
        }
    }

    /// STATEMENTS

    private void statements(Statement node, Container container) {
        if( node == null )
            return;
        for( var stmt : asBlockStatements(node) ) {
            addChild(container, stmt);
            statement(stmt, container);
        }
    }

    private void statement(Statement node, Container parent) {
        if( node instanceof IfStatement is ) {
            var elseBlock = is.getElseBlock();
            var boundary = hasPosition(elseBlock)
                ? boundary(is.getIfBlock(), start(elseBlock))
                : end(is);
            var thenContainer = container(is.getIfBlock(), start(is), boundary, parent);
            statements(is.getIfBlock(), thenContainer);
            is.getBooleanExpression().visit(new ExpressionWalker(thenContainer));
            if( hasPosition(elseBlock) ) {
                if( elseBlock instanceof IfStatement )
                    statement(elseBlock, parent);
                else
                    statements(elseBlock, container(elseBlock, boundary, end(is), parent));
            }
        }
        else if( node instanceof TryCatchStatement tcs ) {
            var catches = tcs.getCatchStatements();
            var tryBlock = tcs.getTryStatement();
            var boundary = !catches.isEmpty() && hasPosition(catches.get(0))
                ? boundary(tryBlock, start(catches.get(0)))
                : end(tcs);
            statements(tryBlock, container(tryBlock, start(tcs), boundary, parent));
            for( int i = 0; i < catches.size(); i++ ) {
                var cs = catches.get(i);
                if( !hasPosition(cs) )
                    continue;
                var end = i + 1 < catches.size() && hasPosition(catches.get(i + 1))
                    ? start(catches.get(i + 1))
                    : end(tcs);
                statements(cs.getCode(), container(cs.getCode(), start(cs), end, parent));
            }
        }
        else {
            var chain = chainRoot(node);
            var container = chain != null
                ? container(node, start(node), end(node), parent, true)
                : parent;
            if( chain != null )
                chainSlots(chain, container);
            node.visit(new ExpressionWalker(container));
        }
    }

    /**
     * The root method chain of a statement, if the formatter is able to wrap
     * it -- only a chain at the root of a statement or assignment value is
     * wrapped, and only when it has at least two links.
     */
    private static MethodCallExpression chainRoot(Statement node) {
        Expression expr =
            node instanceof ExpressionStatement es ? es.getExpression() :
            node instanceof ReturnStatement rs ? rs.getExpression() :
            null;
        if( expr instanceof BinaryExpression be && be.getOperation().isA(Types.ASSIGNMENT_OPERATOR) )
            expr = be.getRightExpression();
        if( !(expr instanceof MethodCallExpression root) )
            return null;
        int depth = 0;
        while( expr instanceof MethodCallExpression mce && !mce.isImplicitThis() ) {
            expr = mce.getObjectExpression();
            depth++;
        }
        return depth >= 2 ? root : null;
    }

    /**
     * Register a slot for each link of a method chain: the gap between the
     * end of the receiver and the start of the method name, which is where
     * the formatter breaks the line when the chain is wrapped. A comment
     * there leads the link, or trails the receiver if it is on the same line.
     */
    private void chainSlots(Expression root, Container container) {
        Expression expr = root;
        while( expr instanceof MethodCallExpression mce && !mce.isImplicitThis() ) {
            var receiver = mce.getObjectExpression();
            var method = mce.getMethod();
            if( hasPosition(receiver) && hasPosition(method) )
                addChild(container, mce, end(receiver), start(method), receiver);
            expr = receiver;
        }
    }

    /**
     * Determine where a block ends and the next one begins. The braces and
     * keywords between them are not in the AST, so the split is made at the
     * end of the line that closes the first block: an end-of-line comment
     * stays with the block it follows, an own-line comment leads the next
     * block.
     */
    private static long boundary(Statement block, long next) {
        if( !hasPosition(block) )
            return next;
        return Math.min(Comment.position(block.getLastLineNumber() + 1, 0), next);
    }

    /// EXPRESSIONS

    /**
     * Walk an expression tree, creating a container for each construct
     * that is printed with its own set of brackets.
     */
    private class ExpressionWalker extends CodeVisitorSupport {

        private Container current;

        ExpressionWalker(Container current) {
            this.current = current;
        }

        @Override
        public void visitClosureExpression(ClosureExpression node) {
            if( !hasPosition(node) ) {
                super.visitClosureExpression(node);
                return;
            }
            var prev = current;
            current = container(node, prev);
            statements(node.getCode(), current);
            current = prev;
        }

        @Override
        public void visitListExpression(ListExpression node) {
            var prev = current;
            if( hasPosition(node) )
                current = container(node, prev);
            for( var element : node.getExpressions() ) {
                addChild(current, element);
                element.visit(this);
            }
            current = prev;
        }

        @Override
        public void visitMapExpression(MapExpression node) {
            var prev = current;
            if( hasPosition(node) )
                current = container(node, prev);
            visitMapEntries(node);
            current = prev;
        }

        @Override
        public void visitMethodCallExpression(MethodCallExpression node) {
            var args = asMethodCallArguments(node);
            // a trailing closure is printed outside the parentheses, so it is
            // not an argument for the purpose of comment attachment
            var lastClosureArg = !args.isEmpty() && args.get(args.size() - 1) instanceof ClosureExpression;
            var parenArgs = lastClosureArg ? args.subList(0, args.size() - 1) : args;
            var prev = current;
            if( !parenArgs.isEmpty() && hasPosition(node) ) {
                // the receiver has its own containers, and the gap before the
                // method name belongs to the enclosing method chain
                var from = !node.isImplicitThis() && hasPosition(node.getMethod())
                    ? start(node.getMethod())
                    : start(node);
                current = container(node, from, end(node), prev);
            }
            if( !node.isImplicitThis() )
                node.getObjectExpression().visit(this);
            for( var arg : parenArgs ) {
                if( arg instanceof NamedArgumentListExpression nale ) {
                    visitMapEntries(nale);
                }
                else {
                    addChild(current, arg);
                    arg.visit(this);
                }
            }
            if( lastClosureArg )
                args.get(args.size() - 1).visit(this);
            current = prev;
        }

        private void visitMapEntries(MapExpression node) {
            for( var entry : node.getMapEntryExpressions() ) {
                addChild(current, entry);
                entry.getValueExpression().visit(this);
            }
        }
    }

    /// CONFIG STATEMENTS

    private void configStatements(List<ConfigStatement> statements, Container parent) {
        for( var stmt : statements ) {
            if( !hasPosition(stmt) )
                continue;
            addChild(parent, stmt);
            if( stmt instanceof ConfigApplyBlockNode cabn ) {
                var container = container(cabn, parent);
                for( var apply : cabn.statements )
                    addChild(container, apply);
            }
            else if( stmt instanceof ConfigAssignNode can ) {
                can.value.visit(new ExpressionWalker(parent));
            }
            else if( stmt instanceof ConfigBlockNode cbn ) {
                configStatements(cbn.statements, container(cbn, parent));
            }
            else if( stmt instanceof ConfigIncludeNode cin ) {
                cin.source.visit(new ExpressionWalker(parent));
            }
        }
    }

}
