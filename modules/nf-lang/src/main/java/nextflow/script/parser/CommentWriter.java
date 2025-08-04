package nextflow.script.parser;

import org.antlr.v4.runtime.ParserRuleContext;
import org.antlr.v4.runtime.tree.TerminalNode;
import org.antlr.v4.runtime.Token;
import org.codehaus.groovy.ast.ASTNode;
import groovy.lang.Tuple2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class CommentWriter {

    public final List<Token> commentTokens;
    private final List<Token> allTokens;

    private final List<Comment> unattachedComments;
    private final List<TokenInfo> tokensWithInfo;
    private final Set<Comment> pendingComments = new HashSet<>();
    private final Map<ASTNode, List<Comment>> attachedComments = new HashMap<>();

    private final Map<ASTNode, List<Comment>> leadingComments = new HashMap<>();
    private final Map<ASTNode, Map<String, List<Comment>>> withinComments = new HashMap<>();
    private final Map<ASTNode, List<Comment>> trailingComments = new HashMap<>();

    public CommentWriter(List<Token> allTokens) {
        this.commentTokens = allTokens.stream().filter(t -> t.getChannel() == ScriptLexer.COMMENT).collect(Collectors.toList());
        this.allTokens = allTokens;

        this.tokensWithInfo = allTokens.stream().map(t -> new TokenInfo(t)).collect(Collectors.toList());
        this.unattachedComments = tokensWithInfo.stream()
            .filter(t -> t.isComment())
            .map(t -> t.getComment())
            .collect(Collectors.toList());
    }

    public class TokenInfo {
        private final Token token;
        private final Comment comment;
        
        public TokenInfo(Token token) {
            this.token = token;
            if (token.getChannel() == ScriptLexer.COMMENT) {
                this.comment = new Comment(token);
            } else {
                this.comment = null;
            }
        }

        public Token getToken() {
            return token;
        }

        public Comment getComment() {
            return comment;
        }
        
        public boolean isComment() {
            return comment != null;
        }
    }


    public String toNFComments() {
        var header = "// CommentWriter with " + commentTokens.size() + " comment tokens";
        var unattachedHeader = "// Unattached comments:";
        var indent = "    ";
        var unattachedText = unattachedComments.stream()
            .map(t -> {
                var line = t.getToken().getLine();
                var column = t.getToken().getCharPositionInLine();
                var text = t.getToken().getText().replaceAll("\n", "\\n");
                var isLeading = t.isLeading();
                var isTrailing = t.isTrailing();
                var attachedTo = t.getAttachedTo();
                var attachedToText = attachedTo != null ? attachedTo.getText() : "null";
                return (
                    "// " + indent + line + ":" + column + "{\n" + 
                    "// " + indent + indent + "content: " + text + "\n" +
                    "// " + indent + indent + "attach:  " + attachedToText + "\n" +
                    "// " + indent + "}\n"
                );
            })
            .collect(Collectors.joining("\n"));
        var pendingHeader = "// Pending comments:";
        var pendingText = pendingComments.stream()
            .map(t -> {
                var line = t.getToken().getLine();
                var column = t.getToken().getCharPositionInLine();
                var text = t.getToken().getText().replaceAll("\n", "\\n");
                return (    
                    "// " + indent + line + ":" + column + "{\n" + 
                    "// " + indent + indent + "content: " + text + "\n" +
                    "// " + indent + "}\n"
                );
            })
            .collect(Collectors.joining("\n"));
        var attachedHeader = "// Attached comments:";
        var attachedText = attachedComments.values()
            .stream()
            .flatMap(List::stream)
            .map(t -> {
                var line = t.getToken().getLine();
                var column = t.getToken().getCharPositionInLine();
                var text = t.getToken().getText().replaceAll("\n", "\\n");
                var isLeading = t.isLeading();
                var isTrailing = t.isTrailing();
                ASTNode attachedTo = t.getAttachedTo();
                var simpleName = attachedTo != null ? attachedTo.getClass().getSimpleName() : null;
                var maybeText = attachedTo != null ? attachedTo.getText(): null;
                var attachToCode = maybeText != null ? maybeText.replaceAll("\n", "\\n") : "null";
                var attachedToText = attachedTo != null ? simpleName + " " + attachToCode : "null";
                return (
                    "// " + indent + line + ":" + column + "{\n" + 
                    "// " + indent + indent + "content    : " + text + "\n" +
                    "// " + indent + indent + "attached to: " + attachedToText + "\n" +
                    "// " + indent + indent + "position   : " + t.getPositionInfo() + "\n" +
                    "// " + indent + indent + "written    : " + t.isWritten() + "\n" +
                    "// " + indent + "}\n"
                );
            })
            .collect(Collectors.joining("\n"));
        return header + "\n" + unattachedHeader + "\n" + unattachedText + "\n" + pendingHeader + "\n" + pendingText + "\n" + attachedHeader + "\n" + attachedText;
    }

    public Map<String, List<Comment>> getAttachedComments(ASTNode node) {
        var within = new ArrayList<Comment>();
        var leading = new ArrayList<Comment>();
        var trailing = new ArrayList<Comment>(); 
        attachedComments
            .getOrDefault(node, new ArrayList<>())
            .forEach(c -> {
                if( c.isWithin() ) {
                    within.add(c);
                } else if( c.isLeading() ) {
                    leading.add(c);
                } else if( c.isTrailing() ) {
                    trailing.add(c);
                }
            });
        return Map.of(
            "leading", leading,
            "within", within,
            "trailing", trailing
        );
    }

    public List<Comment> getLeadingComments(ParserRuleContext ctx) {
        if (unattachedComments.isEmpty()) {
            return new ArrayList<>();
        }

        List<Comment> comments = new ArrayList<>();
        int contextStartIndex = ctx.getStart().getTokenIndex();
        int firstCommentIndex = unattachedComments.get(0).getToken().getTokenIndex();
       
        // Get all comment tokens that are before the context and only separated by newlines
        for (int i = contextStartIndex - 1; i >= firstCommentIndex; i--) {
            TokenInfo token = tokensWithInfo.get(i);
            if (token.isComment()) {
                comments.add(token.getComment());
            } else if (token.getToken().getType() == ScriptLexer.NL) {
               continue; 
            } else {
                break;
            }
        }
        
        // Remove the comments we've processed from unattachedComments
        unattachedComments.removeAll(comments);
        
        return comments;
    }
    
    

    public List<Comment> getCommentsInbetween(
        Token start, Token end, String positionInfo, boolean isLeading, boolean isTrailing
    ) {
        List<Comment> comments = new ArrayList<>();
        Iterator<Comment> it = unattachedComments.iterator();
        int startIdx = start.getStopIndex();
        int endIdx = end.getStartIndex();
        while( it.hasNext() ) {
            Comment comment = it.next();
            if (comment.getStartIndex() >= startIdx && comment.getStopIndex() <= endIdx) {
                System.err.println(comment.getToken().getLine() + ":" + comment.getToken().getCharPositionInLine() + " " + comment.getToken().getText());
                comments.add(comment);
                comment.markAsWithin(positionInfo, isLeading, isTrailing);
                comment.makePending();
                it.remove();
            }
        }
        return comments;
    }

    public List<Comment> getWithinComments(TerminalNode start, TerminalNode end, String positionInfo, boolean isLeading, boolean isTrailing) {
        return getCommentsInbetween(start.getSymbol(), end.getSymbol(), positionInfo, isLeading, isTrailing);
    }

    public List<Comment> getTrailingComments(Token lastCtxToken) {
        // A trailing comment is either on the same line as the last token of the context,
        // or on the next line line after the last token of the context, but in this case it must
        // be separated by at least one blank line from the next context

        List<Comment> comments = new ArrayList<>();
        int contextStopIndex = lastCtxToken.getTokenIndex();
        ListIterator<TokenInfo> it = tokensWithInfo.listIterator(contextStopIndex);
        TokenInfo token = it.next(); 
        if (token.isComment()) 
            comments.add(comment); // If the token directly after the context then it is always trailing

        int nlsBefore = 0;
        List<Comment> candiateComments = new ArrayList<>();
        boolean isLeading = false;
        while( it.hasNext() && nlsBefore < 2 ) {
            TokenInfo token = it.next();
            if (token.isComment()) {
                candidateComments.add(token.getComment());
            } else if (token.getToken().getType() == ScriptLexer.NL) {
                nlsBefore++;
            } else {
                // We encountered a non newline token before there were a blank line 
                // this means that the comment should be leading to the next statement instead
                isLeading = true;
            }
        }
        if (!isLeading) {
            comments.addAll(candiatesComments);
        }
        commments.foreach(c -> {
            c.makePending();
            c.markAsTrailing();
        })
        unattachedComments.removeAll(comments); 
        return comments;
    }

    public void attachComments(ASTNode node, List<Comment> comments) {
        comments.forEach(c -> c.attachTo(node));
    }

    public List<Comment> getLeadingComments(ASTNode node) {
        return leadingComments.getOrDefault(node, new ArrayList<>());
    }

    public Map<String, List<Comment>> getWithinComments(ASTNode node) {
        return withinComments.getOrDefault(node, new HashMap<>());
    }

    public List<Comment> getTrailingComments(ASTNode node) {
        return trailingComments.getOrDefault(node, new ArrayList<>());
    }

    public class Comment {
        private final Token token;
        private boolean isLeading = false;
        private boolean isTrailing = false;
        private boolean isWithin = false;

        private ASTNode attachedTo = null;

        // If the comment is within a node, we need to keep track of where it was originally
        // We use the field below to keep track of this
        private String positionInfo = null;

        private boolean isWritten = false;

        public Comment(Token token) {
            this.token = token;
        }

        public void attachTo(ASTNode node) {
            if (attachedTo != null) {
                throw new IllegalStateException("Comment already attached to " + attachedTo);
            }
            this.attachedTo = node;
            if (pendingComments.contains(this)) {
                pendingComments.remove(this);
            } else {
                System.err.println("Comment not pending: " + token.getText());
            }
            attachedComments.computeIfAbsent(node, k -> new ArrayList<>()).add(this);
            if( isWithin ) {
                withinComments.computeIfAbsent(node, k -> new HashMap<>()).computeIfAbsent(positionInfo, k -> new ArrayList<>()).add(this);
            } else if( isLeading ) {
                leadingComments.computeIfAbsent(node, k -> new ArrayList<>()).add(this);
            } else if( isTrailing ) {
                trailingComments.computeIfAbsent(node, k -> new ArrayList<>()).add(this);
            }
        }

        public void makePending() {
            pendingComments.add(this);
        }

        public Token getToken() {
            return token;
        }

        public ASTNode getAttachedTo() {
            return attachedTo;
        }

        public void markAsLeading() {
            this.isLeading = true;
        }

        public void markAsTrailing() {
            this.isTrailing = true;
        }

        public void markAsWithin(String positionInfo, boolean isLeading, boolean isTrailing) {
            this.isWithin = true;
            this.positionInfo = positionInfo;
            this.isLeading = isLeading;
            this.isTrailing = isTrailing;
        }

        public boolean isLeading() {
            return isLeading;
        }

        public boolean isTrailing() {
            return isTrailing;
        }

        public boolean isWithin() {
            return isWithin;
        }

        public String getPositionInfo() {
            return positionInfo;
        }

        public int getStartIndex() {
            return token.getStartIndex();
        }

        public int getStopIndex() {
            return token.getStopIndex();
        }

        public boolean isWritten() {
            return isWritten;
        }

        public Tuple2<String, Boolean> write() {
            if( isWritten ) {
                throw new IllegalStateException("Comment already written");
            }
            isWritten = true;
            return new Tuple2<>(token.getText(), token.getType() == ScriptLexer.SL_COMMENT);
        }

        

        @Override
        public int hashCode() {
            return token.getTokenIndex();
        }
    }


}
