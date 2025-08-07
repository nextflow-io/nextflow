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
    private final Map<ASTNode, Map<String, List<Comment>>> trailingComments = new HashMap<>();

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

        @Override
        public String toString() {
            return "TokenInfo:" + (isComment() ? "COMMENT" : "SEMANTIC") + ": '" +  token.getText().replace("\n", "\\n") + "'";
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
                var relPos = (
                    (t.isWithin() ? "within " : "")+
                    (t.isLeading() ? "leading " : "") + 
                    (t.isTrailing() ? "trailing " : "")
                );
                return (
                    "// " + indent + line + ":" + column + "{\n" + 
                    "// " + indent + indent + "content    : " + text + "\n" +
                    "// " + indent + indent + "attached to: " + attachedToText + "\n" +
                    "// " + indent + indent + "token      : " + t.getPositionInfo() + "\n" +
                    "// " + indent + indent + "written    : " + t.isWritten() + "\n" +
                    "// " + indent + indent + "rel. pos.  : " + relPos + "\n" +
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

    public List<Comment> processLeadingComments(ParserRuleContext ctx) {
        if (unattachedComments.isEmpty()) {
            return new ArrayList<>();
        }

        List<Comment> comments = new ArrayList<>();
        int contextStartIndex = ctx.getStart().getTokenIndex();
        int firstCommentIndex = unattachedComments.get(0).getToken().getTokenIndex();
       
        // Get all comment tokens that are before the context and only separated by newlines
        for (int i = contextStartIndex - 1; i >= firstCommentIndex; i--) {
            TokenInfo token = tokensWithInfo.get(i);
            if( token.isComment() && !token.getComment().isConsumed() ) {
                comments.add(0, token.getComment());
            } else if( token.getToken().getType() == ScriptLexer.NL ){
               continue; 
            } else {
                break;
            }
        }

        comments.forEach(c -> {
            c.makePending();
            c.markAsLeading(ctx.getStart());
        });
        
        return comments;
    }

    public List<Comment> processInbetweenComments(
        Token start, Token end, String positionInfo, boolean isLeading, boolean isTrailing
    ) {
        List<Comment> comments = new ArrayList<>();
        Iterator<Comment> it = unattachedComments.iterator();
        int startIdx = start.getStopIndex();
        int endIdx = end.getStartIndex();
        while( it.hasNext() ) {
            Comment comment = it.next();
            if (comment.getStartIndex() >= startIdx && comment.getStopIndex() <= endIdx) {
                comments.add(comment);
            }
        }
        comments.forEach(c -> {
            c.markAsWithin(positionInfo, start, end, isLeading, isTrailing);
            c.makePending();

        });
        return comments;
    }

    public List<Comment> processInbetweenComments(TerminalNode start, TerminalNode end, String positionInfo, boolean isLeading, boolean isTrailing) {
        return processInbetweenComments(start.getSymbol(), end.getSymbol(), positionInfo, isLeading, isTrailing);
    }


    public List<Comment> processTrailingComments(Token lastCtxToken, String positionInfo, boolean allowNewlines) {
        // A trailing comment is either on the same line as the last token of the context,
        // or on the next line line after the last token of the context, but in this case it must
        // be separated by at least one blank line from the next context
        List<Comment> comments = new ArrayList<>();
        int contextStopIndex = lastCtxToken.getTokenIndex();
        ListIterator<TokenInfo> it = tokensWithInfo.listIterator(contextStopIndex + 1);
        TokenInfo token = null;
        int nlsBefore = 0;
        // First get tokens that are directly after the context, without any newlines
        while( it.hasNext() ) {
            token = it.next(); 
            if (token.isComment() && !token.getComment().isConsumed()) {
                comments.add(token.getComment());
            } else {
                if (token.getToken().getType() == ScriptLexer.NL) {
                    nlsBefore++;
                }
                break;
            }
        }

        if (allowNewlines) {
            // Get comments that are on the lines following the context
            // but only if they seem to be trailing rather than leading 
            // for the next context
            List<Comment> candidateComments = new ArrayList<>();
            boolean isLeading = false;
            while( it.hasNext() && nlsBefore < 2 ) {
                token = it.next();
                if (token.isComment() && !token.getComment().isConsumed()) {
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
                comments.addAll(candidateComments);
            }
        }

        comments.forEach(c -> {
            c.makePending();
            c.markAsTrailing(lastCtxToken, positionInfo);
        });

        return comments;
    }

    public List<Comment> processTrailingComments(TerminalNode term, String positionInfo, boolean allowNewlines) {
        return processTrailingComments(term.getSymbol(), positionInfo, allowNewlines);
    }

    public List<Comment> processTrailingComments(TerminalNode term, boolean allowNewlines) {
        return processTrailingComments(term.getSymbol(), "STANDARD", allowNewlines);
    }

    public List<Comment> processTrailingComments(Token token, boolean allowNewlines) {
        return processTrailingComments(token, "STANDARD", allowNewlines);
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

    public Map<String, List<Comment>> getTrailingComments(ASTNode node) {
        return trailingComments.getOrDefault(node, new HashMap<>());
    }

    public class Comment {
        private final Token token;
        private boolean isLeading = false;
        private boolean isTrailing = false;
        private boolean isWithin = false;

        private ASTNode attachedTo = null;
        private boolean consumed = false;

        // If the comment is within a node, we need to keep track of where it was originally
        // We use the field below to keep track of this
        private String positionInfo = null;
        private Token precedingToken = null;
        private Token followingToken = null;

        private boolean isWritten = false;

        public Comment(Token token) {
            this.token = token;
        }

        public void attachTo(ASTNode node) {
            if (attachedTo != null) {
                throw new IllegalStateException(this + " already attached to " + attachedTo);
            }
            this.attachedTo = node;
            pendingComments.remove(this);
            attachedComments.computeIfAbsent(node, k -> new ArrayList<>()).add(this);
            if( isWithin ) {
                withinComments.computeIfAbsent(node, k -> new HashMap<>()).computeIfAbsent(positionInfo, k -> new ArrayList<>()).add(this);
            } else if( isLeading ) {
                leadingComments.computeIfAbsent(node, k -> new ArrayList<>()).add(this);
            } else if( isTrailing ) {
                trailingComments.computeIfAbsent(node, k -> new HashMap<>()).computeIfAbsent(positionInfo, k -> new ArrayList<>()).add(this);
            }
        }

        public void makePending() {
            if (consumed) {
                throw new IllegalStateException(this + " is already consumed");
            }
            unattachedComments.remove(this);
            pendingComments.add(this);
            consumed = true;            
        }

        public Token getToken() {
            return token;
        }

        public ASTNode getAttachedTo() {
            return attachedTo;
        }

        public void markAsLeading(Token followingToken) {
            this.isLeading = true;
            this.followingToken = followingToken;
        }

        public void markAsTrailing(Token precedingToken, String positionInfo) {
            this.isTrailing = true;
            this.precedingToken = precedingToken;
            this.positionInfo = positionInfo;
        }

        public void markAsWithin(
            String positionInfo,
            Token start,
            Token end,
            boolean isLeading,
            boolean isTrailing
        ) {
            this.isWithin = true;
            this.positionInfo = positionInfo;
            this.precedingToken = start;
            this.followingToken = end;
            this.isLeading = isLeading;
            this.isTrailing = isTrailing;
        }

        public boolean isConsumed() {
            return consumed;
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
                throw new IllegalStateException(this + " already written");
            }
            isWritten = true;
            return new Tuple2<>(token.getText(), token.getType() == ScriptLexer.SL_COMMENT);
        }

        @Override 
        public String toString() {
            return String.format(
                "Comment:%d:%d:'%s'",
                token.getLine(),
                token.getCharPositionInLine(),
                token.getText().replaceAll("\n", "\\n")
            );
        }
        @Override
        public int hashCode() {
            return token.getTokenIndex();
        }
    }


}
