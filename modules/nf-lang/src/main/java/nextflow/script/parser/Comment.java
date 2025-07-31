package nextflow.script.parser;

import org.antlr.v4.runtime.Token;
import org.codehaus.groovy.ast.ASTNode;

public class Comment {
    private final Token token;
    private final ASTNode attachedTo;
    private final boolean isLeading;
    private final boolean isTrailing;

    public Comment(Token token, ASTNode attachedTo, boolean isLeading, boolean isTrailing) {
        this.token = token;
        this.attachedTo = attachedTo;
        this.isLeading = isLeading;
        this.isTrailing = isTrailing;
    }

    public Token getToken() {
        return token;
    }

    public ASTNode getAttachedTo() {
        return attachedTo;
    }

    public boolean isLeading() {
        return isLeading;
    }

    public boolean isTrailing() {
        return isTrailing;
    }
}
