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
/*
 *  Licensed to the Apache Software Foundation (ASF) under one
 *  or more contributor license agreements.  See the NOTICE file
 *  distributed with this work for additional information
 *  regarding copyright ownership.  The ASF licenses this file
 *  to you under the Apache License, Version 2.0 (the
 *  "License"); you may not use this file except in compliance
 *  with the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an
 *  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 *  KIND, either express or implied.  See the License for the
 *  specific language governing permissions and limitations
 *  under the License.
 */
package nextflow.script.parser;

import groovy.lang.groovydoc.Groovydoc;
import groovy.lang.groovydoc.GroovydocHolder;
import org.antlr.v4.runtime.ParserRuleContext;
import org.codehaus.groovy.GroovyBugError;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.AnnotationNode;

import java.util.regex.Pattern;

import static org.apache.groovy.parser.antlr4.util.StringUtils.matches;
import static org.codehaus.groovy.runtime.DefaultGroovyMethods.asBoolean;

/**
 * Extract Groovydoc comments from source code and
 * add it as metadata to AST nodes.
 */
class GroovydocManager {

    private static final String GROOVYDOC_PREFIX = "/**";
    private static final Pattern SPACES_PATTERN = Pattern.compile("\\s+");

    private final boolean groovydocEnabled;

    public GroovydocManager(boolean groovydocEnabled) {
        this.groovydocEnabled = groovydocEnabled;
    }

    /**
     * Attach doc comment to member node as meta data
     *
     * @param node
     * @param ctx
     */
    public void handle(ASTNode node, ParserRuleContext ctx) {
        if( !groovydocEnabled )
            return;
        if( !asBoolean(node) || !asBoolean(ctx) )
            return;
        if( !(node instanceof GroovydocHolder) )
            return;

        var docCommentNodeText = findDocCommentByNode(ctx);
        if( docCommentNodeText == null )
            return;

        node.putNodeMetaData(GroovydocHolder.DOC_COMMENT, new Groovydoc(docCommentNodeText, (GroovydocHolder) node));
    }

    private String findDocCommentByNode(ParserRuleContext ctx) {
        var parent = ctx.getParent();
        if( !asBoolean(parent) )
            return null;

        String docCommentNodeText = null;
        boolean sameTypeNodeBefore = false;
        for( var child : parent.children ) {
            if( child == ctx ) {
                // if no doc comment ctx found and no siblings of same type before the ctx,
                // try to find doc comment ctx of its parent
                if( docCommentNodeText == null && !sameTypeNodeBefore )
                    return findDocCommentByNode(parent);
                else
                    return docCommentNodeText;
            }

            if( child.getClass() == ctx.getClass() ) {
                docCommentNodeText = null;
                sameTypeNodeBefore = true;
                continue;
            }

            if( !(child instanceof ScriptParser.NlsContext) && !(child instanceof ScriptParser.SepContext) )
                continue;

            // doc comments are treated as NL tokens
            var newlines = child instanceof ScriptParser.NlsContext
                ? ((ScriptParser.NlsContext) child).NL()
                : ((ScriptParser.SepContext) child).NL();

            if( newlines.isEmpty() )
                continue;

            for( int i = newlines.size() - 1; i >= 0; i-- ) {
                var text = newlines.get(i).getText();
                if( matches(text, SPACES_PATTERN) )
                    continue;
                docCommentNodeText = text.startsWith(GROOVYDOC_PREFIX)
                    ? text
                    : null;
                break;
            }
        }

        throw new GroovyBugError("Groovydoc context can not be found: " + ctx.getText());
    }
}
