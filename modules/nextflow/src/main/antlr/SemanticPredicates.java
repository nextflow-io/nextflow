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
package nextflow.antlr;

import java.util.regex.Pattern;

import org.antlr.v4.runtime.CharStream;
import org.codehaus.groovy.GroovyBugError;

import static nextflow.antlr.NextflowParser.ArgumentsPathExprAltContext;
import static nextflow.antlr.NextflowParser.ClosurePathExprAltContext;
import static nextflow.antlr.NextflowParser.ExpressionContext;
import static nextflow.antlr.NextflowParser.PathExprAltContext;
import static org.apache.groovy.parser.antlr4.util.StringUtils.matches;

/**
 * Some semantic predicates for altering the behaviour of the lexer and parser
 */
public class SemanticPredicates {
    private static final Pattern NONSPACES_PATTERN = Pattern.compile("\\S+?");
    private static final Pattern LETTER_AND_LEFTCURLY_PATTERN = Pattern.compile("[a-zA-Z_{]");
    private static final Pattern NONSURROGATE_PATTERN = Pattern.compile("[^\u0000-\u007F\uD800-\uDBFF]");
    private static final Pattern SURROGATE_PAIR1_PATTERN = Pattern.compile("[\uD800-\uDBFF]");
    private static final Pattern SURROGATE_PAIR2_PATTERN = Pattern.compile("[\uDC00-\uDFFF]");

    public static boolean isFollowedByWhiteSpaces(CharStream cs) {
        for (int index = 1, c = cs.LA(index); !('\r' == c || '\n' == c || CharStream.EOF == c); index++, c = cs.LA(index)) {
            if (matches(String.valueOf((char) c), NONSPACES_PATTERN)) {
                return false;
            }
        }

        return true;
    }

    public static boolean isFollowedBy(CharStream cs, char... chars) {
        int c1 = cs.LA(1);

        for (char c : chars) {
            if (c1 == c) {
                return true;
            }
        }

        return false;
    }

    public static boolean isFollowedByJavaLetterInGString(CharStream cs) {
        int c1 = cs.LA(1);

        if ('$' == c1) { // single $ is not a valid identifier
            return false;
        }

        String str1 = String.valueOf((char) c1);

        if (matches(str1, LETTER_AND_LEFTCURLY_PATTERN)) {
            return true;
        }

        if (matches(str1, NONSURROGATE_PATTERN)
                && Character.isJavaIdentifierPart(c1)) {
            return true;
        }

        int c2 = cs.LA(2);
        String str2 = String.valueOf((char) c2);

        if (matches(str1, SURROGATE_PAIR1_PATTERN)
                && matches(str2, SURROGATE_PAIR2_PATTERN)
                && Character.isJavaIdentifierPart(Character.toCodePoint((char) c1, (char) c2))) {

            return true;
        }

        return false;
    }

    /**
     * Check whether following a method name of command expression.
     * Method name should not end with "2: arguments" or "3: closure"
     *
     * @param context the preceding expression
     */
    public static boolean isFollowingArgumentsOrClosure(ExpressionContext context) {
        if (context instanceof PathExprAltContext)
            return false;

        try {
            var pathExpression = (PathExprAltContext) context;
            var pathElement = pathExpression.children.get(0);
            return pathElement instanceof ClosurePathExprAltContext || pathElement instanceof ArgumentsPathExprAltContext;
        } catch (IndexOutOfBoundsException | ClassCastException e) {
            throw new GroovyBugError("Unexpected structure of expression context: " + context, e);
        }
    }

}
