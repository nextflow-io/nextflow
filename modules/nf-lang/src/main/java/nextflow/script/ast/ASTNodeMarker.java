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
package nextflow.script.ast;

/**
 * Additional markers for AST nodes that are used for static analysis.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public enum ASTNodeMarker {
    // denotes a fully-qualified type annotation (ClassNode)
    FULLY_QUALIFIED,

    // denotes that an assignment is an implicit declaration
    IMPLICIT_DECLARATION,

    // the inferred type of an expression
    INFERRED_TYPE,

    // the number of enclosing parentheses around an expression
    INSIDE_PARENTHESES_LEVEL,

    // the comments preceding a statement or declaration
    LEADING_COMMENTS,

    // the verbatim text of a Groovy-style type annotation (ClassNode)
    LEGACY_TYPE,

    // the MethodNode targeted by a MethodCallExpression
    METHOD_TARGET,

    // the MethodNode targeted by a variable expression (PropertyNode)
    METHOD_VARIABLE_TARGET,

    // the starting quote sequence of a string literal or gstring expression
    QUOTE_CHAR,

    // denotes that an expression list has a trailing comma
    TRAILING_COMMA,

    // the trailing comment on the same line as a statement or declaration
    TRAILING_COMMENT,

    // the verbatim text of a string literal or gstring expression
    VERBATIM_TEXT
}
