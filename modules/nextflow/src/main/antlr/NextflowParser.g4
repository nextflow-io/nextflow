/*
 * This file is adapted from the Antlr4 Java grammar which has the following license
 *
 *  Copyright (c) 2013 Terence Parr, Sam Harwell
 *  All rights reserved.
 *  [The "BSD licence"]
 *
 *    http://www.opensource.org/licenses/bsd-license.php
 *
 * Subsequent modifications by the Groovy community have been done under the Apache License v2:
 *
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

/**
 * Grammar specification for the Nextflow scripting language.
 *
 * Based on the official grammar for Groovy:
 * https://github.com/apache/groovy/blob/GROOVY_4_0_X/src/antlr/GroovyParser.g4
 */
parser grammar NextflowParser;

options {
    superClass = AbstractParser;
    tokenVocab = NextflowLexer;
}

@header {
package nextflow.antlr;

import org.apache.groovy.parser.antlr4.GroovySyntaxError;
}

@members {

    @Override
    public int getSyntaxErrorSource() {
        return GroovySyntaxError.PARSER;
    }

    @Override
    public int getErrorLine() {
        Token token = _input.LT(-1);

        if (null == token) {
            return -1;
        }

        return token.getLine();
    }

    @Override
    public int getErrorColumn() {
        Token token = _input.LT(-1);

        if (null == token) {
            return -1;
        }

        return token.getCharPositionInLine() + 1 + token.getText().length();
    }
}


compilationUnit
    :   nls (scriptStatement (sep scriptStatement)* sep?)? EOF
    ;

//
// top-level statements
//
scriptStatement
    :   includeStatement            #includeStmtAlt
    |   processDef                  #processDefAlt
    |   workflowDef                 #workflowDefAlt
    |   functionDef                 #functionDefAlt
    ;

// -- include statement
includeStatement
    :   INCLUDE includeNames FROM stringLiteral
    ;

includeNames
    :   LBRACE includeName (SEMI includeName)* RBRACE
    ;

includeName
    :   name=identifier
    |   name=identifier AS alias=identifier
    ;

// -- process definition
processDef
    :   PROCESS name=identifier
        nls LBRACE
        (sep processDirective)*
        (sep processInputs)?
        (sep processOutputs)?
        (sep processWhen)?
        sep processExec
        sep? RBRACE
    ;

processDirective
    :   identifier argumentList?
    ;

processInputs
    :   PROCESS_INPUT (sep processDirective)+
    ;

processOutputs
    :   PROCESS_OUTPUT (sep processDirective)+
    ;

processWhen
    :   PROCESS_WHEN nls expression
    ;

processExec
    // explicit script/exec body with optional stub
    :   (   PROCESS_SCRIPT
        |   PROCESS_SHELL
        |   PROCESS_EXEC
        )
        sep? blockStatements
        (sep processStub)?

    // implicit script body
    |   blockStatements
    ;

processStub
    :   PROCESS_STUB sep? blockStatements
    ;

// -- workflow definition
workflowDef
    :   WORKFLOW name=identifier?
        nls LBRACE
        sep? workflowBody
        sep? RBRACE
    ;

workflowBody
    // explicit main block with option take/emit blocks
    :   (   sep WORKFLOW_TAKE
            workflowTakes
        )?
        sep WORKFLOW_MAIN
        sep? workflowMain
        (   (sep WORKFLOW_EMIT)
            workflowEmits
        )?

    // implicit main block
    |   workflowMain
    ;

workflowTakes:
    (sep identifier)+
    ;

workflowMain:
    blockStatements
    ;

workflowEmits:
    (sep identifier)+
    ;

// -- function definition
functionDef
    :   (DEF | type | DEF type) nls
        identifier LPAREN formalParameterList? rparen nls
        LBRACE nls blockStatements? RBRACE
    ;


//
// statements
//
statement
    :   ifElseStatement             #ifElseStmtAlt
    |   RETURN expression?          #returnStmtAlt
    |   assertStatement             #assertStmtAlt
    |   variableDeclaration         #variableDeclarationStmtAlt
    |   multipleAssignmentStatement #multipleAssignmentStmtAlt
    |   assignmentStatement         #assignmentStmtAlt
    |   expressionStatement         #expressionStmtAlt
    |   SEMI                        #emptyStmtAlt
    ;

// -- if/else statement
ifElseStatement
    :   IF parExpression nls tb=ifElseBranch (nls ELSE nls fb=ifElseBranch)?
    ;

ifElseBranch
    :   LBRACE nls blockStatements? RBRACE
    |   statement
    ;

blockStatements
    :   statement (sep statement)* nls
    ;

// -- assert statement
assertStatement
    :   ASSERT condition=expression (nls (COLON | COMMA) nls message=expression)?
    ;

// -- variable declaration
variableDeclaration
    :   DEF type? variableDeclarator
    |   DEF typeNamePairs nls ASSIGN nls initializer=expression
    |   type variableDeclarator
    ;

variableDeclarator
    :   identifier (nls ASSIGN nls initializer=expression)?
    ;

typeNamePairs
    :   LPAREN typeNamePair (COMMA typeNamePair)* rparen
    ;

typeNamePair
    :   type? identifier
    ;

// -- assignment statement
// "(a) = [1]" is a special case of multipleAssignmentStatement, it will be handled by assignmentStatement
multipleAssignmentStatement
    :   <assoc=right>
        left=variableNames nls
        op=ASSIGN nls
        right=expression
    ;

assignmentStatement
    :   <assoc=right>
        left=expression nls
        op=(ASSIGN
        |   ADD_ASSIGN
        |   SUB_ASSIGN
        |   MUL_ASSIGN
        |   DIV_ASSIGN
        |   AND_ASSIGN
        |   OR_ASSIGN
        |   XOR_ASSIGN
        |   RSHIFT_ASSIGN
        |   URSHIFT_ASSIGN
        |   LSHIFT_ASSIGN
        |   MOD_ASSIGN
        |   POWER_ASSIGN
        |   ELVIS_ASSIGN
        ) nls
        right=expression
    ;

// -- expression statement
expressionStatement
    :   expression
        (
            { !SemanticPredicates.isFollowingArgumentsOrClosure($expression.ctx) }?
            argumentList
        |
            /* if expression is a method call, no need to have any more arguments */
        )
    ;


//
// expressions
//
expression
    // must come before postfix expression to resolve the ambiguities between casting and call on parentheses expression, e.g. (int)(1 / 2)
    :   LPAREN type rparen castOperandExpression                                            #castExprAlt

    // postfix (++/--)
    |   pathExpression op=(INC | DEC)                                                       #postfixExprAlt

    // qualified name, list/map element, method invocation
    |   pathExpression                                                                      #pathExprAlt

    // bitwise not (~) / logical not (!) (level 1)
    |   op=(BITNOT | NOT) nls expression                                                    #unaryNotExprAlt

    // math power operator (**) (level 2)
    |   left=expression op=POWER nls right=expression                                       #powerExprAlt

    // prefix (++/--) (level 3)
    |   op=(INC | DEC) expression                                                           #prefixExprAlt

    // unary (+/-) (level 3)
    |   op=(ADD | SUB) expression                                                           #unaryAddExprAlt

    // multiplication/division/modulo (level 4)
    |   left=expression nls op=(MUL | DIV | MOD) nls right=expression                       #multDivExprAlt

    // binary addition/subtraction (level 5)
    |   left=expression op=(ADD | SUB) nls right=expression                                 #addSubExprAlt

    // bit shift, range (level 6)
    |   left=expression nls
        ((  dlOp=LT LT
        |   tgOp=GT GT GT
        |   dgOp=GT GT
        )
        |(  riOp=RANGE_INCLUSIVE
        |   reOp=RANGE_EXCLUSIVE_RIGHT
        )) nls
        right=expression                                                                    #shiftExprAlt

    // boolean relational expressions (level 7)
    |   left=expression nls op=AS nls type                                                  #relationalCastExprAlt
    |   left=expression nls op=INSTANCEOF nls type                                          #relationalTypeExprAlt
    |   left=expression nls op=(LE | GE | GT | LT | IN) nls right=expression                #relationalExprAlt

    // equality/inequality (==/!=) (level 8)
    |   left=expression nls
        op=(EQUAL
        |   NOTEQUAL
        |   SPACESHIP
        ) nls
        right=expression                                                                    #equalityExprAlt

    // regex find and match (=~ and ==~) (level 8.5)
    |   left=expression nls op=(REGEX_FIND | REGEX_MATCH) nls right=expression              #regexExprAlt

    // bitwise and (&)  (level 9)
    |   left=expression nls op=BITAND nls right=expression                                  #bitwiseAndExprAlt

    // exclusive or (^)  (level 10)
    |   left=expression nls op=XOR nls right=expression                                     #exclusiveOrExprAlt

    // bitwise or (|)  (level 11)
    |   left=expression nls op=BITOR nls right=expression                                   #bitwiseOrExprAlt

    // logical and (&&)  (level 12)
    |   left=expression nls op=AND nls right=expression                                     #logicalAndExprAlt

    // logical or (||)  (level 13)
    |   left=expression nls op=OR nls right=expression                                      #logicalOrExprAlt

    // ternary, elvis (level 14)
    |   <assoc=right>
        condition=expression nls
        (   QUESTION nls tb=expression nls COLON nls
        |   ELVIS nls
        )
        fb=expression                                                                       #conditionalExprAlt
    ;

castOperandExpression
    :   LPAREN type rparen castOperandExpression            #castCastExprAlt

    |   pathExpression op=(INC | DEC)                       #postfixCastExprAlt
    |   pathExpression                                      #pathCastExprAlt

    // bitwise not (~) / logical not (!)
    |   op=(BITNOT | NOT) nls castOperandExpression         #unaryNotCastExprAlt

    // prefix (++/--)
    |   op=(INC | DEC) castOperandExpression                #prefixCastExprAlt

    // unary (+/-)
    |   op=(ADD | SUB) castOperandExpression                #unaryAddCastExprAlt
    ;

variableNames
    :   LPAREN identifier (COMMA identifier)+ rparen
    ;

// -- path expression
pathExpression
    :   primary pathElement*
    ;

primary
    :   identifier                  #identifierPrmrAlt
    |   literal                     #literalPrmrAlt
    |   gstring                     #gstringPrmrAlt
    |   NEW nls creator             #newPrmrAlt
    |   parExpression               #parenPrmrAlt
    |   closure                     #closurePrmrAlt
    |   list                        #listPrmrAlt
    |   map                         #mapPrmrAlt
    |   builtInType                 #builtInTypePrmrAlt
    ;

pathElement
    // property expression
    :   nls
        (   DOT                 // dot operator
        |   SPREAD_DOT          // spread operator:         x*.y === x?.collect { it.y }
        |   SAFE_DOT            // optional-null operator:  x?.y === (x!=null) ? x.y : null
        )
        nls namePart                                    #propertyPathExprAlt

    // method call expression (with closure)
    |   closure                                         #closurePathExprAlt

    // method call expression
    |   arguments                                       #argumentsPathExprAlt

    // index expression
    |   indexPropertyArgs                               #indexPathExprAlt
    ;

namePart
    :   identifier
    |   stringLiteral
    |   keywords
    ;

indexPropertyArgs
    :   LBRACK expressionList RBRACK
    ;

// -- variable, type identifiers
identifier
    :   Identifier
    |   CapitalizedIdentifier
    |   IN
    ;

// -- primitive literals
literal
    :   IntegerLiteral          #integerLiteralAlt
    |   FloatingPointLiteral    #floatingPointLiteralAlt
    |   stringLiteral           #stringLiteralAlt
    |   BooleanLiteral          #booleanLiteralAlt
    |   NullLiteral             #nullLiteralAlt
    ;

stringLiteral
    :   StringLiteral
    ;

// -- gstring expression
gstring
    :   GStringBegin gstringDqPart* GStringEnd
    |   TdqGStringBegin gstringTdqPart* TdqGStringEnd
    ;

gstringDqPart
    :   GStringText                         #gstringDqTextAlt
    |   GStringPath                         #gstringDqPathAlt
    |   GStringExprStart expression RBRACE  #gstringDqExprAlt
    ;

gstringTdqPart
    :   TdqGStringText                          #gstringTdqTextAlt
    |   TdqGStringPath                          #gstringTdqPathAlt
    |   TdqGStringExprStart expression RBRACE   #gstringTdqExprAlt
    ;

// -- constructor method call
creator
    :   createdName nls arguments
    ;

createdName
    :   primitiveType
    |   qualifiedClassName typeArgumentsOrDiamond?
    ;

typeArgumentsOrDiamond
    :   LT GT
    |   typeArguments
    ;

// -- parenthetical expression
parExpression
    :   LPAREN expression rparen
    ;

// -- closure expression
closure
    :   LBRACE (nls (formalParameterList nls)? ARROW)? nls blockStatements? RBRACE
    ;

formalParameterList
    :   formalParameter (COMMA nls formalParameter)*
    ;

formalParameter
    :   DEF? type? ELLIPSIS? identifier (nls ASSIGN nls expression)?
    ;

// -- list expression
list
    :   LBRACK expressionList? COMMA? RBRACK
    ;

expressionList
    :   expressionListElement (COMMA nls expressionListElement)*
    ;

expressionListElement
    :   MUL? expression
    ;

// -- map expression
map
    :   LBRACK
        (   mapEntryList COMMA?
        |   COLON
        )
        RBRACK
    ;

mapEntryList
    :   mapEntry (COMMA mapEntry)*
    ;

mapEntry
    :   mapEntryLabel COLON nls expression
    |   MUL COLON nls expression
    ;

mapEntryLabel
    :   keywords
    |   primary
    ;

// -- primitive type
builtInType
    :   BuiltInPrimitiveType
    ;

// -- argument list
arguments
    :   LPAREN argumentList? COMMA? rparen
    ;


argumentList
    :   argumentListElement
        (   COMMA nls
            argumentListElement
        )*
    ;

argumentListElement
    :   expressionListElement
    |   namedArg
    ;

namedArg
    :   namedArgLabel COLON nls expression
    |   MUL COLON nls expression
    ;

namedArgLabel
    :   keywords
    |   identifier
    |   literal
    |   gstring
    ;

//
// types
//
type
    :   primitiveType
    |   qualifiedClassName typeArguments?
    ;

primitiveType
    :   BuiltInPrimitiveType
    ;

qualifiedClassName
    :   qualifiedNameElements className
    ;

qualifiedNameElements
    :   (qualifiedNameElement DOT)*
    ;

qualifiedNameElement
    :   identifier
    |   AS
    |   DEF
    |   FROM
    |   IN
    ;

className
    :   CapitalizedIdentifier
    ;

typeArguments
    :   LT nls type (COMMA nls type)* nls GT
    ;


//
// keywords, whitespace
//
keywords
    :   AS
    |   DEF
    |   FROM
    |   IN
    |   INSTANCEOF
    |   RETURN
    |   NullLiteral
    |   BooleanLiteral
    |   BuiltInPrimitiveType
    ;

rparen
    :   RPAREN
    ;

nls
    :   NL*
    ;

sep :   (NL | SEMI)+
    ;
