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
parser grammar ScriptParser;

options {
    superClass = AbstractParser;
    tokenVocab = ScriptLexer;
}

@header {
package nextflow.script.parser;

import org.apache.groovy.parser.antlr4.GroovySyntaxError;

import static nextflow.script.parser.SemanticPredicates.*;
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
    :   nls (scriptDeclarationOrStatement (sep scriptDeclarationOrStatement)* sep?)? EOF
    ;

scriptDeclarationOrStatement
    :   scriptDeclaration
    |   statement
    ;

//
// script declarations
//
scriptDeclaration
    :   featureFlagDeclaration      #featureFlagDeclAlt
    |   includeDeclaration          #includeDeclAlt
    |   importDeclaration           #importDeclAlt
    |   paramDeclaration            #paramDeclAlt
    |   enumDef                     #enumDefAlt
    |   processDef                  #processDefAlt
    |   workflowDef                 #workflowDefAlt
    |   outputDef                   #outputDefAlt
    |   functionDef                 #functionDefAlt
    |   incompleteScriptDeclaration #incompleteScriptDeclAlt
    ;

// -- feature flag declaration
featureFlagDeclaration
    :   featureFlagName nls ASSIGN nls expression
    ;

featureFlagName
    :   NEXTFLOW (DOT identifier)+
    ;

// -- include declaration
includeDeclaration
    :   INCLUDE includeNames FROM stringLiteral
    ;

includeNames
    :   LBRACE nls includeName (sep includeName)* sep? RBRACE
    ;

includeName
    :   name=identifier
    |   name=identifier AS alias=identifier
    ;

// -- import declaration (legacy)
importDeclaration
    :   IMPORT qualifiedClassName
    ;

// -- param declaration
paramDeclaration
    :   PARAMS (DOT identifier)+ nls ASSIGN nls expression
    ;

// -- enum definition
enumDef
    :   ENUM identifier nls LBRACE
        nls enumBody? COMMA?
        nls RBRACE
    ;

enumBody
    :   identifier (nls COMMA nls identifier)*
    ;

// -- process definition
processDef
    :   PROCESS name=identifier nls LBRACE
        body=processBody?
        sep? RBRACE
    ;

processBody
    // explicit script/exec body with optional stub
    :   (sep processDirectives)?
        (sep processInputs)?
        (sep processOutputs)?
        (sep processWhen)?
        sep processExec
        (sep processStub)?

    // explicit "Mahesh" form
    |   (sep processDirectives)?
        (sep processInputs)?
        sep processExec
        (sep processStub)?
        sep processOutputs

    // implicit script/exec body
    |   (sep processDirectives)?
        (sep processInputs)?
        (sep processOutputs)?
        (sep processWhen)?
        sep blockStatements
    ;

processDirectives
    :   statement (sep statement)*
    ;

processInputs
    :   INPUT COLON nls statement (sep statement)*
    ;

processOutputs
    :   OUTPUT COLON nls statement (sep statement)*
    ;

processWhen
    :   WHEN COLON nls expression
    ;

processExec
    :   (SCRIPT | SHELL | EXEC) COLON nls blockStatements
    ;

processStub
    :   STUB COLON nls blockStatements
    ;

// -- workflow definition
workflowDef
    :   WORKFLOW name=identifier? nls LBRACE
        body=workflowBody?
        sep? RBRACE
    ;

workflowBody
    // explicit main block with optional take/emit blocks
    :   (sep TAKE COLON nls workflowTakes)?
        sep MAIN COLON nls workflowMain
        (sep EMIT COLON nls workflowEmits)?
        (sep PUBLISH COLON nls workflowPublishers)?

    // explicit emit block with optional take/main blocks
    |   (sep TAKE COLON nls workflowTakes)?
        (sep MAIN COLON nls workflowMain)?
        sep EMIT COLON nls workflowEmits
        (sep PUBLISH COLON nls workflowPublishers)?

    // implicit main block
    |   sep? workflowMain
    ;

workflowTakes
    :   identifier (sep identifier)*
    ;

workflowMain
    :   blockStatements
    ;

workflowEmits
    :   statement (sep statement)*
    ;

workflowPublishers
    :   statement (sep statement)*
    ;

// -- output definition
outputDef
    :   OUTPUT nls LBRACE
        outputBody?
        sep? RBRACE
    ;

outputBody
    :   sep? outputTargetBody (sep outputTargetBody)*
    ;

outputTargetBody
    : statement
    ;

// -- function definition
functionDef
    :   (DEF | legacyType | DEF legacyType) identifier LPAREN nls (formalParameterList nls)? rparen nls LBRACE
        nls blockStatements? RBRACE
    ;

// -- incomplete script declaration
incompleteScriptDeclaration
    :   identifier (DOT identifier)* DOT?
    ;


//
// statements
//
statement
    :   ifElseStatement                 #ifElseStmtAlt
    |   tryCatchStatement               #tryCatchStmtAlt
    |   RETURN expression?              #returnStmtAlt
    |   THROW expression                #throwStmtAlt
    |   assertStatement                 #assertStmtAlt
    |   variableDeclaration             #variableDeclarationStmtAlt
    |   multipleAssignmentStatement     #multipleAssignmentStmtAlt
    |   assignmentStatement             #assignmentStmtAlt
    |   expressionStatement             #expressionStmtAlt
    |   SEMI                            #emptyStmtAlt
    ;

// -- if/else statement
ifElseStatement
    :   IF parExpression nls tb=statementOrBlock (nls ELSE nls fb=statementOrBlock)?
    ;

statementOrBlock
    :   LBRACE nls blockStatements? RBRACE
    |   statement
    ;

blockStatements
    :   statement (sep statement)* sep?
    ;

// -- try/catch statement
tryCatchStatement
    :   TRY nls statementOrBlock (nls catchClause)*
    ;

catchClause
    :   CATCH LPAREN catchTypes? identifier rparen nls statementOrBlock
    ;

catchTypes
    :   qualifiedClassName (BITOR qualifiedClassName)*
    ;

// -- assert statement
assertStatement
    :   ASSERT condition=expression (nls COLON nls message=expression)?
    ;

// -- variable declaration
variableDeclaration
    :   (DEF | legacyType | DEF legacyType) identifier (nls ASSIGN nls initializer=expression)?
    |   DEF variableNames nls ASSIGN nls initializer=expression
    ;

variableNames
    :   LPAREN identifier (COMMA identifier)+ rparen
    ;

// -- assignment statement
multipleAssignmentStatement
    :   variableNames nls ASSIGN nls expression
    ;

assignmentStatement
    :   target=expression nls
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
        source=expression
    ;

// -- expression statement
expressionStatement
    :   expression
        (
            { isValidDirective($expression.ctx) }? argumentList
        |
            /* only certain expressions can be called as a directive (no parens) */
        )
    ;


//
// expressions
//
expression
    // identifiers, literals, closures, lists, maps, method calls, index/property expressions
    :   primary pathElement*                                                                #pathExprAlt

    // bitwise not (~) / logical not (!) (level 1)
    |   op=(BITNOT | NOT) nls expression                                                    #unaryNotExprAlt

    // math power operator (**) (level 2)
    |   left=expression op=POWER nls right=expression                                       #powerExprAlt

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
    |   left=expression nls op=(INSTANCEOF | NOT_INSTANCEOF) nls type                       #relationalTypeExprAlt
    |   left=expression nls op=(LE | GE | GT | LT | IN | NOT_IN) nls right=expression       #relationalExprAlt

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

    // incomplete expression
    |   expression nls (DOT | SPREAD_DOT | SAFE_DOT)                                        #incompleteExprAlt
    ;

primary
    :   identifier                  #identifierPrmrAlt
    |   literal                     #literalPrmrAlt
    |   gstring                     #gstringPrmrAlt
    |   NEW creator                 #newPrmrAlt
    |   parExpression               #parenPrmrAlt
    |   closure                     #closurePrmrAlt
    |   list                        #listPrmrAlt
    |   map                         #mapPrmrAlt
    |   builtInType                 #builtInTypePrmrAlt
    ;

pathElement
    // property expression
    :   nls
        (   DOT
        |   SPREAD_DOT  // spread dot:  xs*.y == xs?.collect { x -> x.y }
        |   SAFE_DOT    // safe dot:    x?.y == (x != null) ? x.y : null
        )
        namedProperty                                   #propertyPathExprAlt

    // method call expression (with closure)
    |   closure                                         #closurePathExprAlt
    |   closureWithLabels                               #closureWithLabelsPathExprAlt

    // method call expression
    |   arguments                                       #argumentsPathExprAlt

    // index expression
    |   indexPropertyArgs                               #indexPathExprAlt
    ;

namedProperty
    :   identifier
    |   stringLiteral
    |   keywords
    ;

indexPropertyArgs
    :   LBRACK expressionList RBRACK
    ;

// -- variable, function, type identifiers
identifier
    :   Identifier
    |   CapitalizedIdentifier
    |   IN
    |   NEXTFLOW
    |   PARAMS
    |   FROM
    |   PROCESS
    |   EXEC
    |   INPUT
    |   OUTPUT
    |   SCRIPT
    |   SHELL
    |   STUB
    |   WHEN
    |   WORKFLOW
    |   EMIT
    |   MAIN
    |   PUBLISH
    |   TAKE
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
    :   createdName arguments
    ;

createdName
    :   primitiveType
    |   qualifiedClassName typeArguments?
    ;

// -- parenthetical expression
parExpression
    :   LPAREN nls expression nls rparen
    ;

// -- closure expression
closure
    :   LBRACE (nls (formalParameterList nls)? ARROW)? nls blockStatements? RBRACE
    ;

formalParameterList
    :   formalParameter (COMMA nls formalParameter)*
    ;

formalParameter
    :   DEF? legacyType? identifier (nls ASSIGN nls expression)?
    ;

closureWithLabels
    :   LBRACE (nls (formalParameterList nls)? ARROW)? nls blockStatementsWithLabels RBRACE
    ;

blockStatementsWithLabels
    :   statementOrLabeled (sep statementOrLabeled)* sep?
    ;

statementOrLabeled
    :   identifier COLON nls statementOrLabeled
    |   statement
    ;

// -- list expression
list
    :   LBRACK nls expressionList? COMMA? nls RBRACK
    ;

expressionList
    :   expression (nls COMMA nls expression)*
    ;

// -- map expression
map
    :   LBRACK nls mapEntryList COMMA? nls RBRACK
    |   LBRACK COLON RBRACK
    ;

mapEntryList
    :   mapEntry (nls COMMA nls mapEntry)*
    ;

mapEntry
    :   mapEntryLabel COLON expression
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
    :   LPAREN nls argumentList? COMMA? nls rparen
    ;

argumentList
    :   argumentListElement (nls COMMA nls argumentListElement)*
    ;

argumentListElement
    :   expression
    |   namedArg
    ;

namedArg
    :   namedProperty COLON expression
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
    |   IN
    ;

className
    :   CapitalizedIdentifier
    ;

typeArguments
    :   LT type (COMMA type)* GT
    ;

legacyType
    :   type (LBRACK RBRACK)*
    ;


//
// keywords, whitespace
//
keywords
    :   AS
    |   DEF
    |   IMPORT
    |   IN
    |   INSTANCEOF
    |   RETURN
    |   NEXTFLOW
    |   PARAMS
    |   INCLUDE
    |   FROM
    |   PROCESS
    |   EXEC
    |   INPUT
    |   OUTPUT
    |   SCRIPT
    |   SHELL
    |   STUB
    |   WHEN
    |   WORKFLOW
    |   EMIT
    |   MAIN
    |   PUBLISH
    |   TAKE
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
