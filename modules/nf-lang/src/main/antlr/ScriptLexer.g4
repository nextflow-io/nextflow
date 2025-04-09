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
 * https://github.com/apache/groovy/blob/GROOVY_4_0_X/src/antlr/GroovyLexer.g4
 */
lexer grammar ScriptLexer;

options {
    superClass = AbstractLexer;
}

@header {
package nextflow.script.parser;

import java.util.*;
import java.util.regex.Pattern;
import org.antlr.v4.runtime.CharStream;
import org.apache.groovy.parser.antlr4.GroovySyntaxError;

import static nextflow.script.parser.SemanticPredicates.*;
}

@members {
    private boolean errorIgnored;
    private long tokenIndex;
    private int  lastTokenType;
    private int  invalidDigitCount;

    /**
     * Record the index and token type of the current token while emitting tokens.
     */
    @Override
    public void emit(Token token) {
        this.tokenIndex++;

        int tokenType = token.getType();
        if (Token.DEFAULT_CHANNEL == token.getChannel()) {
            this.lastTokenType = tokenType;
        }

        super.emit(token);
    }

    private static final int[] REGEX_CHECK_ARRAY = {
        // DEC,
        // INC,
        // THIS,
        RBRACE,
        RBRACK,
        RPAREN,
        GStringEnd,
        TdqGStringEnd,
        NullLiteral,
        StringLiteral,
        BooleanLiteral,
        IntegerLiteral,
        FloatingPointLiteral,
        Identifier, CapitalizedIdentifier
    };
    static {
        Arrays.sort(REGEX_CHECK_ARRAY);
    }

    private boolean isRegexAllowed() {
        return (Arrays.binarySearch(REGEX_CHECK_ARRAY, this.lastTokenType) < 0);
    }

    @Override
    public int getSyntaxErrorSource() {
        return GroovySyntaxError.LEXER;
    }

    @Override
    public int getErrorLine() {
        return getLine();
    }

    @Override
    public int getErrorColumn() {
        return getCharPositionInLine() + 1;
    }

    private static boolean isJavaIdentifierStartAndNotIdentifierIgnorable(int codePoint) {
        return Character.isJavaIdentifierStart(codePoint) && !Character.isIdentifierIgnorable(codePoint);
    }

    private static boolean isJavaIdentifierPartAndNotIdentifierIgnorable(int codePoint) {
        return Character.isJavaIdentifierPart(codePoint) && !Character.isIdentifierIgnorable(codePoint);
    }

}


//
// §3.10.5 String Literals
//
StringLiteral
    :   DqStringQuotationMark  DqStringCharacter*  DqStringQuotationMark
    |   SqStringQuotationMark  SqStringCharacter*  SqStringQuotationMark
    |   Slash { this.isRegexAllowed() && _input.LA(1) != '*' }?  SlashyStringCharacter+  Slash

    |   TdqStringQuotationMark  TdqStringCharacter*  TdqStringQuotationMark
    |   TsqStringQuotationMark  TsqStringCharacter*  TsqStringQuotationMark
    ;

GStringBegin
    :   DqStringQuotationMark -> pushMode(DQ_GSTRING_MODE)
    ;
TdqGStringBegin
    :   TdqStringQuotationMark -> pushMode(TDQ_GSTRING_MODE)
    ;

mode DQ_GSTRING_MODE;
GStringEnd
    :   DqStringQuotationMark -> popMode
    ;

GStringPath
    :   Dollar IdentifierInGString (Dot IdentifierInGString)*
    ;

GStringText
    :   DqStringCharacter+
    ;

GStringExprStart
    :   '${' -> pushMode(DEFAULT_MODE)
    ;

mode TDQ_GSTRING_MODE;
TdqGStringEnd
    :   TdqStringQuotationMark -> popMode
    ;

TdqGStringPath
    :   Dollar IdentifierInGString (Dot IdentifierInGString)*
    ;

TdqGStringText
    :   TdqStringCharacter+
    ;

TdqGStringExprStart
    :   '${' -> pushMode(DEFAULT_MODE)
    ;

mode DEFAULT_MODE;
// character in the double quotation string. e.g. "a"
fragment
DqStringCharacter
    :   ~["\r\n\\$]
    |   EscapeSequence
    ;

// character in the single quotation string. e.g. 'a'
fragment
SqStringCharacter
    :   ~['\r\n\\]
    |   EscapeSequence
    ;

// character in the triple double quotation string. e.g. """a"""
fragment
TdqStringCharacter
    :   ~["\\$]
    |   DqStringQuotationMark { _input.LA(1) != '"' || _input.LA(2) != '"' || _input.LA(3) == '"' && (_input.LA(4) != '"' || _input.LA(5) != '"') }?
    |   EscapeSequence
    ;

// character in the triple single quotation string. e.g. '''a'''
fragment
TsqStringCharacter
    :   ~['\\]
    |   SqStringQuotationMark { _input.LA(1) != '\'' || _input.LA(2) != '\'' || _input.LA(3) == '\'' && (_input.LA(4) != '\'' || _input.LA(5) != '\'') }?
    |   EscapeSequence
    ;

// character in the slashy string. e.g. /a/
fragment
SlashyStringCharacter
    :   SlashEscape
    |   Dollar { !isFollowedByJavaLetterInGString(_input) }?
    |   ~[/$\u0000]
    ;


// Groovy keywords
AS              : 'as';
DEF             : 'def';
IN              : 'in';
// TRAIT           : 'trait';
// THREADSAFE      : 'threadsafe'; // reserved keyword

// the reserved type name of Java10
// VAR             : 'var';

//
// §3.9 Keywords
//
BuiltInPrimitiveType
    :   BOOLEAN
    |   CHAR
    |   BYTE
    |   SHORT
    |   INT
    |   LONG
    |   FLOAT
    |   DOUBLE
    ;

// ABSTRACT      : 'abstract';
ASSERT        : 'assert';

fragment
BOOLEAN       : 'boolean';

// BREAK         : 'break';
// YIELD         : 'yield';

fragment
BYTE          : 'byte';

// CASE          : 'case';
CATCH         : 'catch';

fragment
CHAR          : 'char';

// CLASS         : 'class';
// CONST         : 'const';
// CONTINUE      : 'continue';
// DEFAULT       : 'default';
// DO            : 'do';

fragment
DOUBLE        : 'double';

ELSE          : 'else';
ENUM          : 'enum';
// EXTENDS       : 'extends';
// FINAL         : 'final';
// FINALLY       : 'finally';

fragment
FLOAT         : 'float';

// FOR           : 'for';
IF            : 'if';
// GOTO          : 'goto';
// IMPLEMENTS    : 'implements';
IMPORT        : 'import';
INSTANCEOF    : 'instanceof';

fragment
INT           : 'int';

// INTERFACE     : 'interface';

fragment
LONG          : 'long';

// NATIVE        : 'native';
NEW           : 'new';
// NON_SEALED    : 'non-sealed';
// PACKAGE       : 'package';
// PERMITS       : 'permits';
// PRIVATE       : 'private';
// PROTECTED     : 'protected';
// PUBLIC        : 'public';
// RECORD        : 'record';
RETURN        : 'return';
// SEALED        : 'sealed';

fragment
SHORT         : 'short';

// STATIC        : 'static';
// STRICTFP      : 'strictfp';
// SUPER         : 'super';
// SWITCH        : 'switch';
// SYNCHRONIZED  : 'synchronized';
// THIS          : 'this';
THROW         : 'throw';
// THROWS        : 'throws';
// TRANSIENT     : 'transient';
TRY           : 'try';
// VOID          : 'void';
// VOLATILE      : 'volatile';
// WHILE         : 'while';

// -- feature flag, param declarations
NEXTFLOW        : 'nextflow';
PARAMS          : 'params';

// -- include declaration
INCLUDE         : 'include';
FROM            : 'from';

// -- process definition
PROCESS         : 'process';
EXEC            : 'exec';
INPUT           : 'input';
OUTPUT          : 'output';
SCRIPT          : 'script';
SHELL           : 'shell';
STUB            : 'stub';
WHEN            : 'when';

// -- workflow definition
WORKFLOW        : 'workflow';
EMIT            : 'emit';
MAIN            : 'main';
PUBLISH         : 'publish';
TAKE            : 'take';


//
// §3.10.1 Integer Literals
//
IntegerLiteral
    :   (   DecimalIntegerLiteral
        |   HexIntegerLiteral
        |   OctalIntegerLiteral
        |   BinaryIntegerLiteral
        )
        (Underscore { require(errorIgnored, "Number ending with underscores is invalid", -1, true); })?

    // !!! Error Alternative !!!
    |   Zero ([0-9] { invalidDigitCount++; })+ { require(errorIgnored, "Invalid octal number", -(invalidDigitCount + 1), true); } IntegerTypeSuffix?
    ;

fragment
Zero
    :   '0'
    ;

fragment
DecimalIntegerLiteral
    :   DecimalNumeral IntegerTypeSuffix?
    ;

fragment
HexIntegerLiteral
    :   HexNumeral IntegerTypeSuffix?
    ;

fragment
OctalIntegerLiteral
    :   OctalNumeral IntegerTypeSuffix?
    ;

fragment
BinaryIntegerLiteral
    :   BinaryNumeral IntegerTypeSuffix?
    ;

fragment
IntegerTypeSuffix
    :   [lLiIgG]
    ;

fragment
DecimalNumeral
    :   Zero
    |   NonZeroDigit (Digits? | Underscores Digits)
    ;

fragment
Digits
    :   Digit (DigitOrUnderscore* Digit)?
    ;

fragment
Digit
    :   Zero
    |   NonZeroDigit
    ;

fragment
NonZeroDigit
    :   [1-9]
    ;

fragment
DigitOrUnderscore
    :   Digit
    |   Underscore
    ;

fragment
Underscores
    :   Underscore+
    ;

fragment
Underscore
    :   '_'
    ;

fragment
HexNumeral
    :   Zero [xX] HexDigits
    ;

fragment
HexDigits
    :   HexDigit (HexDigitOrUnderscore* HexDigit)?
    ;

fragment
HexDigit
    :   [0-9a-fA-F]
    ;

fragment
HexDigitOrUnderscore
    :   HexDigit
    |   Underscore
    ;

fragment
OctalNumeral
    :   Zero Underscores? OctalDigits
    ;

fragment
OctalDigits
    :   OctalDigit (OctalDigitOrUnderscore* OctalDigit)?
    ;

fragment
OctalDigit
    :   [0-7]
    ;

fragment
OctalDigitOrUnderscore
    :   OctalDigit
    |   Underscore
    ;

fragment
BinaryNumeral
    :   Zero [bB] BinaryDigits
    ;

fragment
BinaryDigits
    :   BinaryDigit (BinaryDigitOrUnderscore* BinaryDigit)?
    ;

fragment
BinaryDigit
    :   [01]
    ;

fragment
BinaryDigitOrUnderscore
    :   BinaryDigit
    |   Underscore
    ;


//
// §3.10.2 Floating-Point Literals
//
FloatingPointLiteral
    :   (   DecimalFloatingPointLiteral
        |   HexadecimalFloatingPointLiteral
        )
        (Underscore { require(errorIgnored, "Number ending with underscores is invalid", -1, true); })?
    ;

fragment
DecimalFloatingPointLiteral
    :   Digits? Dot Digits ExponentPart? FloatTypeSuffix?
    |   Digits ExponentPart FloatTypeSuffix?
    |   Digits FloatTypeSuffix
    ;

fragment
ExponentPart
    :   ExponentIndicator SignedInteger
    ;

fragment
ExponentIndicator
    :   [eE]
    ;

fragment
SignedInteger
    :   Sign? Digits
    ;

fragment
Sign
    :   [+\-]
    ;

fragment
FloatTypeSuffix
    :   [fFdDgG]
    ;

fragment
HexadecimalFloatingPointLiteral
    :   HexSignificand BinaryExponent FloatTypeSuffix?
    ;

fragment
HexSignificand
    :   HexNumeral Dot?
    |   Zero [xX] HexDigits? Dot HexDigits
    ;

fragment
BinaryExponent
    :   BinaryExponentIndicator SignedInteger
    ;

fragment
BinaryExponentIndicator
    :   [pP]
    ;

fragment
Dot :   '.'
    ;


//
// §3.10.3 Boolean Literals
//
BooleanLiteral
    :   'true'
    |   'false'
    ;


//
// §3.10.6 Escape Sequences for Character and String Literals
//
fragment
EscapeSequence
    :   Backslash [btnfrs"'\\]
    |   OctalEscape
    |   UnicodeEscape
    |   DollarEscape
    |   LineEscape
    ;

fragment
OctalEscape
    :   Backslash OctalDigit
    |   Backslash OctalDigit OctalDigit
    |   Backslash ZeroToThree OctalDigit OctalDigit
    ;

// Groovy allows 1 or more u's after the backslash
fragment
UnicodeEscape
    :   Backslash 'u' HexDigit HexDigit HexDigit HexDigit
    ;

fragment
ZeroToThree
    :   [0-3]
    ;

// Groovy Escape Sequences
fragment
DollarEscape
    :   Backslash Dollar
    ;

fragment
LineEscape
    :   Backslash LineTerminator
    ;

fragment
LineTerminator
    :   '\r'? '\n' | '\r'
    ;

fragment
SlashEscape
    :   Backslash Slash
    ;

fragment
Backslash
    :   '\\'
    ;

fragment
Slash
    :   '/'
    ;

fragment
Dollar
    :   '$'
    ;

fragment
DqStringQuotationMark
    :   '"'
    ;

fragment
SqStringQuotationMark
    :   '\''
    ;

fragment
TdqStringQuotationMark
    :   '"""'
    ;

fragment
TsqStringQuotationMark
    :   '\'\'\''
    ;


//
// §3.10.7 The Null Literal
//
NullLiteral
    :   'null'
    ;

//
// Groovy Operators
//
RANGE_INCLUSIVE         : '..';
// RANGE_EXCLUSIVE_LEFT    : '<..';
RANGE_EXCLUSIVE_RIGHT   : '..<';
// RANGE_EXCLUSIVE_FULL    : '<..<';
SPREAD_DOT              : '*.';
SAFE_DOT                : '?.';
// SAFE_INDEX              : '?[' { this.enterParen();     } -> pushMode(DEFAULT_MODE);
// SAFE_CHAIN_DOT          : '??.';
ELVIS                   : '?:';
// METHOD_POINTER          : '.&';
// METHOD_REFERENCE        : '::';
REGEX_FIND              : '=~';
REGEX_MATCH             : '==~';
POWER                   : '**';
SPACESHIP               : '<=>';
// IDENTICAL               : '===';
// NOT_IDENTICAL           : '!==';
ARROW                   : '->';

// !internalPromise will be parsed as !in ternalPromise, so semantic predicates are necessary
NOT_INSTANCEOF      : '!instanceof' { isFollowedBy(_input, ' ', '\t', '\r', '\n') }?;
NOT_IN              : '!in'         { isFollowedBy(_input, ' ', '\t', '\r', '\n', '[', '(', '{') }?;


//
// §3.11 Separators
//
LPAREN          : '('  /* { this.enterParen();     } */ -> pushMode(DEFAULT_MODE);
RPAREN          : ')'  /* { this.exitParen();      } */ -> popMode;

LBRACE          : '{'  /* { this.enterParen();     } */ -> pushMode(DEFAULT_MODE);
RBRACE          : '}'  /* { this.exitParen();      } */ -> popMode;

LBRACK          : '['  /* { this.enterParen();     } */ -> pushMode(DEFAULT_MODE);
RBRACK          : ']'  /* { this.exitParen();      } */ -> popMode;

SEMI            : ';';
COMMA           : ',';
DOT             : Dot;


//
// §3.12 Operators
//
ASSIGN          : '=';
GT              : '>';
LT              : '<';
NOT             : '!';
BITNOT          : '~';
QUESTION        : '?';
COLON           : ':';
EQUAL           : '==';
LE              : '<=';
GE              : '>=';
NOTEQUAL        : '!=';
AND             : '&&';
OR              : '||';
// INC             : '++';
// DEC             : '--';
ADD             : '+';
SUB             : '-';
MUL             : '*';
DIV             : Slash;
BITAND          : '&';
BITOR           : '|';
XOR             : '^';
MOD             : '%';

ADD_ASSIGN      : '+=';
SUB_ASSIGN      : '-=';
MUL_ASSIGN      : '*=';
DIV_ASSIGN      : '/=';
AND_ASSIGN      : '&=';
OR_ASSIGN       : '|=';
XOR_ASSIGN      : '^=';
MOD_ASSIGN      : '%=';
LSHIFT_ASSIGN   : '<<=';
RSHIFT_ASSIGN   : '>>=';
URSHIFT_ASSIGN  : '>>>=';
ELVIS_ASSIGN    : '?=';
POWER_ASSIGN    : '**=';


//
// §3.8 Identifiers (must appear after all keywords in the grammar)
//
CapitalizedIdentifier
    :   JavaLetter { Character.isUpperCase(_input.LA(-1)) }? JavaLetterOrDigit*
    ;

Identifier
    :   JavaLetter JavaLetterOrDigit*
    ;

fragment
IdentifierInGString
    :   JavaLetterInGString JavaLetterOrDigitInGString*
    ;

fragment
JavaLetter
    :   [a-zA-Z$_] // these are the "java letters" below 0x7F
    |   // covers all characters above 0x7F which are not a surrogate
        ~[\u0000-\u007F\uD800-\uDBFF]
        { isJavaIdentifierStartAndNotIdentifierIgnorable(_input.LA(-1)) }?
    |   // covers UTF-16 surrogate pairs encodings for U+10000 to U+10FFFF
        [\uD800-\uDBFF] [\uDC00-\uDFFF]
        { Character.isJavaIdentifierStart(Character.toCodePoint((char) _input.LA(-2), (char) _input.LA(-1))) }?
    ;

fragment
JavaLetterInGString
    :   JavaLetter { _input.LA(-1) != '$' }?
    ;

fragment
JavaLetterOrDigit
    :   [a-zA-Z0-9$_] // these are the "java letters or digits" below 0x7F
    |   // covers all characters above 0x7F which are not a surrogate
        ~[\u0000-\u007F\uD800-\uDBFF]
        { isJavaIdentifierPartAndNotIdentifierIgnorable(_input.LA(-1)) }?
    |   // covers UTF-16 surrogate pairs encodings for U+10000 to U+10FFFF
        [\uD800-\uDBFF] [\uDC00-\uDFFF]
        { Character.isJavaIdentifierPart(Character.toCodePoint((char) _input.LA(-2), (char) _input.LA(-1))) }?
    ;

fragment
JavaLetterOrDigitInGString
    :   JavaLetterOrDigit { _input.LA(-1) != '$' }?
    ;

fragment
ShCommand
    :   ~[\r\n\uFFFF]*
    ;

// Additional symbols not defined in the lexical specification
// AT : '@';
// ELLIPSIS : '...';

// Whitespace, line escape and comments
WS  : ([ \t]+ | LineEscape+) -> skip
    ;

// Inside (...) and [...] but not {...}, ignore newlines.
NL  : LineTerminator   /* { this.ignoreTokenInsideParens(); } */
    ;

// Multiple-line comments (including groovydoc comments)
ML_COMMENT
    :   '/*' .*? '*/'       /* { this.ignoreMultiLineCommentConditionally(); } */ -> type(NL)
    ;

// Single-line comments
SL_COMMENT
    :   '//' ~[\r\n\uFFFF]* /* { this.ignoreTokenInsideParens(); } */             -> type(NL)
    ;

// Script-header comments.
// The very first characters of the file may be "#!".  If so, ignore the first line.
SH_COMMENT
    :   '#!' { require(errorIgnored || 0 == this.tokenIndex, "Shebang comment should appear at the first line", -2, true); } ShCommand (LineTerminator '#!' ShCommand)* -> type(NL)
    ;

// Unexpected characters will be handled by groovy parser later.
UNEXPECTED_CHAR
    :   . { require(errorIgnored, "Unexpected character: '" + getText().replace("'", "\\'") + "'", -1, false); }
    ;
