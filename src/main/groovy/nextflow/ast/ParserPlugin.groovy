/*
 * Copyright (c) 2012, the authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.ast

import java.lang.reflect.Field

import antlr.LexerSharedInputState
import antlr.RecognitionException
import antlr.Token
import antlr.TokenStreamException
import antlr.TokenStreamRecognitionException
import groovy.util.logging.Slf4j
import org.codehaus.groovy.antlr.AntlrParserPlugin
import org.codehaus.groovy.antlr.GroovySourceToken
import org.codehaus.groovy.antlr.SourceBuffer
import org.codehaus.groovy.antlr.UnicodeEscapingReader
import org.codehaus.groovy.antlr.UnicodeLexerSharedInputState
import org.codehaus.groovy.antlr.parser.GroovyLexer
import org.codehaus.groovy.antlr.parser.GroovyRecognizer
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.ParserPlugin
import org.codehaus.groovy.control.ParserPluginFactory
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.Reduction
import org.codehaus.groovy.syntax.SyntaxException

/**
 * This parser plugin is required to handle the 'as' keyword in the Nextflow DSL
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class SourceModifierParserPlugin extends AntlrParserPlugin {

    @groovy.transform.PackageScope
    Closure<CustomLexer> lexerFactory

    @groovy.transform.PackageScope
    CustomLexer lexer

    Reduction parseCST(SourceUnit sourceUnit, Reader reader) throws CompilationFailedException {
        super.parseCST(sourceUnit, reader)
    }

    protected void transformCSTIntoAST(SourceUnit sourceUnit, Reader reader, SourceBuffer sourceBuffer) throws CompilationFailedException {
        ast = null;

        setController(sourceUnit);

        UnicodeEscapingReader unicodeReader = new UnicodeEscapingReader(reader, sourceBuffer);
        UnicodeLexerSharedInputState inputState = new UnicodeLexerSharedInputState(unicodeReader);
        // create the lexer using the factory closure when available
        lexer = lexerFactory ? lexerFactory.call(inputState) : new CustomLexer(inputState)
        unicodeReader.setLexer(lexer);
        // create the parser
        GroovyRecognizer parser = GroovyRecognizer.make(lexer);
        parser.setSourceBuffer(sourceBuffer);
        setTokenNames(parser.getTokenNames())
        parser.setFilename(sourceUnit.getName());

        // start parsing at the compilationUnit rule
        try {
            parser.compilationUnit();
        }
        catch (TokenStreamRecognitionException tsre) {
            RecognitionException e = tsre.recog;
            SyntaxException se = new SyntaxException(e.getMessage(), e, e.getLine(), e.getColumn());
            se.setFatal(true);
            sourceUnit.addError(se);
        }
        catch (RecognitionException e) {
            SyntaxException se = new SyntaxException(e.getMessage(), e, e.getLine(), e.getColumn());
            se.setFatal(true);
            sourceUnit.addError(se);
        }
        catch (TokenStreamException e) {
            sourceUnit.addException(e);
        }

        ast = parser.getAST();
    }

    protected void setTokenNames(String[] tokenNames) {

        Class clazz = this.class.superclass
        Field field = clazz.getDeclaredField('tokenNames')
        field.setAccessible(true)
        field.set(this, tokenNames)
    }

}

@Slf4j
class CustomLexer extends GroovyLexer {

    CustomLexer(LexerSharedInputState state) { super(state)  }

    CustomLexer(InputStream stream) {  super(stream)  }

    CustomLexer(Reader stream) { super(stream) }


    Token nextToken() throws TokenStreamException {
        def tkn = super.nextToken()
        return previous = fix(tkn)
    }

    int currentLevel

    int openLevel

    String block

    Token previous

    boolean process


    protected replace( GroovySourceToken token ) {

        def newToken = new GroovySourceToken(86)
        newToken.with {
            setText('_as')
            setColumn(token.column)
            setColumnLast( token.columnLast )
            setLine(token.line)
            setLineLast( token.lineLast )
        }

        log.trace "* trace: '${token.text}' type: ${token.type} "
        return newToken
    }

    protected Token fix(Token token) {
        switch(token?.type) {

            case LCURLY: // {
                currentLevel ++
                if( process && !openLevel ) {
                    log.trace "Opening 'process' definition [$currentLevel]"
                    openLevel = currentLevel
                }
                break;

            case RCURLY: // }
                if( process && currentLevel == openLevel ) {
                    log.trace "Closing 'process' definition [$currentLevel]"
                    process = false
                    openLevel = 0
                }

                // decrement
                currentLevel --
                break;


            case RPAREN: // )
            case RBRACK: // ]
                currentLevel --
                break

            case LPAREN: // (
            case LBRACK: // [
                currentLevel ++
                break;

            case IDENT:
                if( token.text == 'process' && previous?.type == NLS ) {
                    log.trace "Found 'process' definition [$currentLevel]"
                    process = true
                }
                break

            case LITERAL_as:
                if( block in ['input','output','share'] && process && currentLevel == openLevel ) {
                    token = replace(token as GroovySourceToken)
                }
                break

            case COLON:
                if( previous?.text in ['input','output','share','exec','script']) {
                    block = previous.text
                    log.trace "Entering process block: '$block'"
                }
                break
        }

        return token;
    }

}

class SourcePreProcessor extends ParserPluginFactory {

    Closure<CustomLexer> lexerFactory

    ParserPlugin createParserPlugin() {
        def result = new SourceModifierParserPlugin()
        // inject the lexer factory
        if( lexerFactory )
            result.lexerFactory = lexerFactory
        // return the plugin
        return result
    }
}