/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.io.FileType
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.formatter.ConfigFormattingVisitor
import nextflow.config.parser.ConfigParserPluginFactory
import nextflow.exception.AbortOperationException
import nextflow.script.formatter.FormattingOptions
import nextflow.script.formatter.ScriptFormattingVisitor
import nextflow.script.parser.ScriptParserPluginFactory
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.ErrorCollector
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Format Nextflow scripts and config files")
class CmdFormat extends CmdBase {

    @Parameter(description = 'List of paths to format')
    List<String> args = []

    @Parameter(names=['-spaces'], description = 'Number of spaces to indent')
    int spaces

    @Parameter(names = ['-tabs'], description = 'Indent with tabs')
    Boolean tabs

    @Parameter(names = ['-harshil-alignment'], description = 'Use Harshil alignment')
    Boolean harhsilAlignment

    private ScriptParser scriptParser

    private ConfigParser configParser

    private FormattingOptions formattingOptions

    @Override
    String getName() { 'format' }

    @Override
    void run() {
        if( !args )
            throw new AbortOperationException("No input files specified")

        if( spaces && tabs )
            throw new AbortOperationException("Cannot specify both `-spaces` and `-tabs`")

        if( !spaces && !tabs )
            spaces = 4

        scriptParser = new ScriptParser()
        configParser = new ConfigParser()
        formattingOptions = new FormattingOptions(spaces, !tabs, harhsilAlignment, false)

        for( final arg : args ) {
            final file = new File(arg)
            if( file.isFile() ) {
                format(file)
                continue
            }

            file.eachFileRecurse(FileType.FILES, this.&format)
        }
    }

    void format(File file) {
        final name = file.getName()
        if( name.endsWith('.nf') )
            formatScript(file)
        else if( name.endsWith('.config') )
            formatConfig(file)
    }

    private void formatScript(File file) {
        final sourceUnit = scriptParser.parse(file)
        if( sourceUnit.getErrorCollector().hasErrors() ) {
            printErrors(sourceUnit)
        }
        else {
            log.debug "Formatting script ${file}"
            final formatter = new ScriptFormattingVisitor(sourceUnit, formattingOptions)
            formatter.visit()
            file.text = formatter.toString()
        }
    }

    private void formatConfig(File file) {
        final sourceUnit = configParser.parse(file)
        if( sourceUnit.getErrorCollector().hasErrors() ) {
            printErrors(sourceUnit)
        }
        else {
            log.debug "Formatting config ${file}"
            final formatter = new ConfigFormattingVisitor(sourceUnit, formattingOptions)
            formatter.visit()
            file.text = formatter.toString()
        }
    }

    private void printErrors(SourceUnit sourceUnit) {
        final errorMessages = sourceUnit.getErrorCollector().getErrors()
        System.err.println()
        System.err.println "Errors found in '${sourceUnit.getName()}':"
        for( final message : errorMessages ) {
            if( message instanceof SyntaxErrorMessage ) {
                final cause = message.getCause()
                System.err.println "Line ${cause.getStartLine()}, column ${cause.getStartColumn()}: ${cause.getOriginalMessage()}"
            }
        }
    }
}

import org.antlr.v4.runtime.RecognitionException
import org.codehaus.groovy.GroovyBugError
import org.codehaus.groovy.control.CompilationFailedException

@CompileStatic
class ScriptParser {

    private final CompilerConfiguration config
    private final GroovyClassLoader classLoader

    ScriptParser() {
        config = getConfig()
        classLoader = new GroovyClassLoader(ClassLoader.getSystemClassLoader().getParent(), config, true)
    }

    SourceUnit parse(File file) {
        final sourceUnit = new SourceUnit(
                file,
                config,
                classLoader,
                new LazyErrorCollector(config))
        try {
            sourceUnit.parse()
            sourceUnit.buildAST()
        }
        catch( RecognitionException e ) {
        }
        catch( CompilationFailedException e ) {
        }
        catch( GroovyBugError | Exception e ) {
            System.err.println "Unexpected exception while compiling source files: ${e}"
        }
        return sourceUnit
    }

    private static CompilerConfiguration getConfig() {
        final config = new CompilerConfiguration()
        config.setPluginFactory(new ScriptParserPluginFactory())
        return config
    }
}

@CompileStatic
class ConfigParser {

    private final CompilerConfiguration config
    private final GroovyClassLoader classLoader

    ConfigParser() {
        config = getConfig()
        classLoader = new GroovyClassLoader(ClassLoader.getSystemClassLoader().getParent(), config, true)
    }

    SourceUnit parse(File file) {
        final sourceUnit = new SourceUnit(
                file,
                config,
                classLoader,
                new LazyErrorCollector(config))
        try {
            sourceUnit.parse()
            sourceUnit.buildAST()
        }
        catch( RecognitionException e ) {
        }
        catch( CompilationFailedException e ) {
        }
        catch( GroovyBugError | Exception e ) {
            System.err.println "Unexpected exception while compiling source files: ${e}"
        }
        return sourceUnit
    }

    private static CompilerConfiguration getConfig() {
        final config = new CompilerConfiguration()
        config.setPluginFactory(new ConfigParserPluginFactory())
        return config
    }
}

@CompileStatic
class LazyErrorCollector extends ErrorCollector {

    LazyErrorCollector(CompilerConfiguration config) {
        super(config)
    }

    @Override
    protected void failIfErrors() {
    }
}
