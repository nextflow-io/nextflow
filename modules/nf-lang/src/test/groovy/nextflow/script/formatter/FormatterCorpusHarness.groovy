package nextflow.script.formatter

import java.nio.file.Files
import java.nio.file.Path

import nextflow.config.control.ConfigParser
import nextflow.config.formatter.ConfigFormattingVisitor
import nextflow.script.control.ScriptParser
import nextflow.script.control.ScriptResolveVisitor
import nextflow.script.dsl.Types
import spock.lang.Specification
import test.TestUtils

/**
 * Scratch harness (not for CI): formats a corpus of realistic scripts and
 * configs under every option combination and asserts that no comment is
 * lost or altered and that formatting is idempotent.
 */
class FormatterCorpusHarness extends Specification {

    static List<FormattingOptions> optionCombos() {
        def combos = []
        for( def spaces : [true, false] )
            for( def harshil : [true, false] )
                for( def mahesh : [true, false] )
                    for( def sort : [true, false] )
                        for( def maxLen : [0, 40, 120] )
                            combos << new FormattingOptions(4, spaces, harshil, mahesh, sort, maxLen)
        return combos
    }

    String formatScript(FormattingOptions options, String contents) {
        def scriptParser = new ScriptParser()
        def source = scriptParser.parse('main.nf', contents)
        new ScriptResolveVisitor(source, scriptParser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, Collections.emptyList()).visit()
        if( TestUtils.hasSyntaxErrors(source) )
            return null // skip files that don't parse
        def formatter = new ScriptFormattingVisitor(source, options, contents)
        formatter.visit()
        return formatter.toString()
    }

    String formatConfig(FormattingOptions options, String contents) {
        def parser = new ConfigParser()
        def source = parser.parse('nextflow.config', contents)
        if( source.getErrorCollector().hasErrors() )
            return null
        def formatter = new ConfigFormattingVisitor(source, options, contents)
        formatter.visit()
        return formatter.toString()
    }

    def 'corpus: comment preservation and idempotence' () {
        given:
        // the test working directory is the module directory
        def root = Path.of(System.getProperty('corpus.root', '../..')).toAbsolutePath().normalize()
        def scripts = []
        def configs = []
        ['tests', 'docs/snippets', 'validation'].each { dir ->
            def base = root.resolve(dir)
            if( !Files.exists(base) )
                return
            base.toFile().eachFileRecurse { f ->
                if( !f.isFile() )
                    return
                if( f.name.endsWith('.nf') )
                    scripts << f
                else if( f.name.endsWith('.config') )
                    configs << f
            }
        }
        def combos = optionCombos()
        int checked = 0
        int skipped = 0
        def failures = []

        when:
        for( def f : scripts + configs ) {
            def isConfig = f.name.endsWith('.config')
            def text = f.text
            def first = null
            try {
                first = isConfig ? formatConfig(combos[0], text) : formatScript(combos[0], text)
            }
            catch( Exception e ) {
                failures << "${f}: [parse/format crashed] ${e}"
                continue
            }
            if( first == null ) {
                skipped++
                continue
            }
            for( def options : combos ) {
                try {
                    def formatted = isConfig ? formatConfig(options, text) : formatScript(options, text)
                    if( formatted == null ) {
                        failures << "${f} ${options}: formatted output is null"
                        continue
                    }
                    if( CommentReattacher.commentTexts(text, isConfig) != CommentReattacher.commentTexts(formatted, isConfig) ) {
                        failures << "${f} ${options}: COMMENTS LOST OR ALTERED"
                        continue
                    }
                    def second = isConfig ? formatConfig(options, formatted) : formatScript(options, formatted)
                    if( second == null ) {
                        failures << "${f} ${options}: REFORMAT FAILED TO PARSE"
                        continue
                    }
                    if( second != formatted ) {
                        failures << "${f} ${options}: NOT IDEMPOTENT"
                        continue
                    }
                    checked++
                }
                catch( Exception e ) {
                    failures << "${f} ${options}: [format crashed] ${e}"
                }
            }
        }
        println "corpus harness: ${checked} file-option checks passed, ${skipped} files skipped (parse errors), ${failures.size()} failures"
        failures.each { println "FAILURE: $it" }

        then:
        failures.isEmpty()
    }

}
