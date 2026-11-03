/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.script.control

import nextflow.script.ast.ProcessNodeV2
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 * @see nextflow.script.control.ImplicitStagers
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ImplicitStagersTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    def 'should generate implicit stagers for process inputs'() {
        when:
        final stagers = stagerTargets("""
            nextflow.enable.types = true

            record Vcf {
                id: String
                vcf: Path
                tbi: Path?
                extra: List<Path>
            }

            record Cohort {
                id: String
                vcfs: List<Vcf>
            }

            process CHECK {
                input:
                ${declaration}

                script:
                'true'
            }
            """)

        then:
        stagers == expected

        where:
        declaration                             || expected
        'x: Integer'                            || []
        'vcf: Path'                             || ['vcf']
        'vcf: Path?'                            || ['vcf']
        'vcfs: List<Path>'                      || ['vcfs']
        'vcfs: Set<Path>'                       || ['vcfs']
        'record(id: String, vcf: Path)'         || ['vcf']
        'vcf: Vcf'                              || ['vcf?.vcf', 'vcf?.tbi', 'vcf?.extra']
        'vcf: Vcf?'                             || ['vcf?.vcf', 'vcf?.tbi', 'vcf?.extra']
        'vcfs: List<Vcf>'                       || [
            'vcfs*.vcf?.findAll({ java.lang.Object p -> ... })',
            'vcfs*.tbi?.findAll({ java.lang.Object p -> ... })',
            'vcfs*.extra?.collectMany({ java.lang.Object e -> ... })'
        ]
        'record(id: String, vcfs: List<Vcf>)'  || [
            'vcfs*.vcf?.findAll({ java.lang.Object p -> ... })',
            'vcfs*.tbi?.findAll({ java.lang.Object p -> ... })',
            'vcfs*.extra?.collectMany({ java.lang.Object e -> ... })'
        ]
        'cohorts: List<Cohort>'                 || [
            'cohorts*.vcfs?.collectMany({ java.lang.Object e -> ... })*.vcf?.findAll({ java.lang.Object p -> ... })',
            'cohorts*.vcfs?.collectMany({ java.lang.Object e -> ... })*.tbi?.findAll({ java.lang.Object p -> ... })',
            'cohorts*.vcfs?.collectMany({ java.lang.Object e -> ... })*.extra?.collectMany({ java.lang.Object e -> ... })'
        ]
    }

    private List<String> stagerTargets(String contents) {
        scriptParser.compiler().getSources().clear()
        final source = scriptParser.parse('main.nf', contents.stripIndent())
        scriptParser.analyze()
        final errors = TestUtils.getErrors(source)
        assert errors.isEmpty() : errors.collect { it.getOriginalMessage() + ' @' + it.getStartLine() }.join('; ')

        // the stagers are the second statement of the process body
        final process = (ProcessNodeV2) source.getAST().getProcesses().first()
        final result = new ProcessToGroovyVisitorV2(source).transform(process)
        final processCall = (MethodCallExpression) ((ExpressionStatement) result).getExpression()
        final body = (BlockStatement) ((ClosureExpression) processCall.getArguments().getExpression(1)).getCode()
        final stagers = (BlockStatement) body.getStatements().get(1)
        return stagers.getStatements().collect { st ->
            final call = (MethodCallExpression) ((ExpressionStatement) st).getExpression()
            assert call.getMethodAsString() == 'stageAs'
            final closure = (ClosureExpression) call.getArguments().getExpression(0)
            return closureBody(closure).getText()
        }
    }

    private static closureBody(ClosureExpression closure) {
        final code = closure.getCode()
        final st = code instanceof BlockStatement ? code.getStatements().get(0) : code
        return ((ExpressionStatement) st).getExpression()
    }

}
