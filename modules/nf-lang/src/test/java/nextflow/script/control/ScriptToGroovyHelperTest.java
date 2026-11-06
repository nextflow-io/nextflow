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
package nextflow.script.control;

import nextflow.script.ast.ProcessNodeV1;
import nextflow.script.ast.ScriptNode;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.ConstructorCallExpression;
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.control.SourceUnit;
import org.junit.jupiter.api.Test;

import static nextflow.test.TestUtils.getErrors;
import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies the helpers used to lower scripts to Groovy.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptToGroovyHelperTest {

    private static SourceUnit parse(String contents) {
        var source = new ScriptParser().parse("main.nf", contents);
        assertThat(getErrors(source)).isEmpty();
        return source;
    }

    @Test
    void collectParamsAndTaskExtReferencesFromProcessScript() {
        var source = parse("""
            process hello {
                input:
                val x

                script:
                \"""
                echo ${params.x} ${task.ext.args} ${x}
                \"""
            }
            """);
        var process = (ProcessNodeV1) ((ScriptNode) source.getAST()).getProcesses().get(0);

        var refs = (ListExpression) new ScriptToGroovyHelper(source).getVariableRefs(process.exec);
        var names = refs.getExpressions().stream()
            .map(ref -> ((TupleExpression) ((ConstructorCallExpression) ref).getArguments()).getExpression(0))
            .map(name -> ((ConstantExpression) name).getValue())
            .toList();

        assertThat(names).containsExactlyInAnyOrder("params.x", "task.ext.args");
    }

    @Test
    void getSourceTextOfProcessBody() {
        var source = parse("""
            process hello {
                script:
                \"""
                echo 'hello!'
                \"""
            }
            """);
        var process = (ProcessNodeV1) ((ScriptNode) source.getAST()).getProcesses().get(0);

        assertThat(new ScriptToGroovyHelper(source).getSourceText(process.exec))
            .isEqualTo("    \"\"\"\n    echo 'hello!'\n    \"\"\"\n");
    }

}
