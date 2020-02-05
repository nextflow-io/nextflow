/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.ast

import java.nio.file.Paths

import nextflow.processor.TaskPath
import nextflow.util.BlankSeparatedList
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskCmdXformVisitorTest extends Specification {

    def 'should replace special values in task command' () {

        given:
        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskCmdXform))
        and:
        def context = [
                X: 'foo',
                Y: new TaskPath(Paths.get('f o o.txt')),
                Z: new BlankSeparatedList(1, 2, 3)
                ]
        def shell = new GroovyShell(new Binding(context), config)

        when:
        def result = shell.evaluate('''
            return "$X, $Y, $Z"
            ''')
        
        then:
        result == 'foo, f\\ o\\ o.txt, 1 2 3'
    }
}
