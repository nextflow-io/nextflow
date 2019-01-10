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

package nextflow.executor

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashTemplateEngineTest extends Specification {

    def 'should strip comments from template' () {

        given:
        def engine = new BashTemplateEngine()
        def template = '''\
            ## comment
            #!/bin/bash
            # NEXTFLOW TASK: foo
            line 1
              ## comment
            ##
            line 2
            line 3
            '''.stripIndent()

        expect:
        engine.render(template, [:]) == '''\
                #!/bin/bash
                # NEXTFLOW TASK: foo
                line 1
                line 2
                line 3
                '''.stripIndent()

    }

}
