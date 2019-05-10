/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import spock.lang.Specification

/**
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 */
class WrFileCopyStrategyTest extends Specification {

    // *** write tests

    // def 'should return task env' () {

    //     given:
    //     def strategy = new TesFileCopyStrategy()

    //     when:
    //     def script = strategy.getEnvScript([ALPHA:'xx', BETA:'yy'],false)
    //     then:
    //     script == '''\
    //         export ALPHA="xx"
    //         export BETA="yy"
    //         '''.stripIndent()

    //     when:
    //     script = strategy.getEnvScript([ALPHA:'xx', BETA:'yy'], true)
    //     then:
    //     thrown(UnsupportedOperationException)
    // }
}
