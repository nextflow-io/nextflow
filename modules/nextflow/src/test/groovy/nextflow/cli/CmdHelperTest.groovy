/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cli

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdHelperTest extends Specification {

    def 'should normalise filter string' () {

        expect:
        CmdHelper.fixEqualsOp(null) == null
        CmdHelper.fixEqualsOp('') == ''

        CmdHelper.fixEqualsOp('cpu=1')  == 'cpu==1'
        CmdHelper.fixEqualsOp('cpu==1')  == 'cpu==1'
        CmdHelper.fixEqualsOp('cpu>=1')  == 'cpu>=1'
        CmdHelper.fixEqualsOp('cpu<=1')  == 'cpu<=1'
        CmdHelper.fixEqualsOp('cpu!=1')  == 'cpu!=1'

        CmdHelper.fixEqualsOp('cpu = 1 && mem >= 10')  == 'cpu == 1 && mem >= 10'
        CmdHelper.fixEqualsOp('cpu == 1 && mem >= 10')  == 'cpu == 1 && mem >= 10'

        CmdHelper.fixEqualsOp('=')  == '='
        CmdHelper.fixEqualsOp('==')  == '=='
        CmdHelper.fixEqualsOp('x=')  == 'x='
        CmdHelper.fixEqualsOp('=x')  == '=x'

    }

}
