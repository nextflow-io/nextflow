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

package nextflow.util

import java.nio.file.Path
import java.nio.file.Paths

import ch.qos.logback.classic.Level
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.cli.CliOptions
import nextflow.extension.OpCall
import nextflow.script.BaseScript
import nextflow.script.ScriptBinding
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LoggerHelperTest extends Specification {


    def 'should parse error line'() {

        given:
        final pwd = System.getProperty("user.dir")
        Map<String,Path> names = [
                'Script_f1bbc0ef': Paths.get('/some/path/main.nf'),
                'Script_1b751fe9': Paths.get('/other/path/module.nf'),
                'Script_12345678': Paths.get("$pwd/foo/script.nf")
        ]

        expect:
        LoggerHelper.getErrorLine(LINE, names) == EXPECTED
        
        where:
        EXPECTED                        | LINE
        null                            | 'at nextflow.script.ScriptRunner.run(ScriptRunner.groovy:289)'
        null                            | 'at nextflow.script.BaseScript.run(BaseScript.groovy:151)'
        ['/some/path/main.nf', '63']    | 'at Script_f1bbc0ef.runScript(Script_f1bbc0ef:63)'
        null                            | 'at Script_1b751fe9$_runScript_closure1.doCall(Script_1b751fe9)'
        ['/other/path/module.nf', '10'] | 'at Script_1b751fe9$_runScript_closure1.doCall(Script_1b751fe9:10)'
        ['foo/script.nf', '55']         | 'at Script_12345678.runScript(Script_12345678:55)'
        null                            | 'at Script_12345678.runScript(Script_xxxxxxxx:63)'
    }


    def 'should create LoggerHelper object' () {

        given:
        def logger = new LoggerHelper(Mock(CliOptions))
        when:
        logger.setDaemon(true)
        then:
        logger.daemon


    }

    def 'should format error message' () {

        given:
        def message =
                """
                startup failed:
                _nf_script_c9a99616: 3: Unknown process block definition: `outpu` @ line 3, column 4.
                stdout() into (A,B,C)
                """
                .stripIndent().leftTrim()


        expect:
        LoggerHelper.formatStartupErrorMessage(message) ==
                """
                Unknown process block definition: `outpu` @ line 3, column 4.
                stdout() into (A,B,C)
                """
                .stripIndent().leftTrim()
    }

    def 'should compare logger levels' () {

        expect:
        assert Level.INFO.levelInt > Level.DEBUG.levelInt
    }


    def 'should check if hash log prefix' () {

        expect:
        LoggerHelper.isHashLogPrefix(STR) == EXPECTED

        where:
        STR             | EXPECTED
        null            | false
        'abc'           | false
        '[xx'           | false
        '[12/abcdef]'   | true
        '[11/xbcdef]'   | false
        '[11/abcdez]'   | false
        '[pp/abcdef]'   | false
    }

    def 'should check hex digit' () {
        expect:
        LoggerHelper.isHex(CHAR as char) == EXPECTED
        where:
        CHAR    | EXPECTED
        '0'     | true
        '5'     | true
        '9'     | true
        'a'     | true
        'f'     | true
        'g'     | false
        'z'     | false
    }

    @Unroll
    def 'should get method error string' () {

        when:
        def ex = new MissingMethodException(METHOD, CLAZZ)
        def result = LoggerHelper.getMissingMethodMessage(ex)
        then:
        result == MSG

        where:
        CLAZZ           | METHOD        | MSG
        ScriptBinding   | 'printx'      | 'Unknown method invocation `printx`'
        BaseScript      | 'printq'      | 'Unknown method invocation `printq` -- Did you mean?\n  print\n  printf'
        DataflowQueue   | 'view'        | 'Unknown method invocation `view` on channel type'
    }

    @Unroll
    def 'should return type message' () {
        expect:
        LoggerHelper.fmtType(TYPE) == MSG

        where:
        TYPE             | MSG
        DataflowQueue    | 'channel type'
        DataflowVariable | 'channel type'
        OpCall           | 'operator type'
        String           | 'String type'
        new DataflowQueue() | 'channel object'
        'abc'           | 'String object'
    }
}
