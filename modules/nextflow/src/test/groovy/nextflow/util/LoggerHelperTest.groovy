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

import spock.lang.Specification

import ch.qos.logback.classic.Level
import nextflow.cli.CliOptions
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LoggerHelperTest extends Specification {


    def 'should parse error line'() {

        expect:
        LoggerHelper.getErrorLine('at pfam3d.run(pfam3d.nf:189) ~[na:na]') == ['pfam3d.nf','189']
        LoggerHelper.getErrorLine('at pfam3d.run(JavaClass.java:189) ~[na:na]') == null
        LoggerHelper.getErrorLine('at pfam3d.run(pfam3d.nf:189) ~[na:na]','pfam3d.nf') == ['pfam3d.nf','189']
        LoggerHelper.getErrorLine('at pfam3d.run(pfam3d.nf:189) ~[na:na]','hola') == null
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
}
