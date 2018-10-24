/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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

}
