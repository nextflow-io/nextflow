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

package nextflow.ga4gh.tes.executor

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TesFileCopyStrategyTest extends Specification {

    def 'should return task env' () {

        given:
        def strategy = new TesFileCopyStrategy()

        when:
        def script = strategy.getEnvScript([ALPHA:'xx', BETA:'yy'])
        then:
        script == '''\
            export ALPHA="xx"
            export BETA="yy"
            '''.stripIndent()

        when:
        script = strategy.getEnvScript([ALPHA:'xx', BETA:'yy'], 'foo')
        then:
        script == '''\
            foo() {
            cat << EOF
            export ALPHA="xx"
            export BETA="yy"
            EOF
            }
            '''.stripIndent()
    }
}
