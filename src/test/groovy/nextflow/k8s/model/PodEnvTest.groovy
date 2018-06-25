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

package nextflow.k8s.model

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
    class PodEnvTest extends Specification {

    def 'should return env spec' () {
        expect:
        PodEnv.value('ALPHA', 'aaa').toSpec() == [name:'ALPHA', value:'aaa']
    }

    def 'should create env secret spec' () {
        expect:
        PodEnv.secret('ALPHA', 'data/key-1').toSpec() == [
                name: 'ALPHA',
                valueFrom: [secretKeyRef:[name:'data', key:'key-1']]
        ]

        PodEnv.secret('ALPHA', 'data').toSpec() == [
                name: 'ALPHA',
                valueFrom: [secretKeyRef:[name:'data', key:'ALPHA']]
        ]

    }

    def 'should create env config spec' () {
        expect:
        PodEnv.config('ALPHA', 'data/key-1').toSpec() == [
                name: 'ALPHA',
                valueFrom: [configMapKeyRef:[name:'data', key:'key-1']]
        ]
    }

}
