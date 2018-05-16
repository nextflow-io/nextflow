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

package nextflow.k8s.client

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sResponseExceptionTest extends Specification {

    def 'should create response from valid json' () {

        when:
        def resp = new K8sResponseException(
                'Request /this/that failed',
                new K8sResponseJson('{"foo":"one","bar":"two"}'))
        then:
        resp.getMessage() == '''
                            Request /this/that failed
                            
                              {
                                  "foo": "one",
                                  "bar": "two"
                              }
                            '''.stripIndent().leftTrim()
    }


    def 'should create response from error message' () {

        when:
        def resp = new K8sResponseException(
                'Request /this/that failed',
                new K8sResponseJson('Oops.. it crashed badly'))
        then:
        resp.getMessage() == 'Request /this/that failed -- Oops.. it crashed badly'
    }
}
