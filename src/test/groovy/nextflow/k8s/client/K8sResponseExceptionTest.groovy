/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
