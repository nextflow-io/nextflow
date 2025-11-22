/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.k8s.model

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PodVolumeClaimTest extends Specification {


    def 'should create a pod volume claim' () {
        when:
        def vol1 = new PodVolumeClaim('foo', '/bar')

        then:
        vol1.claimName == 'foo'
        vol1.mountPath == '/bar'
        vol1.readOnly == false

        when:
        def vol2 = new PodVolumeClaim(volumeClaim: 'alpha', mountPath: '/gamma')
        then:
        vol2.claimName == 'alpha'
        vol2.mountPath == '/gamma'
        vol2.readOnly == false

        when:
        def vol3 = new PodVolumeClaim('aaa', '/bbb', null, true)

        then:
        vol3.claimName == 'aaa'
        vol3.mountPath == '/bbb'
        vol3.readOnly == true

        when:
        def vol4 = new PodVolumeClaim(volumeClaim: 'ccc', mountPath: '/ddd', readOnly: true)

        then:
        vol4.claimName == 'ccc'
        vol4.mountPath == '/ddd'
        vol4.readOnly == true

    }

    def 'should sanitize paths' () {

        expect :
        new PodVolumeClaim('foo','/data/work//').mountPath == '/data/work'
        new PodVolumeClaim('foo','//').mountPath == '/'
        new PodVolumeClaim('foo','/data').mountPath == '/data'

        when:
        new PodVolumeClaim('foo','data')
        then:
        thrown(IllegalArgumentException)
    }

}
