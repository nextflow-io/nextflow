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
class PodVolumeClaimTest extends Specification {


    def 'should create a pod volume claim' () {
        when:
        def vol1 = new PodVolumeClaim('foo', '/bar')

        then:
        vol1.claimName == 'foo'
        vol1.mountPath == '/bar'

        when:
        def vol2 = new PodVolumeClaim(volumeClaim: 'alpha', mountPath: '/gamma')
        then:
        vol2.claimName == 'alpha'
        vol2.mountPath == '/gamma'

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
