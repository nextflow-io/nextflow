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
class PodMountSecretTest extends Specification {

    def 'should create mount for configmap' () {

        when:
        def opt = new PodMountSecret(mountPath: '/etc/some/name', secret: 'here' )
        then:
        opt.mountPath == '/etc/some/name'
        opt.fileName == null
        opt.secretName == 'here'
        opt.secretKey == null

        when:
        opt = new PodMountSecret(mountPath: '/etc/some/name', secret: 'here/there.txt' )
        then:
        opt.mountPath == '/etc/some'
        opt.fileName == 'name'
        opt.secretName == 'here'
        opt.secretKey == 'there.txt'

    }
}
