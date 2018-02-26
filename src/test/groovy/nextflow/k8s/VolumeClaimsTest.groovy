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

package nextflow.k8s

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class VolumeClaimsTest extends Specification {

    def 'should find volume by path' () {

        when:
        def volumes = new VolumeClaims( [ foo: [mountPath: '/data/work'], bar: [mountPath: '/other/path'] ] )
        then:
        volumes.findVolumeByPath('/data/work') == 'foo'
        volumes.findVolumeByPath('/other/path') == 'bar'
        volumes.findVolumeByPath('/data') == null
    }

    def 'should collect mount paths' () {

        when:
        def volumes = new VolumeClaims( [ foo: [mountPath: '/data/work'], bar: [mountPath: '/other/path'] ] )
        then:
        volumes.getMountPaths() == ['/data/work', '/other/path']

    }

    def 'should sanitize paths' () {
        given:
        def VOLS = [
                vol1: [mountPath: '/data/work//'],
                vol2: [mountPath: '//'],
                vol3: [mountPath: '/data']]

        when:
        new VolumeClaims(VOLS)
        then:
        VOLS.vol1.mountPath == '/data/work'
        VOLS.vol2.mountPath == '/'
        VOLS.vol3.mountPath == '/data'
    }

    def 'should add a new volume' () {

        given:
        def volClaims = new VolumeClaims()

        when:
        volClaims.add('vol1', '/some/path/')
        then:
        volClaims.get('vol1').mountPath == '/some/path'

        when:
        volClaims.add('vol2', '/other/path')
        then:
        volClaims.get('vol2').mountPath == '/other/path'
        volClaims.getFirstMount() == '/some/path'

        when:
        volClaims.add('vol3:/this/that')
        then:
        volClaims.get('vol3').mountPath == '/this/that'
    }

    def 'should add a volumes only if not exists' () {
        given:
        def volClaims = new VolumeClaims()
        volClaims.add('vol1', '/some/path/')

        when:
        volClaims.addAllSkipExisting( [vol1: [mountPath: '/foo'], vol2: [mountPath:'/bar']] )

        then:
        volClaims.size() == 2 
        volClaims.vol1.mountPath == '/some/path'
        volClaims.vol2.mountPath == '/bar'
    }
}
