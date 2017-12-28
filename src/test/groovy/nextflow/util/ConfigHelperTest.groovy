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

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigHelperTest extends Specification {

    @Unroll
    def "should parse string value: #str" () {

        expect:
        ConfigHelper.parseValue(str) == value

        where:
        str         | value
        'hola'      | 'hola'
        '1'         | 1
        "${Long.MAX_VALUE}" | Long.MAX_VALUE
        'True'      | true
        'False'     | false
        "10.2"      | 10.2
        '5sec'      | Duration.of('5sec')
        'live_in_3d'| 'live_in_3d'

    }

    def testResolveClasspaths() {

        given:
        def path1 = Files.createTempDirectory('path1')
        path1.resolve('file1').text = 'File 1'
        path1.resolve('file2.jar').text = 'File 2'
        path1.resolve('dir').mkdir()
        path1.resolve('dir/file3').text = 'File 3'
        path1.resolve('dir/file4').text = 'File 4'

        def path2 = Files.createTempDirectory('path2')
        path2.resolve('file5').text = 'File 5'
        path2.resolve('file6.jar').text = 'File 6'

        def path3 = Paths.get('/some/file')

        when:
        def list = ConfigHelper.resolveClassPaths([path1, path2, path3])
        then:
        list.size() == 4
        list.contains( path1 )
        list.contains( path1.resolve('file2.jar') )
        list.contains( path2 )
        list.contains( path2.resolve('file6.jar') )

        cleanup:
        path1?.deleteDir()
        path2?.deleteDir()

    }


}
