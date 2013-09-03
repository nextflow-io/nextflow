/*
 * Copyright (c) 2012, the authors.
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

package nextflow

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowTest extends Specification {


    def testFile() {

        expect:
        Nextflow.file('file.log') == new File('file.log').canonicalFile
        Nextflow.file('relative/file.test') == new File( new File('.').canonicalFile, 'relative/file.test')
        Nextflow.file('/user/home/file.log') == new File('/user/home/file.log')
        Nextflow.file('~') == new File( System.getProperty('user.home') )
        Nextflow.file('~/file.test') == new File( System.getProperty('user.home'), 'file.test' )
        Nextflow.file('~file.name') == new File('~file.name').canonicalFile

    }

    def testFileWithWildcards() {

        setup:
        new File('hola1').text = 'abc'
        new File('helo2').text = 'abc'
        new File('hello3').text = 'abc'

        expect:
        Nextflow.file('ciao*') == []
        Nextflow.file('hel*').sort() == [ new File('helo2').canonicalFile, new File('hello3').canonicalFile ].sort()
        Nextflow.file('hol??') == [ new File('hola1').canonicalFile  ]


        cleanup:
        new File('hola1').delete()
        new File('helo2').delete()
        new File('hello3').delete()

    }
}
