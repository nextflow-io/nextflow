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

package nextflow.file.igfs

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgPathTest extends Specification{


    def testCreate() {
        given:
        def fs = Mock(IgFileSystem)

        when:
        def path = new IgPath(fs, 'some/path/file.name')
        then:
        path.toString() == 'some/path/file.name'
        path.getParent().toString() == 'some/path'
        path.getFileName().toString() == 'file.name'
        path.getName(0).toString() == 'some'
        path.getName(1).toString() == 'path'
        path.getName(2).toString() == 'file.name'
        path.toAbsolutePath().toString() == '/some/path/file.name'
    }

    def testIterator() {

        when:
        def path = new IgPath(Mock(IgFileSystem), '/some/path/file.name')
        then:
        path.iterator().collect { it.toString() } == ['some','path','file.name']

    }

    def testEqualsAndHashCode() {

        given:
        def fs1 = Mock(IgFileSystem);
        def fs2 = Mock(IgFileSystem);

        when:
        def path1 = new IgPath(fs1, '/some/path/file_x')
        def path2 = new IgPath(fs1, '/some/path/file_z')
        def path3 = new IgPath(fs1, '/some/path/file_x')
        def path4 = new IgPath(fs2, '/some/path/file_x')

        then:
        !path1.equals(path2)
        path1.equals(path3)
        !path1.equals(path4)

        path1.hashCode() != path2.hashCode()
        path1.hashCode() == path3.hashCode()


    }



}
