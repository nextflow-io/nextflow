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

package nextflow.splitter

import java.nio.file.Paths

import nextflow.Session
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractTextSplitterTest extends Specification {

    def testGetCollectorBaseFile() {

        given:
        def folder = TestHelper.createInMemTempDir()
        def session = new Session()
        session.workDir = folder
        def result
        def splitter

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.options(file: folder, by: 2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path == folder.resolve('chunk')
        result.hash != null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.options(file: folder, by: 2, elem:2)
        splitter.multiSplit = true
        result = splitter.getCollectorBaseFile()
        then:
        result.path == folder.resolve('chunk_2')
        result.hash != null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.options(file: Paths.get('/some/file.txt'), by: 2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path == Paths.get('/some/file.txt')
        result.hash != null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.options(file: true, by: 2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path.name == 'chunk'
        result.path.toString().startsWith( folder.toString() )
        result.hash == null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.options(file: 'chunk_name', by:2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path.name == 'chunk_name'
        result.path.toString().startsWith( folder.toString() )
        result.hash == null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.options(file: 'chunk_name', by:2, elem:3)
        splitter.multiSplit = true
        result = splitter.getCollectorBaseFile()
        then:
        result.path.name == 'chunk_name_3'
        result.path.toString().startsWith( folder.toString() )
        result.hash == null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.sourceFile = Paths.get('/some/file.txt')
        splitter.options(file: true, by: 2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path.name == 'file.txt'
        result.path.toString().startsWith( folder.toString() )
        result.hash == null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.sourceFile = Paths.get('/some/file.fasta.gz')
        splitter.options(file: true, by: 2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path.name == 'file.fasta'
        result.path.toString().startsWith( folder.toString() )
        result.hash == null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.sourceFile = Paths.get('/some/file.fa.gz')
        splitter.options(file: 'my_file_name', by: 2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path.name == 'my_file_name'
        result.path.toString().startsWith( folder.toString() )
        result.hash == null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.sourceFile = Paths.get('/some/file.fa.gz')
        splitter.options(file: folder, by: 2)
        result = splitter.getCollectorBaseFile()
        then:
        result.path.name == 'file.fa'
        result.path.toString().startsWith( folder.toString() )
        result.hash != null

        when:
        splitter = [:] as AbstractTextSplitter
        splitter.sourceFile = Paths.get('/some/file.fasta.gz')
        splitter.options(file: true)
        result = splitter.createCollector()
        then:
        result != null
    }

    def 'should change cache hash when using a different session' () {

        given:
        def folder = TestHelper.createInMemTempDir()
        def result
        def splitter = [:] as AbstractTextSplitter
        splitter.options(file: folder, by: 2)

        when:
        new Session()
        result = splitter.getCollectorBaseFile()
        then:
        result.path == folder.resolve('chunk')
        result.hash != null

        when:
        // copy current hash
        def hash = result.hash
        // invoke with using the SAME session
        result = splitter.getCollectorBaseFile()
        then:
        result.path == folder.resolve('chunk')
        result.hash != null
        // hash is the SAME
        result.hash.toString() == hash.toString()

        when:
        // create NEW session
        new Session()
        // invoke with using DIFFERENT session
        result = splitter.getCollectorBaseFile()
        then:
        result.path == folder.resolve('chunk')
        result.hash != null
        // hashes are DIFFERENT
        result.hash.toString() != hash.toString()

    }

}
