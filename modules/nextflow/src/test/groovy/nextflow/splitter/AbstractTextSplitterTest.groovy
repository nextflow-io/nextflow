/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
