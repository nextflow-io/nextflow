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

package nextflow.util
import java.nio.file.Files

import nextflow.Session
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RemoteSessionTest extends Specification {

    def testZipUnzipFolder() {

        given:
        def path = Files.createTempDirectory('zip-test')
        path.resolve('file1').text = 'File 1'
        path.resolve('file2').text = 'File 2'
        path.resolve('dir').mkdir()
        path.resolve('dir/file3').text = 'File 3'
        path.resolve('dir/file4').text = 'File 4'

        when:
        def remote = [:]  as RemoteSession
        def bytes = remote.zip(path)
        then:
        bytes.size()>0

        when:
        def target = remote.unzip(bytes)
        then:
        target.isDirectory()
        target.list().sort() == ['file1','file2','dir'].sort()
        target.resolve('dir').list().sort() == ['file3','file4'].sort()
        target.resolve('file1').text == 'File 1'
        target.resolve('file2').text == 'File 2'
        target.resolve('dir/file3').text == 'File 3'
        target.resolve('dir/file4').text == 'File 4'

        cleanup:
        path?.deleteDir()
        target?.deleteDir()

    }


    def testGetClassPath() {

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

        def session = new Session()
        session.@libDir = [ path1, path2 ]

        when:
        def remote = new RemoteSession(session)
        then:
        remote.getLibDir().size() == 2
        remote.getClasspath().size() == 4
        remote.isLibInitialized

        cleanup:
        remote?.close()
        path1?.deleteDir()
        path2?.deleteDir()
    }

}
