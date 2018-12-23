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

import spock.lang.Specification

import java.nio.file.Files

import nextflow.Session
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
        session.classesDir = Files.createTempDirectory('test')

        when:
        def remote = new RemoteSession(session)
        then:
        remote.localPaths.size() == 3
        remote.getClasspath().size() == 5
        remote.isDeserialized

        cleanup:
        remote?.close()
        path1?.deleteDir()
        path2?.deleteDir()
        session?.classesDir?.deleteDir()
    }

    def 'should zip classpath and lib paths' () {
        given:
        def folder = Files.createTempDirectory('test')
        def cp = folder.resolve('classpath'); cp.mkdir(); Files.createFile(cp.resolve('Main.class'))
        def lib1 = folder.resolve('lib1'); lib1.mkdir(); Files.createFile(lib1.resolve('Lib1.class'))
        def lib2 = folder.resolve('lib2'); lib2.mkdir(); Files.createFile(lib2.resolve('Lib2.class'))

        def session = Mock(Session)

        when:
        def remote = new RemoteSession(session)
        def classpath = remote.getClasspath()

        then:
        1 * session.getClassesDir() >> cp
        1 * session.getLibDir() >> [ lib1, lib2 ]
        classpath.size() == 3
        classpath[0].resolve('Main.class').exists()
        classpath[1].resolve('Lib1.class').exists()
        classpath[2].resolve('Lib2.class').exists()

        when:
        remote.close()
        then:
        !classpath[0].exists()
        !classpath[1].exists()
        !classpath[2].exists()

        cleanup:
        folder?.deleteDir()
    }

}
