/*
 * Copyright 2013-2024, Seqera Labs
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
import java.nio.file.Path

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ClassLoaderFactoryTest extends Specification {

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

        def path3 = Path.of('/some/file')

        when:
        def list = ClassLoaderFactory.resolveClassPaths([path1, path2, path3])
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
