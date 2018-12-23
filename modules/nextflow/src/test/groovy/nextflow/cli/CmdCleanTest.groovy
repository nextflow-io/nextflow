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

package nextflow.cli

import spock.lang.Specification
import java.nio.file.Files

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 */
class CmdCleanTest extends Specification {

    def 'empty folder' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cleaner = new CmdClean()

        when:
        def result = cleaner.deleteFolder(folder,false)


        then:
        assert(result == true)
        assert(!folder.exists())

        cleanup:
        def exists = folder?.exists()
        if (exists)
            folder?.deleteDir()

    }

    def 'empty folder, keep logs' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cleaner = new CmdClean()

        when:
        def result = cleaner.deleteFolder(folder,true)

        then:
        assert(result == true)
        assert(folder.exists())

        cleanup:
        def exists = folder?.exists()
        if (exists)
            folder?.deleteDir()

    }


    def 'non-empty folder (.command.file name) keep logs' () {

        given:
        def folder = Files.createTempDirectory('test')
        def file = Files.createFile(folder.resolve('.command.test'))
        def cleaner = new CmdClean()

        when:
        def result = cleaner.deleteFolder(file,true)

        then:
        assert(result == true)
        assert(file.getName().equals('.command.test'))

        cleanup:
        def fileExists = file?.exists()
        if (fileExists)
            file.delete()
        def folderExists = folder?.exists()
        if (folderExists)
            folder?.deleteDir()
    }





}
