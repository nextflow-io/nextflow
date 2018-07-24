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
