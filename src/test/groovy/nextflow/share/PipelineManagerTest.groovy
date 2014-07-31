/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.share
import java.nio.file.Files

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PipelineManagerTest extends Specification {

    def testList() {

        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()

        def manager = new PipelineManager().setRoot(folder.toFile())

        when:
        def list = manager.list()
        then:
        list.sort() == ['cbcrg/pipe1','cbcrg/pipe2','ncbi/blast']

        cleanup:
        folder?.deleteDir()

    }


    def testPull() {

        given:
        def folder = Files.createTempDirectory('test')
        def manager = new PipelineManager().setRoot(folder.toFile()).setPipeline('pditommaso/awesome-pipeline')

        when:
        manager.pull()
        then:
        folder.resolve('pditommaso/awesome-pipeline/.git').isDirectory()

        when:
        manager.pull()
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()

    }

    def testClone() {

        given:
        def dir = Files.createTempDirectory('test')
        def manager = new PipelineManager('pditommaso/awesome-pipeline')

        when:
        manager.clone(dir.toFile())

        then:
        dir.resolve('README.md').exists()
        dir.resolve('.git').isDirectory()

        cleanup:
        dir?.deleteDir()

    }

}
