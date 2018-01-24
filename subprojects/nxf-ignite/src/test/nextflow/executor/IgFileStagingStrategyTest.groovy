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

package nextflow.executor
import java.nio.file.Files

import nextflow.processor.TaskBean
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgFileStagingStrategyTest extends Specification {


    def 'should copy files to the target directory'() {

        given:
        def targetDir = Files.createTempDirectory('target')
        def scratchDir = Files.createTempDirectory('source')

        def strategy = [:] as IgFileStagingStrategy
        scratchDir.resolve('file1').text = 'file 1'
        scratchDir.resolve('file_x_1').text = 'x 1'
        scratchDir.resolve('file_x_2').text = 'x 2'
        scratchDir.resolve('file_z_3').text = 'z 3'
        scratchDir.resolve('dir_a').resolve('dir_b').mkdirs()
        scratchDir.resolve('dir_a').resolve('dir_b').resolve('alpha.txt').text = 'aaa'
        scratchDir.resolve('dir_a').resolve('dir_b').resolve('beta.txt').text = 'bbb'

        /*
         * copy a single file
         */
        when:
        strategy.copyToTargetDir('file1', scratchDir, targetDir)
        then:
        targetDir.resolve('file1').text == 'file 1'

        /*
         * copy multiple files with wildcards
         */
        when:
        strategy.copyToTargetDir('file_x*', scratchDir, targetDir)
        then:
        targetDir.resolve('file_x_1').text == 'x 1'
        targetDir.resolve('file_x_2').text == 'x 2'
        !targetDir.resolve('file_z_3').exists()

        /*
         * copy a directory maintaining relative paths
         */
        when:
        strategy.copyToTargetDir('dir_a/dir_b', scratchDir, targetDir)
        then:
        targetDir.resolve('dir_a/dir_b/alpha.txt').text == 'aaa'
        targetDir.resolve('dir_a/dir_b/beta.txt').text == 'bbb'

        cleanup:
        targetDir?.deleteDir()
        scratchDir?.deleteDir()
    }


    def 'should stage input files and unstage output files'() {

        given:
        def sessionId = UUID.randomUUID()
        def sourcePath = Files.createTempDirectory('stage-test')
        def sourceFile1 = sourcePath.resolve('file1')
        def sourceFile2 = sourcePath.resolve('file2')

        def targetPath = Files.createTempDirectory('target-path')

        sourceFile1.text = 'Content for file1'
        sourceFile2.text = 'Content for file2'

        def task = new TaskBean(
                inputFiles: [ file1:sourceFile1, file2: sourceFile2 ],
                outputFiles: ['x','y','z'],
                targetDir: targetPath
                )

        when:
        def ggTask = new IgFileStagingStrategy(sessionId: sessionId, task: task)
        then:
        ggTask.sessionId == sessionId

        when:
        ggTask.stage()
        then:
        ggTask.localWorkDir != null
        ggTask.localWorkDir.resolve('file1').text == sourceFile1.text
        ggTask.localWorkDir.resolve('file2').text == sourceFile2.text
        Files.isSymbolicLink(ggTask.localWorkDir.resolve('file1'))
        Files.isSymbolicLink(ggTask.localWorkDir.resolve('file2'))

        when:
        ggTask.localWorkDir.resolve('x').text = 'File x'
        ggTask.localWorkDir.resolve('y').text = 'File y'
        ggTask.localWorkDir.resolve('z').text = 'File z'
        ggTask.unstage()
        then:
        targetPath.exists()
        targetPath.resolve('x').text == 'File x'
        targetPath.resolve('y').text == 'File y'
        targetPath.resolve('z').text == 'File z'

        cleanup:
        sourcePath?.deleteDir()
        targetPath?.deleteDir()
        ggTask?.localWorkDir?.deleteDir()

    }


}
