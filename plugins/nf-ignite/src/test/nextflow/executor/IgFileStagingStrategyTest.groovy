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
