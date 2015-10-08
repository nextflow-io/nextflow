/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Paths

import nextflow.Session
import nextflow.file.FileHolder
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.TokenVar
import org.apache.ignite.logger.slf4j.Slf4jLogger
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgExecutorTest extends Specification {


    def testCopyToTarget() {

        given:
        def targetDir = Files.createTempDirectory('target')
        def scratchDir = Files.createTempDirectory('source')

        def task = [:] as IgBaseTask
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
        task.copyToTargetDir('file1', scratchDir, targetDir)
        then:
        targetDir.resolve('file1').text == 'file 1'

        /*
         * copy multiple files with wildcards
         */
        when:
        task.copyToTargetDir('file_x*', scratchDir, targetDir)
        then:
        targetDir.resolve('file_x_1').text == 'x 1'
        targetDir.resolve('file_x_2').text == 'x 2'
        !targetDir.resolve('file_z_3').exists()

        /*
         * copy a directory maintaining relative paths
         */
        when:
        task.copyToTargetDir('dir_a/dir_b', scratchDir, targetDir)
        then:
        targetDir.resolve('dir_a/dir_b/alpha.txt').text == 'aaa'
        targetDir.resolve('dir_a/dir_b/beta.txt').text == 'bbb'

        cleanup:
        targetDir?.deleteDir()
        scratchDir?.deleteDir()
    }


    def testGgBashTask() {

        given:
        def sessionId = UUID.randomUUID()
        def sourcePath = Files.createTempDirectory('stage-test')
        def sourceFile1 = sourcePath.resolve('file1')
        def sourceFile2 = sourcePath.resolve('file2')

        def targetPath = Files.createTempDirectory('target-path')

        sourceFile1.text = 'Content for file1'
        sourceFile2.text = 'Content for file2'

        def processor = Mock(TaskProcessor)
        processor.getProcessEnvironment() >> [ALPHA: 1, BETA:2 ]
        processor.getSession() >> new Session()
        processor.getConfig() >> [:]

        def binding = new Binding()
        def task = new TaskRun(
                id: 123,
                name: 'TestRun',
                workDir: Paths.get('/some/path'),
                processor: processor,
                config: new TaskConfig(storeDir: targetPath, shell: '/bin/zsh'),
                script: 'echo Hello world!')
        def s1 = new FileInParam(binding, []).bind( new TokenVar('x') )
        task.setInput(s1, [ new FileHolder(sourceFile1), new FileHolder(sourceFile2) ])

        task.setOutput( new FileOutParam(binding, []).bind('x'), 'file1' )
        task.setOutput( new FileOutParam(binding, []).bind('y'), 'file2' )
        task.setOutput( new FileOutParam(binding, []).bind('z'), 'file3' )

        when:
        def ggTask = new IgBashTask(task, sessionId)
        ggTask.log = new Slf4jLogger()

        then:
        ggTask.taskId == 123
        ggTask.name == 'TestRun'
        ggTask.workDir == Paths.get('/some/path')
        ggTask.targetDir == targetPath
        ggTask.inputFiles == [ file1:sourceFile1, file2: sourceFile2 ]
        ggTask.outputFiles == task.getOutputFilesNames()

        ggTask.environment == [ALPHA: 1, BETA:2 ]
        ggTask.shell == ['/bin/zsh']
        ggTask.script == 'echo Hello world!'
        ggTask.sessionId == sessionId

        when:
        ggTask.stage()
        then:
        ggTask.scratchDir != null
        ggTask.scratchDir.resolve('file1').text == sourceFile1.text
        ggTask.scratchDir.resolve('file2').text == sourceFile2.text
        Files.isSymbolicLink(ggTask.scratchDir.resolve('file1'))
        Files.isSymbolicLink(ggTask.scratchDir.resolve('file2'))

        when:
        ggTask.scratchDir.resolve('x').text = 'File x'
        ggTask.scratchDir.resolve('y').text = 'File y'
        ggTask.scratchDir.resolve('z').text = 'File z'
        ggTask.unstage()
        then:
        targetPath.exists()
        targetPath.resolve('x').text == 'File x'
        targetPath.resolve('y').text == 'File y'
        targetPath.resolve('z').text == 'File z'

        cleanup:
        sourcePath?.deleteDir()
        targetPath?.deleteDir()
        ggTask?.scratchDir?.deleteDir()

    }



}
