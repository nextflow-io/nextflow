/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Global
import nextflow.Session
import nextflow.processor.TaskBean
import nextflow.file.FilePorter
import spock.lang.Specification

/**
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * based on TesBashBuilderTest by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WrBashBuilderTest extends Specification {

    def 'test bash wrapper' () {
        given:
        def folder = Files.createTempDirectory('test')

        when:
        def bash = new WrBashBuilder([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
        ] as TaskBean)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                        .stripIndent().leftTrim()

        folder.resolve('.command.run').text.contains('nxf_main')

        cleanup:
        folder.deleteDir()
    }

    def 'test bash wrapper with outputs in S3' () {
        given:
        def folder = 's3://bucket/work' as Path
        def bash = new WrBashBuilder([
                name: 'Hello 2',
                workDir: folder,
                targetDir: folder,
                script: 'echo Hello world!',
                outputFiles: ['test.bam','test.bai'],
        ] as TaskBean)

        when:
        def binding = bash.makeBinding()

        then:
        binding.touch_file == 'touch .mnt/bucket/work/.command.begin'
        binding.unstage_outputs == '''\
                  cp -fRL test.bam .mnt/bucket/work/ || true
                  cp -fRL test.bai .mnt/bucket/work/ || true
                  '''.stripIndent().rightTrim()
    }

    // *** I can't figure out how to get this test to work without it trying to
    //     actually upload the input files to S3
    //
    // def setupSpec() {
    //     new Session()
    // }
    // def 'test bash wrapper with inputs in S3' () {
    //     given:
    //     def s = (Session)Global.session
    //     def porter = Spy(new FilePorter(s))
    //     s.filePorter = porter
    //     def folder = 's3://bucket/work' as Path
    //     def bash = new WrBashBuilder([
    //             name: 'Hello 3',
    //             workDir: '/local/work' as Path,
    //             targetDir: folder,
    //             script: 'echo Hello world!',
    //             inputFiles: ['sample_1.fq': 's3://input/data/sample_1.fq' as Path, 'sample_2.fq': '/input/data/sample_2.fq' as Path],
    //     ] as TaskBean)

    //     when:
    //     def binding = bash.makeBinding()
    //     porter.stageForeignFiles(_) >> { def input -> return input }

    //     then:
    //     binding.stage_inputs == '''\
    //             rm -f sample_1.fq
    //             rm -f sample_2.fq
    //             ln -s .inputs/input/data/sample_1.fq sample_1.fq
    //             ln -s /input/data/sample_2.fq sample_2.fq
    //             '''.stripIndent().rightTrim()
    // }

}
