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

package nextflow.processor

import java.nio.file.FileSystems
import java.nio.file.Files

import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.nio.S3Path
import nextflow.file.FileHelper
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PublishDirS3Test extends Specification {

    def 'should change mode to `copy`' () {

        given:
        def processor = [:] as TaskProcessor
        processor.name = 'foo'

        def targetDir = FileHelper.asPath( 's3://bucket/work' )
        def publisher = new PublishDir(mode:'symlink', path: targetDir, sourceFileSystem: FileSystems.default)

        when:
        publisher.validatePublishMode()
        then:
        publisher.mode == PublishDir.Mode.COPY
    }

    def 'should tag files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def source = folder.resolve('hello.txt'); source.text = 'Hello'

        and:
        def processor = [:] as TaskProcessor
        processor.name = 'foo'
        and:
        def targetDir = FileHelper.asPath( 's3://bucket/work' )
        def publisher = new PublishDir(tags: [FOO:'this',BAR:'that'], path: targetDir, sourceFileSystem: FileSystems.default)
        def spy = Spy(publisher)

        when:
        spy.apply1(source, true)
        then:
        1 * spy.safeProcessFile(source, _) >> { sourceFile, s3File ->
            assert s3File instanceof S3Path
            assert (s3File as S3Path).getTagsList().find{ it.getKey()=='FOO'}.value == 'this'
            assert (s3File as S3Path).getTagsList().find{ it.getKey()=='BAR'}.value == 'that'
        }

        cleanup:
        folder?.deleteDir()
    }

}
