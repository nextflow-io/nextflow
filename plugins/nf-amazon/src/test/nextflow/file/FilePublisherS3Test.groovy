/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.file

import java.nio.file.Files

import nextflow.Session
import nextflow.cloud.aws.nio.S3Path
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilePublisherS3Test extends Specification {

    def 'should change mode to `copy`' () {

        given:
        def session = Mock(Session) { getConfig() >> [:] }
        def source = Files.createTempDirectory('test').resolve('hello.txt')
        def target = FileHelper.asPath( 's3://bucket/work/hello.txt' )
        def publisher = new FilePublisher(session, [mode: 'symlink'])

        expect:
        publisher.resolveMode(source, target) == FilePublisher.Mode.COPY

        cleanup:
        source?.parent?.deleteDir()
    }

    def 'should tag files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def source = folder.resolve('hello.txt'); source.text = 'Hello'
        and:
        def session = Mock(Session) { getConfig() >> [:] }
        def target = FileHelper.asPath( 's3://bucket/work/hello.txt' )
        def publisher = new FilePublisher(session, [tags: [FOO:'this',BAR:'that']])

        when:
        publisher.applyFileAttributes(source, target)
        then:
        target instanceof S3Path
        (target as S3Path).getTagsList().find{ it.key()=='FOO'}.value() == 'this'
        (target as S3Path).getTagsList().find{ it.key()=='BAR'}.value() == 'that'

        cleanup:
        folder?.deleteDir()
    }

}
