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
 *
 */

package nextflow.data.cid

import java.nio.file.Files
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileTime
import java.time.Instant

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CidObserverTest extends Specification {

    def 'should save task run' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        def session = Mock(Session) { getConfig()>>config }
        def observer = new CidObserver()
        observer.onFlowCreate(session)
        and:
        def hash = HashCode.fromInt(123456789)
        and:
        def task = Mock(TaskRun) {
            getId() >> TaskId.of(100)
            getName() >> 'foo'
            getHash() >> hash
        }
        when:
        observer.storeTaskRun(task)
        then:
        folder.resolve(hash.toString()).text == '{"id":100,"name":"foo","hash":"15cd5b07","annotations":null}'

        cleanup:
        folder?.deleteDir()
    }

    def 'should save task output' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        def session = Mock(Session) { getConfig()>>config }
        def observer = Spy(new CidObserver())
        observer.onFlowCreate(session)
        and:
        def workDir = folder.resolve('12/34567890')
        Files.createDirectories(workDir)
        and:
        def outFile = workDir.resolve('foo/bar/file.bam')
        Files.createDirectories(outFile.parent)
        outFile.text = 'some data'
        and:
        def hash = HashCode.fromInt(123456789)
        and:
        def task = Mock(TaskRun) {
            getId() >> TaskId.of(100)
            getName() >> 'foo'
            getHash() >> hash
            getWorkDir() >> workDir
        }
        and:
        def ts1 = Instant.ofEpochMilli(1737914400)
        def ts2 = Instant.ofEpochMilli(1737914500)
        def attrs = Mock(BasicFileAttributes) {
            size() >> 100
            creationTime() >> FileTime.from(ts1)
            lastModifiedTime() >> FileTime.from(ts2)
        }
        and:
        observer.readAttributes(outFile) >> attrs

        when:
        observer.storeTaskOutput(task, outFile)
        then:
        folder.resolve("${hash}/foo/bar/file.bam").text
            == '{"uri":"cid://15cd5b07/foo/bar/file.bam","size":100,"createdAt":1737914400,"modifiedAt":1737914500,"annotations":null}'

        cleanup:
        folder?.deleteDir()
    }

}
