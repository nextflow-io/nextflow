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

import groovy.json.JsonOutput
import nextflow.util.CacheHelper

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
        def config = [cid:[store:[location:folder.toString()]]]
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
        folder.resolve(".meta/${hash.toString()}/.data.json").text == JsonOutput.prettyPrint('{"type":"Task","id":100,"name":"foo","hash":"15cd5b07","inputs": null,"annotations":null}')

        cleanup:
        folder?.deleteDir()
    }

    def 'should save task output' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [cid:[store:[location:folder.toString()]]]
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
        def fileHash = CacheHelper.hasher(outFile).hash().toString()
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
        def attrs = Files.readAttributes(outFile, BasicFileAttributes)
        def expectedString = '{"type":"Output",' +
            '"path":"' + outFile.toString() + '",' +
            '"hash":"'+ fileHash + '",' +
            '"source":"cid://15cd5b07",' +
            '"size":'+attrs.size() + ',' +
            '"createdAt":' + attrs.creationTime().toMillis() + ',' +
            '"modifiedAt":'+ attrs.lastModifiedTime().toMillis() + ',' +
            '"annotations":null}'

        and:
        observer.readAttributes(outFile) >> attrs

        when:
        observer.storeTaskOutput(task, outFile)
        then:
        folder.resolve(".meta/${hash}/foo/bar/file.bam/.data.json").text
            == JsonOutput.prettyPrint(expectedString)

        cleanup:
        folder?.deleteDir()
    }

}
