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

package nextflow.processor

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.ScriptRuntimeException
import nextflow.file.FilePublisher
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PublishDirTest extends Specification {

    def setup() {
        Global.session = Mock(Session) { getConfig()>>[:] }
    }

    def 'should create a publish dir obj'() {
        given:
        def publish

        when:
        publish = PublishDir.create(path: '/data')
        then:
        publish.path == Paths.get('/data')

        when:
        publish = PublishDir.create(path: 'data')
        then:
        publish.path == Paths.get('data').complete()

        when:
        publish = PublishDir.create( path: Paths.get('data') )
        then:
        publish.path == Paths.get('data').complete()

        when:
        publish = PublishDir.create( [path: '/some/dir', overwrite: true, pattern: '*.bam', mode: 'link'] )
        then:
        publish.path == Paths.get('/some/dir')
        publish.mode == FilePublisher.Mode.LINK
        publish.pattern == '*.bam'
        publish.overwrite
        publish.enabled

        when:
        publish = PublishDir.create( [path: '/some/data', mode: 'copy', enabled: false] )
        then:
        publish.path == Paths.get('/some/data')
        publish.mode == FilePublisher.Mode.COPY
        publish.pattern == null
        publish.overwrite == null
        !publish.enabled

        when:
        publish = PublishDir.create( [path: '/some/data', mode: 'copy', enabled: 'false'] )
        then:
        publish.path == Paths.get('/some/data')
        publish.mode == FilePublisher.Mode.COPY
        publish.pattern == null
        publish.overwrite == null
        !publish.enabled

        when:
        publish = PublishDir.create( [path:'this/folder', overwrite: false, pattern: '*.txt', mode: 'copy'] )
        then:
        publish.path == Paths.get('this/folder').complete()
        publish.mode == FilePublisher.Mode.COPY
        publish.pattern == '*.txt'
        publish.overwrite == false

    }

    def 'should create publish dir with extended params' () {
        given:
        def publish

        when:
        publish = PublishDir.create(tags: ['foo','bar'])
        then:
        publish.@tags == ['foo','bar']

        when:
        publish = PublishDir.create(contentType: 'text/json')
        then:
        publish.@contentType == 'text/json'

        when:
        publish = PublishDir.create(storageClass: 'xyz')
        then:
        publish.@storageClass == 'xyz'
    }

    def 'should create symlinks for output files' () {
        given:
        def folder = Files.createTempDirectory('nxf')
        folder.resolve('work-dir').mkdir()
        folder.resolve('work-dir/file1.txt').text = 'aaa'
        folder.resolve('work-dir/file2.bam').text = 'bbb'
        folder.resolve('work-dir/file3.fastq').text = 'ccc'
        folder.resolve('work-dir/file4.temp').text = 'zzz'

        def workDir = folder.resolve('work-dir')
        def publishDir = folder.resolve('pub-dir')
        def task = new TaskRun(workDir: workDir, config: new TaskConfig(), name: 'foo')

        when:
        def outputs = [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('file3.fastq')
        ] as Set
        def publisher = new PublishDir(path: publishDir)
        publisher.apply(outputs, task)

        then:
        publishDir.resolve('file1.txt').exists()
        publishDir.resolve('file2.bam').exists()
        publishDir.resolve('file3.fastq').exists()
        !publishDir.resolve('file3.temp').exists()

        publishDir.resolve('file1.txt').isLink()
        publishDir.resolve('file2.bam').isLink()
        publishDir.resolve('file3.fastq').isLink()


        when:
        outputs = [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('file3.fastq')
        ] as Set
        publishDir.deleteDir()
        publisher = new PublishDir(path: publishDir, pattern: '*.bam')
        publisher.apply( outputs, task )

        then:
        !publishDir.resolve('file1.txt').exists()
        publishDir.resolve('file2.bam').exists()
        !publishDir.resolve('file3.fastq').exists()
        !publishDir.resolve('file3.temp').exists()

        cleanup:
        folder?.deleteDir()

    }

    def 'should copy output files' () {

        given:
        def session = new Session()
        def folder = Files.createTempDirectory('nxf')
        folder.resolve('work-dir').mkdir()
        folder.resolve('work-dir/file1.txt').text = 'aaa'
        folder.resolve('work-dir/file2.bam').text = 'bbb'
        folder.resolve('work-dir/file3.fastq').text = 'ccc'
        folder.resolve('work-dir/dir-x').mkdir()
        folder.resolve('work-dir/dir-x').resolve('file.1').text = 'xxx'
        folder.resolve('work-dir/dir-x').resolve('file.2').text = 'yyy'
        folder.resolve('work-dir/dir-y').mkdir()
        folder.resolve('work-dir/dir-y').resolve('file.3').text = '333'
        folder.resolve('work-dir/dir-y').resolve('file.4').text = '444'

        def workDir = folder.resolve('work-dir')
        def publishDir = folder.resolve('pub-dir')
        def task = new TaskRun(workDir: workDir, config: new TaskConfig(), name: 'foo')

        when:
        def outputs = [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('dir-x'),
                workDir.resolve('dir-y/file.3')
        ] as Set
        def publisher = new PublishDir(path: publishDir, mode: 'copy')
        publisher.apply( outputs, task )
        and:
        session.@publishPoolManager.shutdown(false)

        then:
        publishDir.resolve('file1.txt').text == 'aaa'
        publishDir.resolve('file2.bam').text == 'bbb'
        publishDir.resolve('dir-x').isDirectory()
        publishDir.resolve('dir-x/file.1').text == 'xxx'
        publishDir.resolve('dir-x/file.2').text == 'yyy'
        publishDir.resolve('dir-y/file.3').text == '333'

        !publishDir.resolve('file3.fastq').exists()
        !publishDir.resolve('file3.temp').exists()
        !publishDir.resolve('dir-y/file.4').exists()

        !publishDir.resolve('file1.txt').isLink()
        !publishDir.resolve('file2.bam').isLink()
        !publishDir.resolve('file3.fastq').isLink()
        !publishDir.resolve('dir-x').isLink()

        cleanup:
        folder?.deleteDir()

    }

    def 'should apply saveAs closure' () {

        given:
        def folder = Files.createTempDirectory('nxf')
        folder.resolve('work-dir').mkdir()
        folder.resolve('pub-dir1').mkdir()
        folder.resolve('pub-dir2').mkdir()
        folder.resolve('work-dir/file1.txt').text = 'aaa'
        folder.resolve('work-dir/file2.bam').text = 'bbb'
        folder.resolve('work-dir/file3.fastq').text = 'ccc'
        folder.resolve('work-dir/file4.temp').text = 'zzz'

        def workDir = folder.resolve('work-dir')
        def target1 = folder.resolve('pub-dir1')
        def target2 = folder.resolve('pub-dir2')
        def task = new TaskRun(workDir: workDir, config: new TaskConfig(), name: 'foo')

        when:
        def outputs = [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('file3.fastq'),
                workDir.resolve('file4.temp')
        ] as Set
        def rule = { String it ->
            if( it == 'file1.txt' ) return 'file_one.txt'
            if( it !='file4.temp' ) return target2.resolve(it)
            return null
        }
        def publisher = new PublishDir(path: target1, saveAs: rule)
        publisher.apply( outputs, task )

        then:
        target1.resolve('file_one.txt').text == 'aaa'
        target2.resolve('file2.bam').text == 'bbb'
        target2.resolve('file3.fastq').text == 'ccc'
        !target1.resolve('file4.temp').exists()
        !target2.resolve('file4.temp').exists()
        workDir.resolve('file4.temp').exists()

        cleanup:
        folder.deleteDir()

    }

    def 'should check null path' () {
        given:
        def pub = new PublishDir()

        expect:
        !pub.checkNull('ciao')
        pub.checkNull('null')
        pub.checkNull('null-x')
        pub.checkNull(' null')
    }

    def 'should not apply disable rule' () {

        given:
        def folder = Files.createTempDirectory('nxf')
        folder.resolve('work-dir').mkdir()
        folder.resolve('work-dir/file1.txt').text = 'aaa'

        def workDir = folder.resolve('work-dir')
        def publishDir = folder.resolve('pub-dir')
        def task = new TaskRun(workDir: workDir, config: Mock(TaskConfig))

        when:
        def outputs = [
                workDir.resolve('file1.txt'),
        ] as Set
        def publisher = new PublishDir(path: publishDir, enabled: false)
        publisher.apply(outputs, task)

        then:
        !publishDir.resolve('file1.txt').exists()

        cleanup:
        folder?.deleteDir()

    }

    def 'should only return not overlapping paths' () {
        given:
        def publisher = new PublishDir()

        expect:
        publisher.dedupPaths( GIVEN.collect{Paths.get(it)} ).collect{it.toString()} == EXPECT

        where:
        GIVEN                                               | EXPECT
        ['/foo/bar.txt','/foo/bar.txt', '/foo']             | ['/foo']
        ['/foo/bar.txt','/foo/bar.txt', '/foo/delta']       | ['/foo/bar.txt', '/foo/delta']
        ['/foo/x1','/foo/x1/y', '/bar/x2']                  | ['/foo/x1','/bar/x2']
    }

    def 'should set failOnError via env variable' () {
        given:
        SysEnv.push(ENV)

        when:
        def publish = new PublishDir()
        then:
        publish.failOnError == EXPECTED
        cleanup:
        SysEnv.pop()

        where:
        ENV                                         | EXPECTED
        [:]                                         | true
        [NXF_PUBLISH_FAIL_ON_ERROR: 'true']         | true
        [NXF_PUBLISH_FAIL_ON_ERROR: 'false']        | false
    }


    def 'should use rellink as default mode when the stage in mode is rellink' () {
        given:
        def folder = Files.createTempDirectory('nxf')
        def workDir = folder.resolve('work-dir'); workDir.mkdir()
        workDir.resolve('file1.txt').text = 'aaa'
        def publishDir = folder.resolve('pub-dir')
        def task = new TaskRun(workDir: workDir, config: new TaskConfig(stageInMode: 'rellink'), name: 'foo')

        when:
        new PublishDir(path: publishDir).apply([workDir.resolve('file1.txt')] as Set, task)
        then:
        Files.isSymbolicLink(publishDir.resolve('file1.txt'))
        !Files.readSymbolicLink(publishDir.resolve('file1.txt')).isAbsolute()
        publishDir.resolve('file1.txt').text == 'aaa'

        cleanup:
        folder?.deleteDir()
    }

    def 'should reject an invalid publish mode' () {
        when:
        PublishDir.create(path: '/data', mode: 'nope')
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.startsWith "Invalid publish mode 'nope'"
    }

}
