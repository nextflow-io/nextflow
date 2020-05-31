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

package nextflow.processor

import spock.lang.Specification

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths
import java.util.concurrent.TimeUnit

import nextflow.Session
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PublishDirTest extends Specification {

    def 'should create a publish dir obj'() {

        PublishDir publish

        when:
        publish = PublishDir.create(path: '/data')
        then:
        publish.path == Paths.get('/data')

        when:
        publish =  PublishDir.create(path: 'data')
        then:
        publish.path == Paths.get('data').complete()

        when:
        publish =  PublishDir.create( path: Paths.get('data') )
        then:
        publish.path == Paths.get('data').complete()

        when:
        publish =  PublishDir.create( [path: '/some/dir', overwrite: true, pattern: '*.bam', mode: 'link'] )
        then:
        publish.path == Paths.get('/some/dir')
        publish.mode == PublishDir.Mode.LINK
        publish.pattern == '*.bam'
        publish.overwrite
        publish.enabled

        when:
        publish =  PublishDir.create( [path: '/some/data', mode: 'copy', enabled: false] )
        then:
        publish.path == Paths.get('/some/data')
        publish.mode == PublishDir.Mode.COPY
        publish.pattern == null
        publish.overwrite == null
        !publish.enabled

        when:
        publish =  PublishDir.create( [path:'this/folder', overwrite: false, pattern: '*.txt', mode: 'copy'] )
        then:
        publish.path == Paths.get('this/folder').complete()
        publish.mode == PublishDir.Mode.COPY
        publish.pattern == '*.txt'
        publish.overwrite == false

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
        def task = new TaskRun(workDir: workDir, config: Mock(TaskConfig))

        when:
        def outputs =  [
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
        def task = new TaskRun(workDir: workDir, config: Mock(TaskConfig))

        when:
        def outputs = [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('dir-x'),
                workDir.resolve('dir-y/file.3')
        ] as Set
        def publisher = new PublishDir(path: publishDir, mode: 'copy')
        publisher.apply( outputs, task )

        session.fileTransferThreadPool.shutdown()
        session.fileTransferThreadPool.awaitTermination(5, TimeUnit.SECONDS)

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
        def task = new TaskRun(workDir: workDir, config: Mock(TaskConfig))

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
        def publisher = new PublishDir(path: target1, saveAs: rule )
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

    def 'should default mode to `symlink`' () {

        given:
        def processor = [:] as TaskProcessor
        processor.name = 'foo'

        def targetDir = Paths.get('/scratch/dir')
        def publisher = new PublishDir(path: targetDir, sourceFileSystem: FileSystems.default)

        when:
        publisher.validatePublishMode()
        then:
        publisher.mode == PublishDir.Mode.SYMLINK
    }


    def 'should change mode to `copy` when the target is a foreign file system' () {

        given:
        def workDirFileSystem = TestHelper.createInMemTempDir().fileSystem
        def processor = [:] as TaskProcessor
        processor.name = 'foo'

        def targetDir = TestHelper.createInMemTempDir()
        def publisher = new PublishDir(mode:'symlink', path: targetDir, sourceFileSystem: workDirFileSystem)

        when:
        publisher.validatePublishMode()
        then:
        publisher.mode == PublishDir.Mode.COPY
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
        def outputs =  [
                workDir.resolve('file1.txt'),
        ] as Set
        def publisher = new PublishDir(path: publishDir, enabled: false)
        publisher.apply(outputs, task)

        then:
        !publishDir.resolve('file1.txt').exists()

        cleanup:
        folder?.deleteDir()

    }
}
