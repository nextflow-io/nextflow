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

package nextflow.processor

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths
import java.util.concurrent.TimeUnit

import nextflow.Global
import nextflow.file.FileHelper
import spock.lang.Specification

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

        when:
        publish =  PublishDir.create( [path: '/some/data', mode: 'copy'] )
        then:
        publish.path == Paths.get('/some/data')
        publish.mode == PublishDir.Mode.COPY
        publish.pattern == null
        publish.overwrite == null


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

        final workDir = folder.resolve('work-dir')
        final publishDir = folder.resolve('pub-dir')
        final task = new TaskRun(workDir: workDir, config: Mock(TaskConfig))

        when:
        def outputs =  [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('file3.fastq')
        ]
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
        ]
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
        Global.session = null
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

        final workDir = folder.resolve('work-dir')
        final publishDir = folder.resolve('pub-dir')
        final task = new TaskRun(workDir: workDir, config: Mock(TaskConfig))

        when:
        def outputs = [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('dir-x'),
                workDir.resolve('dir-y/file.3')
        ]
        def publisher = new PublishDir(path: publishDir, mode: 'copy')
        publisher.apply( outputs, task )

        PublishDir.executor.shutdown()
        PublishDir.executor.awaitTermination(5, TimeUnit.SECONDS)

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

        final workDir = folder.resolve('work-dir')
        final target1 = folder.resolve('pub-dir1')
        final target2 = folder.resolve('pub-dir2')
        final task = new TaskRun(workDir: workDir, config: Mock(TaskConfig))

        when:
        def outputs = [
                workDir.resolve('file1.txt'),
                workDir.resolve('file2.bam'),
                workDir.resolve('file3.fastq'),
                workDir.resolve('file4.temp')
        ]
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
}
