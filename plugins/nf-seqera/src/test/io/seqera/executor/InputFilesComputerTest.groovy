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
 *
 */

package io.seqera.executor

import java.nio.file.Files
import java.nio.file.Path

import nextflow.file.FileHolder
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Tests for InputFilesComputer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class InputFilesComputerTest extends Specification {

    @TempDir
    Path tempDir

    def 'should return null for null task'() {
        expect:
        InputFilesComputer.compute((TaskRun) null) == null
    }

    def 'should return null for task with no input files'() {
        given:
        def task = Mock(TaskRun) {
            getInputFiles() >> []
        }

        expect:
        InputFilesComputer.compute(task) == null
    }

    def 'should return null for empty file list'() {
        expect:
        InputFilesComputer.compute([]) == null
    }

    def 'should return null for null file list'() {
        expect:
        InputFilesComputer.compute((List) null) == null
    }

    def 'should compute metrics for single small file'() {
        given:
        def file = tempDir.resolve('small.txt')
        Files.write(file, 'hello'.bytes)
        def files = [new FileHolder(file)]

        when:
        def metrics = InputFilesComputer.compute(files)

        then:
        metrics.count == 1
        metrics.totalBytes == 5
        metrics.bins.size() == 5
        metrics.bins[0].range == '<=1MB'
        metrics.bins[0].count == 1
        metrics.bins[1].count == 0
        metrics.bins[2].count == 0
        metrics.bins[3].count == 0
        metrics.bins[4].count == 0
    }

    def 'should compute metrics for multiple files in different bins'() {
        given:
        // Create files of different sizes
        def smallFile = tempDir.resolve('small.txt')
        Files.write(smallFile, new byte[500])  // 500 bytes - <=1MB bin

        def mediumFile = tempDir.resolve('medium.dat')
        Files.write(mediumFile, new byte[2 * 1024 * 1024])  // 2MB - <=10MB bin

        def largeFile = tempDir.resolve('large.dat')
        Files.write(largeFile, new byte[50 * 1024 * 1024])  // 50MB - <=100MB bin

        def files = [
            new FileHolder(smallFile),
            new FileHolder(mediumFile),
            new FileHolder(largeFile)
        ]

        when:
        def metrics = InputFilesComputer.compute(files)

        then:
        metrics.count == 3
        metrics.totalBytes == 500 + (2 * 1024 * 1024) + (50 * 1024 * 1024)
        metrics.bins[0].count == 1  // <=1MB
        metrics.bins[1].count == 1  // <=10MB
        metrics.bins[2].count == 1  // <=100MB
        metrics.bins[3].count == 0  // <=1GB
        metrics.bins[4].count == 0  // >1GB
    }

    def 'should compute directory size recursively'() {
        given:
        def dir = tempDir.resolve('mydir')
        Files.createDirectory(dir)
        Files.write(dir.resolve('file1.txt'), new byte[100])
        Files.write(dir.resolve('file2.txt'), new byte[200])

        def subDir = dir.resolve('subdir')
        Files.createDirectory(subDir)
        Files.write(subDir.resolve('file3.txt'), new byte[300])

        def files = [new FileHolder(dir)]

        when:
        def metrics = InputFilesComputer.compute(files)

        then:
        metrics.count == 1
        metrics.totalBytes == 600  // 100 + 200 + 300
        metrics.bins[0].count == 1  // <=1MB
    }

    def 'should follow symlinks'() {
        given:
        def realFile = tempDir.resolve('real.txt')
        Files.write(realFile, new byte[1000])

        def symlink = tempDir.resolve('link.txt')
        Files.createSymbolicLink(symlink, realFile)

        def files = [new FileHolder(symlink)]

        when:
        def metrics = InputFilesComputer.compute(files)

        then:
        metrics.count == 1
        metrics.totalBytes == 1000
    }

    def 'should handle non-existent file gracefully'() {
        given:
        def missingFile = tempDir.resolve('does-not-exist.txt')
        def files = [new FileHolder(missingFile)]

        when:
        def metrics = InputFilesComputer.compute(files)

        then:
        metrics.count == 1
        metrics.totalBytes == 0  // Unable to read, returns 0
    }

    def 'should compute from TaskRun'() {
        given:
        def file = tempDir.resolve('task-input.txt')
        Files.write(file, new byte[2048])

        def task = Mock(TaskRun) {
            getInputFiles() >> [new FileHolder(file)]
        }

        when:
        def metrics = InputFilesComputer.compute(task)

        then:
        metrics.count == 1
        metrics.totalBytes == 2048
    }

    def 'should have correct bin labels'() {
        given:
        def file = tempDir.resolve('tiny.txt')
        Files.write(file, 'x'.bytes)
        def files = [new FileHolder(file)]

        when:
        def metrics = InputFilesComputer.compute(files)

        then:
        metrics.bins*.range == ['<=1MB', '<=10MB', '<=100MB', '<=1GB', '>1GB']
    }
}
