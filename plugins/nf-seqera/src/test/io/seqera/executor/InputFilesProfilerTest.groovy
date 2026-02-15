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
 * Tests for InputFilesProfiler
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class InputFilesProfilerTest extends Specification {

    @TempDir
    Path tempDir

    def 'should return null for null task'() {
        expect:
        InputFilesProfiler.compute((TaskRun) null) == null
    }

    def 'should return null for task with no input files'() {
        given:
        def task = Mock(TaskRun) {
            getInputFiles() >> []
        }

        expect:
        InputFilesProfiler.compute(task) == null
    }

    def 'should return null for empty file list'() {
        expect:
        InputFilesProfiler.compute([]) == null
    }

    def 'should return null for null file list'() {
        expect:
        InputFilesProfiler.compute((List) null) == null
    }

    def 'should compute metrics for single small file'() {
        given:
        def file = tempDir.resolve('small.txt')
        Files.write(file, 'hello'.bytes)
        def files = [new FileHolder(file)]

        when:
        def metrics = InputFilesProfiler.compute(files)

        then:
        metrics.count == 1
        metrics.totalBytes == 5
        metrics.maxFileBytes == 5
        metrics.minFileBytes == 5
    }

    def 'should compute metrics for multiple files'() {
        given:
        def smallFile = tempDir.resolve('small.txt')
        Files.write(smallFile, new byte[500])

        def mediumFile = tempDir.resolve('medium.dat')
        Files.write(mediumFile, new byte[2000])

        def largeFile = tempDir.resolve('large.dat')
        Files.write(largeFile, new byte[50000])

        def files = [
            new FileHolder(smallFile),
            new FileHolder(mediumFile),
            new FileHolder(largeFile)
        ]

        when:
        def metrics = InputFilesProfiler.compute(files)

        then:
        metrics.count == 3
        metrics.totalBytes == 500 + 2000 + 50000
        metrics.maxFileBytes == 50000
        metrics.minFileBytes == 500
    }

    def 'should count files in directory recursively'() {
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
        def metrics = InputFilesProfiler.compute(files)

        then:
        metrics.count == 3  // 3 actual files inside the directory
        metrics.totalBytes == 600  // 100 + 200 + 300
        metrics.maxFileBytes == 600
        metrics.minFileBytes == 600
    }

    def 'should count files and directory contents together'() {
        given:
        def file1 = tempDir.resolve('input.fq')
        Files.write(file1, new byte[5000])

        def dir = tempDir.resolve('index')
        Files.createDirectory(dir)
        Files.write(dir.resolve('a.bin'), new byte[100])
        Files.write(dir.resolve('b.bin'), new byte[200])

        def files = [new FileHolder(file1), new FileHolder(dir)]

        when:
        def metrics = InputFilesProfiler.compute(files)

        then:
        metrics.count == 3  // 1 regular file + 2 files in directory
        metrics.totalBytes == 5300
        metrics.maxFileBytes == 5000
        metrics.minFileBytes == 300
    }

    def 'should follow symlinks'() {
        given:
        def realFile = tempDir.resolve('real.txt')
        Files.write(realFile, new byte[1000])

        def symlink = tempDir.resolve('link.txt')
        Files.createSymbolicLink(symlink, realFile)

        def files = [new FileHolder(symlink)]

        when:
        def metrics = InputFilesProfiler.compute(files)

        then:
        metrics.count == 1
        metrics.totalBytes == 1000
        metrics.maxFileBytes == 1000
        metrics.minFileBytes == 1000
    }

    def 'should handle non-existent file gracefully'() {
        given:
        def missingFile = tempDir.resolve('does-not-exist.txt')
        def files = [new FileHolder(missingFile)]

        when:
        def metrics = InputFilesProfiler.compute(files)

        then:
        metrics.count == 1
        metrics.totalBytes == 0
    }

    def 'should compute from TaskRun'() {
        given:
        def file = tempDir.resolve('task-input.txt')
        Files.write(file, new byte[2048])

        def task = Mock(TaskRun) {
            getInputFiles() >> [new FileHolder(file)]
        }

        when:
        def metrics = InputFilesProfiler.compute(task)

        then:
        metrics.count == 1
        metrics.totalBytes == 2048
        metrics.maxFileBytes == 2048
        metrics.minFileBytes == 2048
    }
}
