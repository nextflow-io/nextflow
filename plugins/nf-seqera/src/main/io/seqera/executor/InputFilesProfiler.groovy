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

package io.seqera.executor

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.sched.api.schema.v1a1.InputFilesMetrics
import nextflow.file.FileHolder
import nextflow.processor.TaskRun

/**
 * Computes input files metrics for a task.
 * Follows symlinks and recursively computes directory sizes.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class InputFilesProfiler {

    /**
     * Compute input files metrics for a task.
     *
     * @param task The task to compute metrics for
     * @return InputFilesMetrics or null if no input files
     */
    static InputFilesMetrics compute(TaskRun task) {
        final files = task?.inputFiles
        if( !files || files.isEmpty() )
            return null

        return compute0(files)
    }

    /**
     * Compute input files metrics from a list of file holders.
     *
     * @param files List of FileHolder objects
     * @return InputFilesMetrics or null if list is empty
     */
    static InputFilesMetrics compute(List<FileHolder> files) {
        if( !files || files.isEmpty() )
            return null

        return compute0(files)
    }

    private static InputFilesMetrics compute0(List<FileHolder> files) {
        int totalCount = 0
        long totalBytes = 0
        long maxFileBytes = Long.MIN_VALUE
        long minFileBytes = Long.MAX_VALUE

        for( FileHolder fh : files ) {
            final long[] result = getFileStats(fh.storePath)
            final count = (int) result[0]
            final size = result[1]
            totalCount += count
            totalBytes += size
            if( size > maxFileBytes )
                maxFileBytes = size
            if( size < minFileBytes )
                minFileBytes = size
        }

        return new InputFilesMetrics()
            .count(totalCount)
            .totalBytes(totalBytes)
            .maxFileBytes(maxFileBytes)
            .minFileBytes(minFileBytes)
    }

    /**
     * Get file stats for a path: file count and total size.
     * For regular files, count is 1. For directories, count is the number of files within.
     *
     * @param path The path to measure
     * @return A two-element array: [fileCount, totalSize]
     */
    private static long[] getFileStats(Path path) {
        if( path == null )
            return new long[]{0, 0}

        try {
            if( Files.isDirectory(path) ) {
                return computeDirStats(path)
            }
            // Files.size() follows symlinks by default
            return new long[]{1, Files.size(path)}
        }
        catch( Exception e ) {
            log.warn "Unable to determine size for input file: ${path} - ${e.message}"
            return new long[]{1, 0}
        }
    }

    /**
     * Recursively compute file count and total size of a directory.
     *
     * @param dir The directory path
     * @return A two-element array: [fileCount, totalSize]
     */
    private static long[] computeDirStats(Path dir) {
        final long[] result = [0L, 0L]  // [count, size]

        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                result[0]++
                result[1] += attrs.size()
                return FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFileFailed(Path file, IOException exc) {
                log.warn "Unable to access file during size computation: ${file} - ${exc.message}"
                return FileVisitResult.CONTINUE
            }
        })

        return result
    }
}
