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

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.sched.api.schema.v1a1.InputFilesBin
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
class InputFilesComputer {

    static final long KB = 1024L
    static final long MB = KB * 1024
    static final long GB = MB * 1024

    private static final List<Long> THRESHOLDS = [1*MB, 10*MB, 100*MB, 1*GB]
    private static final List<String> LABELS = ['<=1MB', '<=10MB', '<=100MB', '<=1GB', '>1GB']

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
        final sizes = new ArrayList<Long>(files.size())

        for( FileHolder fh : files ) {
            sizes.add(getSize(fh.storePath))
        }

        final binCounts = new int[5]
        long totalBytes = 0

        for( long size : sizes ) {
            totalBytes += size
            binCounts[findBinIndex(size)]++
        }

        final bins = new ArrayList<InputFilesBin>(5)
        for( int i = 0; i < LABELS.size(); i++ ) {
            bins.add(new InputFilesBin().range(LABELS[i]).count(binCounts[i]))
        }

        return new InputFilesMetrics()
            .count(sizes.size())
            .totalBytes(totalBytes)
            .bins(bins)
    }

    /**
     * Find the bin index for a given file size.
     */
    private static int findBinIndex(long size) {
        for( int i = 0; i < THRESHOLDS.size(); i++ ) {
            if( size <= THRESHOLDS[i] )
                return i
        }
        return 4  // >1GB
    }

    /**
     * Get the size of a path, following symlinks and recursively summing directories.
     *
     * @param path The path to measure
     * @return Size in bytes, or 0 if unable to determine
     */
    private static long getSize(Path path) {
        if( path == null )
            return 0

        try {
            if( Files.isDirectory(path) ) {
                return computeDirSize(path)
            }
            // Files.size() follows symlinks by default
            return Files.size(path)
        }
        catch( Exception e ) {
            log.warn "Unable to determine size for input file: ${path} - ${e.message}"
            return 0
        }
    }

    /**
     * Recursively compute the total size of a directory.
     *
     * @param dir The directory path
     * @return Total size in bytes
     */
    private static long computeDirSize(Path dir) {
        final long[] total = [0L]

        Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
                total[0] += attrs.size()
                return FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFileFailed(Path file, IOException exc) {
                log.warn "Unable to access file during size computation: ${file} - ${exc.message}"
                return FileVisitResult.CONTINUE
            }
        })

        return total[0]
    }
}
