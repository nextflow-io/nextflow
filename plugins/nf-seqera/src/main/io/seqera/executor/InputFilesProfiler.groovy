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
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.sched.api.schema.v1a1.InputFilesMetrics
import nextflow.file.FileHolder
import nextflow.processor.TaskRun

/**
 * Computes input files metrics for a task.
 * Follows symlinks and recursively computes directory sizes.
 *
 * Stats are cached per path, because the same reference files are
 * usually shared by many tasks of the same run. When an executor is
 * given, the files of a task are statted in parallel on it.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class InputFilesProfiler {

    /** Maximum number of paths whose stats are kept in the cache */
    static final long DEFAULT_CACHE_SIZE = 50_000

    /** Keyed by URI rather than Path, to avoid retaining heavyweight cloud path objects */
    private final Cache<URI,long[]> cache

    /** Paths whose stat failure has already been warned about, so a failing path is warned once */
    private final Set<URI> warned

    /** Executor to stat the files of a task in parallel, or null to stat them sequentially */
    private final ExecutorService executor

    InputFilesProfiler() {
        this(null)
    }

    InputFilesProfiler(ExecutorService executor) {
        this(executor, DEFAULT_CACHE_SIZE)
    }

    InputFilesProfiler(ExecutorService executor, long cacheSize) {
        this.executor = executor
        this.cache = CacheBuilder.newBuilder()
            .maximumSize(cacheSize)
            .build()
        this.warned = Collections.newSetFromMap(CacheBuilder.newBuilder()
            .maximumSize(cacheSize)
            .<URI,Boolean>build()
            .asMap())
    }

    /**
     * Compute input files metrics for a task.
     *
     * @param task The task to compute metrics for
     * @return InputFilesMetrics or null if no input files
     */
    InputFilesMetrics compute(TaskRun task) {
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
    InputFilesMetrics compute(List<FileHolder> files) {
        if( !files || files.isEmpty() )
            return null

        return compute0(files)
    }

    private InputFilesMetrics compute0(List<FileHolder> files) {
        int totalCount = 0
        long totalBytes = 0
        long maxFileBytes = Long.MIN_VALUE
        long minFileBytes = Long.MAX_VALUE

        for( long[] result : statAll(files) ) {
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

    private List<long[]> statAll(List<FileHolder> files) {
        final List<Callable<long[]>> tasks = files.collect { FileHolder fh -> { getFileStats(fh.storePath) } as Callable<long[]> }
        if( executor == null || tasks.size() == 1 )
            return tasks.collect { it.call() }
        // invokeAll waits interruptibly and cancels pending stats on interrupt
        return executor.invokeAll(tasks).collect { it.get() }
    }

    /**
     * Get file stats for a path: file count and total size.
     * For regular files, count is 1. For directories, count is the number of files within.
     * Results are cached; failures are not, so a transient error is retried by the next task.
     *
     * @param path The path to measure
     * @return A two-element array: [fileCount, totalSize]
     */
    private long[] getFileStats(Path path) {
        if( path == null )
            return new long[]{0, 0}

        final uri = path.toUri()
        try {
            return cache.get(uri, { statPath(path) } as Callable<long[]>)
        }
        catch( Exception e ) {
            final cause = e.cause ?: e
            // @Slf4j guards the debug message, so it is only built when debug logging is enabled
            if( warned.add(uri) )
                log.warn "Unable to determine size for input file: ${path} - ${cause.message}"
            else
                log.debug "Unable to determine size for input file: ${path} - ${cause.message}"
            return new long[]{1, 0}
        }
    }

    /**
     * Read file stats for a path from its file system, bypassing the cache.
     */
    protected long[] statPath(Path path) {
        if( Files.isDirectory(path) ) {
            return computeDirStats(path)
        }
        // Files.size() follows symlinks by default
        return new long[]{1, Files.size(path)}
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
