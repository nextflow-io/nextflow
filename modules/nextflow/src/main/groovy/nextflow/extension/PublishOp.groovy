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
 */

package nextflow.extension

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.file.FileHelper
import nextflow.processor.PublishDir
import nextflow.util.CsvWriter
/**
 * Publish a workflow output.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PublishOp {

    private Session session

    private DataflowReadChannel source

    private Map opts

    private String path

    private Closure pathResolver

    private IndexOpts indexOpts

    private List indexRecords = []

    private volatile boolean complete

    PublishOp(Session session, DataflowReadChannel source, Map opts) {
        this.session = session
        this.source = source
        this.opts = opts
        this.path = opts.path as String
        if( opts.pathResolver instanceof Closure )
            this.pathResolver = opts.pathResolver as Closure
        if( opts.index )
            this.indexOpts = new IndexOpts(session.outputDir, opts.index as Map)
    }

    boolean getComplete() { complete }

    PublishOp apply() {
        final events = new HashMap(2)
        events.onNext = this.&onNext
        events.onComplete = this.&onComplete
        DataflowHelper.subscribeImpl(source, events)
        return this
    }

    /**
     * For each incoming value, perform the following:
     *
     * 1. Publish any files contained in the value
     * 2. Append a record to the index file for the value (if enabled)
     *
     * @param value
     */
    protected void onNext(value) {
        log.trace "Publish operator received: $value"

        // evaluate dynamic path
        final targetResolver = getTargetDir(value)
        if( targetResolver == null )
            return

        // emit workflow publish event
        session.notifyWorkflowPublish(value)

        // create publisher
        final overrides = targetResolver instanceof Closure
            ? [saveAs: targetResolver]
            : [path: targetResolver]
        final publisher = PublishDir.create(opts + overrides)

        // publish files
        final result = collectFiles([:], value)
        for( final entry : result ) {
            final sourceDir = entry.key
            final files = entry.value
            publisher.apply(files, sourceDir)
        }

        // append record to index file
        if( indexOpts ) {
            final record = indexOpts.mapper != null ? indexOpts.mapper.call(value) : value
            final normalized = normalizePaths(record, targetResolver)
            log.trace "Normalized record for index file: ${normalized}"
            indexRecords << normalized
        }
    }

    /**
     * Compute the target directory for a published value:
     *
     * @param value
     * @return Path | Closure<Path>
     */
    protected Object getTargetDir(value) {
        // if the publish path is a string, resolve it against
        // the base output directory
        final outputDir = session.outputDir
        if( pathResolver == null )
            return outputDir.resolve(path)

        // if the publish path is a closure, invoke it on the
        // published value
        final resolvedPath = pathResolver.call(value)

        // if the resolved path is null, don't publish it
        if( resolvedPath == null )
            return null

        // if the resolved publish path is a string, resolve it
        // against the base output directory
        if( resolvedPath instanceof CharSequence )
            return outputDir.resolve(resolvedPath.toString())

        // if the resolved publish path is a closure, use the closure
        // to transform each published file and resolve it against
        // the base output directory
        if( resolvedPath instanceof Closure )
            return { file -> outputDir.resolve(resolvedPath.call(file) as String) }

        throw new ScriptRuntimeException("Output `path` directive should return a string or closure, but instead returned a ${resolvedPath.class.name}")
    }

    /**
     * Once all values have been published, write the
     * index file (if enabled).
     */
    protected void onComplete(nope) {
        if( indexOpts && indexRecords.size() > 0 ) {
            log.trace "Saving records to index file: ${indexRecords}"
            final indexPath = indexOpts.path
            final ext = indexPath.getExtension()
            indexPath.parent.mkdirs()
            if( ext == 'csv' ) {
                new CsvWriter(header: indexOpts.header, sep: indexOpts.sep).apply(indexRecords, indexPath)
            }
            else if( ext == 'json' ) {
                indexPath.text = DumpHelper.prettyPrint(indexRecords)
            }
            else {
                log.warn "Invalid extension '${ext}' for index file '${indexPath}' -- should be 'csv' or 'json'"
            }
            session.notifyFilePublish(indexPath)
        }

        log.trace "Publish operator complete"
        this.complete = true
    }

    /**
     * Extract files from a received value for publishing.
     * Files external to the work directory are not published.
     *
     * @param result
     * @param value
     */
    protected Map<Path,Set<Path>> collectFiles(Map<Path,Set<Path>> result, value) {
        if( value instanceof Path ) {
            final sourceDir = getTaskDir(value)
            if( sourceDir != null ) {
                if( sourceDir !in result )
                    result[sourceDir] = new HashSet(10)
                result[sourceDir] << value
            }
        }
        else if( value instanceof Collection ) {
            for( final el : value )
                collectFiles(result, el)
        }
        else if( value instanceof Map ) {
            for( final entry : value.entrySet() )
                collectFiles(result, entry.value)
        }
        return result
    }

    /**
     * Transform a value (i.e. path, collection, or map) by
     * normalizing any paths within the value.
     *
     * @param value
     * @param targetResolver
     */
    protected Object normalizePaths(value, targetResolver) {
        if( value instanceof Path ) {
            return List.of(value.getBaseName(), normalizePath(value, targetResolver))
        }

        if( value instanceof Collection ) {
            return value.collect { el ->
                if( el instanceof Path )
                    return normalizePath(el, targetResolver)
                if( el instanceof Collection<Path> )
                    return normalizePaths(el, targetResolver)
                return el
            }
        }

        if( value instanceof Map ) {
            return value
                .findAll { k, v -> v != null }
                .collectEntries { k, v ->
                    if( v instanceof Path )
                        return Map.entry(k, normalizePath(v, targetResolver))
                    if( v instanceof Collection<Path> )
                        return Map.entry(k, normalizePaths(v, targetResolver))
                    return Map.entry(k, v)
                }
        }

        throw new IllegalArgumentException("Index file record must be a list, map, or file: ${value} [${value.class.simpleName}]")
    }

    /**
     * Convert a work directory path to the corresponding
     * publish destination.
     *
     * @param path
     * @param targetResolver
     */
    private Path normalizePath(Path path, targetResolver) {
        // if the source file does not reside in the work directory,
        // return it directly without any normalization
        final sourceDir = getTaskDir(path)
        if( sourceDir == null )
            return path

        // if the target resolver is a closure, use it to transform
        // the source filename to the target path
        if( targetResolver instanceof Closure<Path> )
            return (targetResolver.call(path.getName()) as Path).normalize()

        // if the target resolver is a directory, resolve the source
        // filename against it
        if( targetResolver instanceof Path )
            return targetResolver.resolve(sourceDir.relativize(path)).normalize()

        throw new IllegalStateException()
    }

    /**
     * Try to infer the parent task directory to which a path belongs. It
     * should be a directory starting with a session work dir and having
     * at lest two sub-directories, e.g. work/ab/cdef/etc
     *
     * @param path
     */
    protected Path getTaskDir(Path path) {
        if( path == null )
            return null
        return getTaskDir0(path, session.workDir.resolve('tmp'))
            ?: getTaskDir0(path, session.workDir)
            ?: getTaskDir0(path, session.bucketDir)
    }

    private Path getTaskDir0(Path file, Path base) {
        if( base == null )
            return null
        if( base.fileSystem != file.fileSystem )
            return null
        final len = base.nameCount
        if( file.startsWith(base) && file.getNameCount() > len+2 )
            return base.resolve(file.subpath(len,len+2))
        return null
    }

    static class IndexOpts {
        Path path
        Closure mapper
        def /* boolean | List<String> */ header = false
        String sep = ','

        IndexOpts(Path targetDir, Map opts) {
            this.path = targetDir.resolve(opts.path as String)

            if( opts.mapper )
                this.mapper = opts.mapper as Closure
            if( opts.header != null )
                this.header = opts.header
            if( opts.sep )
                this.sep = opts.sep as String
        }
    }

}
