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

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Global
import nextflow.Session
import nextflow.processor.PublishDir
import nextflow.util.CsvWriter
/**
 * Publish files from a source channel.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PublishOp {

    private DataflowReadChannel source

    private Map opts

    private Path targetDir

    private Closure pathAs

    private IndexOpts indexOpts

    private List indexRecords = []

    private volatile boolean complete

    private Session getSession() { Global.session as Session }

    PublishOp(DataflowReadChannel source, Map opts) {
        this.source = source
        this.opts = opts
        this.targetDir = opts.path as Path
        if( opts.pathAs instanceof Closure )
            this.pathAs = opts.pathAs as Closure
        if( opts.index )
            this.indexOpts = new IndexOpts(targetDir, opts.index as Map)
    }

    boolean getComplete() { complete }

    PublishOp apply() {
        final events = new HashMap(2)
        events.onNext = this.&onNext
        events.onComplete = this.&onComplete
        DataflowHelper.subscribeImpl(source, events)
        return this
    }

    protected void onNext(value) {
        log.trace "Publish operator received: $value"

        // evaluate dynamic path
        final path = pathAs != null
            ? pathAs.call(value)
            : targetDir
        if( path == null )
            return

        // emit workflow publish event
        session.notifyWorkflowPublish(value)

        // create publisher
        final overrides = path instanceof Closure
            ? [saveAs: path]
            : [path: path]
        final publisher = PublishDir.create(opts + overrides)

        // publish files
        final result = collectFiles([:], value)
        for( final entry : result ) {
            final sourceDir = entry.key
            final files = entry.value
            publisher.apply(files, sourceDir)
        }

        // create record for index file
        if( indexOpts ) {
            final record = indexOpts.mapper != null ? indexOpts.mapper.call(value) : value
            final normalized = normalizePaths(record)
            log.trace "Normalized record for index file: ${normalized}"
            indexRecords << normalized
        }
    }

    protected void onComplete(nope) {
        if( indexOpts && indexRecords.size() > 0 ) {
            log.trace "Saving records to index file: ${indexRecords}"
            final ext = indexOpts.path.getExtension()
            if( ext == 'csv' )
                new CsvWriter(header: indexOpts.header, sep: indexOpts.sep).apply(indexRecords, indexOpts.path)
            else if( ext == 'json' )
                indexOpts.path.text = DumpHelper.prettyPrint(indexRecords)
            else
                log.warn "Invalid extension '${ext}' for index file '${indexOpts.path}' -- should be 'csv' or 'json'"
            session.notifyFilePublish(indexOpts.path)
        }

        log.trace "Publish operator complete"
        this.complete = true
    }

    /**
     * Extract files from a received value for publishing.
     *
     * @param result
     * @param value
     */
    protected Map<Path,Set<Path>> collectFiles(Map<Path,Set<Path>> result, value) {
        if( value instanceof Path ) {
            final sourceDir = getTaskDir(value)
            if( sourceDir !in result )
                result[sourceDir] = new HashSet(10)
            result[sourceDir] << value
        }
        else if( value instanceof Collection ) {
            for( final el : value )
                collectFiles(result, el)
        }
        return result
    }

    /**
     * Normalize the paths in a record by converting
     * work directory paths to publish paths.
     *
     * @param value
     */
    protected Object normalizePaths(value) {
        if( value instanceof Path ) {
            return List.of(value.getBaseName(), normalizePath(value))
        }

        if( value instanceof Collection ) {
            return value.collect { el ->
                if( el instanceof Path )
                    return normalizePath(el)
                if( el instanceof Collection<Path> )
                    return normalizePaths(el)
                return el
            }
        }

        if( value instanceof Map ) {
            return value.collectEntries { k, v ->
                if( v instanceof Path )
                    return List.of(k, normalizePath(v))
                if( v instanceof Collection<Path> )
                    return List.of(k, normalizePaths(v))
                return List.of(k, v)
            }
        }

        throw new IllegalArgumentException("Index file record must be a list, map, or file: ${value} [${value.class.simpleName}]")
    }

    private Path normalizePath(Path path) {
        final sourceDir = getTaskDir(path)
        return targetDir.resolve(sourceDir.relativize(path)).normalize()
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
