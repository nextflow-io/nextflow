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
import nextflow.processor.PublishDir
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.WorkflowOutputEvent
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

    private String name

    private DataflowReadChannel source

    private Map publishOpts

    private String path

    private Closure pathResolver

    private IndexOpts indexOpts

    private List indexRecords = []

    private volatile boolean complete

    PublishOp(Session session, String name, DataflowReadChannel source, Map opts) {
        this.session = session
        this.name = name
        this.source = source
        this.publishOpts = opts
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

        // create publisher
        final overrides = new LinkedHashMap()
        if( targetResolver instanceof Closure )
            overrides.saveAs = targetResolver
        else
            overrides.path = targetResolver

        final publisher = PublishDir.create(publishOpts + overrides)

        // publish files
        final result = collectFiles([:], value)
        for( final entry : result ) {
            final sourceDir = entry.key
            final files = entry.value
            publisher.apply(files, sourceDir)
        }

        // append record to index
        final normalized = normalizePaths(value, targetResolver)
        log.trace "Normalized record for index file: ${normalized}"
        indexRecords << normalized
    }

    /**
     * Compute the target directory for a published value.
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
        final dsl = new PublishDsl()
        final cl = (Closure)pathResolver.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        final resolvedPath = cl.call(value)

        // if the closure contained publish statements, use
        // the resulting mapping to create a saveAs closure
        final mapping = dsl.build()
        if( mapping instanceof Map<String,String> )
            return { filename -> outputDir.resolve(mapping[filename]) }

        // if the resolved publish path is a string, resolve it
        // against the base output directory
        if( resolvedPath instanceof CharSequence )
            return outputDir.resolve(resolvedPath.toString())

        throw new ScriptRuntimeException("Invalid output `path` directive -- it should either return a string or use the `>>` operator to publish files")
    }

    private class PublishDsl {
        private Map<String,String> mapping = null

        void publish(Object source, String target) {
            if( source == null )
                return
            if( source instanceof Path ) {
                publish0(source, target)
            }
            else if( source instanceof Collection<Path> ) {
                if( !target.endsWith('/') )
                    throw new ScriptRuntimeException("Invalid publish target '${target}' -- should be a directory (end with a `/`) when publishing a collection of files")
                for( final path : source )
                    publish0(path, target)
            }
            else {
                throw new ScriptRuntimeException("Publish source should be a file or collection of files, but received a ${source.class.name}")
            }
        }

        private void publish0(Path source, String target) {
            log.trace "Publishing ${source} to ${target}"
            if( mapping == null )
                mapping = [:]
            final filename = getTaskDir(source).relativize(source).toString()
            final resolved = target.endsWith('/')
                ? target + filename
                : target
            mapping[filename] = resolved
        }

        Map<String,String> build() {
            return mapping
        }
    }

    /**
     * Once all channel values have been published, publish the final
     * workflow output and index file (if enabled).
     */
    protected void onComplete(nope) {
        // publish individual record if source is a value channel
        final value = CH.isValue(source)
            ? indexRecords.first()
            : indexRecords

        // publish workflow output
        final indexPath = indexOpts ? indexOpts.path : null
        session.notifyWorkflowOutput(new WorkflowOutputEvent(name, value, indexPath))

        // write value to index file
        if( indexOpts ) {
            final ext = indexPath.getExtension()
            indexPath.parent.mkdirs()
            if( ext == 'csv' ) {
                new CsvWriter(header: indexOpts.header, sep: indexOpts.sep).apply(indexRecords, indexPath)
            }
            else if( ext == 'json' ) {
                indexPath.text = DumpHelper.prettyPrintJson(value)
            }
            else if( ext == 'yaml' || ext == 'yml' ) {
                indexPath.text = DumpHelper.prettyPrintYaml(value)
            }
            else {
                log.warn "Invalid extension '${ext}' for index file '${indexPath}' -- should be CSV, JSON, or YAML"
            }
            session.notifyFilePublish(new FilePublishEvent(null, indexPath, publishOpts.labels as List))
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
            return normalizePath(value, targetResolver)
        }

        if( value instanceof Collection ) {
            return value.collect { el ->
                if( el instanceof Path )
                    return normalizePath(el, targetResolver)
                if( el instanceof Collection<Path> )
                    return normalizePaths(el, targetResolver)
                if( el instanceof Map )
                    return normalizePaths(el, targetResolver)
                return el
            }
        }

        if( value instanceof Map ) {
            return value.collectEntries { k, v ->
                if( v instanceof Path )
                    return Map.entry(k, normalizePath(v, targetResolver))
                if( v instanceof Collection<Path> )
                    return Map.entry(k, normalizePaths(v, targetResolver))
                if( v instanceof Map )
                    return Map.entry(k, normalizePaths(v, targetResolver))
                return [k, v]
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
        if( targetResolver instanceof Path ) {
            // note: make sure to convert the relative path to as a string to prevent
            // an exception when mixing different path providers e.g. local fs and remove cloud
            // thrown by {@link Path#resolve) method
            final relPath = sourceDir.relativize(path).toString()
            return targetResolver.resolve(relPath).normalize()
        }

        throw new IllegalStateException("Unexpected targetResolver argument: ${targetResolver}")
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
        def /* boolean | List<String> */ header = false
        String sep = ','

        IndexOpts(Path targetDir, Map opts) {
            this.path = targetDir.resolve(opts.path as String)

            if( opts.header != null )
                this.header = opts.header
            if( opts.sep )
                this.sep = opts.sep as String
        }
    }

}
