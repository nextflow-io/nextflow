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

package nextflow.extension

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.file.FilePublisher
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

    private boolean includeInputs

    private List publishedValues = []

    private DataflowVariable target

    @Lazy
    private FilePublisher publisher = { new FilePublisher(session, publishOpts) }()

    PublishOp(Session session, String name, DataflowReadChannel source, Map opts) {
        this.session = session
        this.name = name
        this.source = source
        this.publishOpts = opts
        this.path = opts.path as String
        if( opts.pathResolver instanceof Closure )
            this.pathResolver = opts.pathResolver as Closure
        if( opts.index )
            this.indexOpts = new IndexOpts(opts.index as Map)
        this.includeInputs = opts.includeInputs as boolean
    }

    DataflowVariable apply() {
        this.target = new DataflowVariable()
        final events = new HashMap(2)
        events.onNext = { value ->
            safeExecute { onNext(value) }
        }
        events.onComplete = {
            safeExecute { onComplete() }
        }
        DataflowHelper.subscribeImpl(source, events)
        return target
    }

    /**
     * Perform an action. If an exception is raised, bind the exception
     * to the target, abort the run, and don't perform any more actions.
     *
     * @param action
     */
    private void safeExecute(Runnable action) {
        if( target.isError() )
            return
        try {
            action.run()
        } catch( Throwable e ) {
            // bind the error before aborting, since the abort
            // can interrupt the current thread
            target.bindError(e)
            log.error("@unknown", e)
            session.abort(e)
        }
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
        log.trace "Received value for workflow output '${name}': ${value}"

        // if output directory is disabled, report output files by
        // their work directory path instead of publishing them
        if( session.outputDir == null ) {
            publishedValues << value
            return
        }

        // resolve the target path of every file in the value
        final targets = resolveTargets(value)

        // publish the files
        checkTargetConflicts(targets)
        for( final entry : targets )
            publisher.publish(entry.key, entry.value)

        // publish value to workflow output
        final normalizedValue = normalizeValue(value, targets)

        log.trace "Published value to workflow output '${name}': ${normalizedValue}"
        publishedValues << normalizedValue
    }

    /**
     * Report two different files in the same value being published to the same
     * target, which would otherwise silently publish whichever file is written first.
     *
     * Files from different values can be published to the same target, in which
     * case the last file wins (e.g. a `versions.yml` file emitted by every task).
     *
     * @param targets
     */
    protected void checkTargetConflicts(Map<Path,Path> targets) {
        final sources = new HashMap<Path,Path>(targets.size())
        for( final entry : targets ) {
            final source = entry.key
            final targetPath = entry.value
            final previous = sources.putIfAbsent(targetPath, source)
            if( previous != null )
                throw new ScriptRuntimeException("Publish target '${targetPath.toUriString()}' for workflow output '${name}' is used by more than one file -- offending files: ${previous.toUriString()}, ${source.toUriString()}")
        }
    }

    /**
     * Resolve the target path of every file in a published value.
     *
     * @param value
     * @return Mapping of source file to absolute target path
     */
    protected Map<Path,Path> resolveTargets(value) {
        // if the publish path is a string, resolve every file in the value
        // against it
        if( pathResolver == null )
            return collectTargets(value, getTargetDir(path))

        // if the publish path is a closure, invoke it on the published value
        final dsl = new PublishDsl()
        final cl = (Closure)pathResolver.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        final resolvedPath = cl.call(value)

        // if the resolved publish path is a string, resolve it
        // against the base output directory
        if( resolvedPath instanceof CharSequence )
            return collectTargets(value, getTargetDir(resolvedPath.toString()))

        // if the closure returned a map of source -> target pairs,
        // treat it the same as a set of publish statements
        if( resolvedPath instanceof Map )
            for( final entry : resolvedPath.entrySet() )
                dsl.publish(entry.key, entry.value as String)

        // if the closure contained publish statements, resolve each declared
        // target against the base output directory
        final mapping = dsl.build()
        if( mapping != null ) {
            final result = new LinkedHashMap<Path,Path>(mapping.size())
            for( final entry : mapping )
                result.put(entry.key, getTargetDir(entry.value))
            return result
        }

        throw new ScriptRuntimeException("Invalid `path` directive for workflow output '${name}' -- expected a string, a map, or publish statements, but received: ${resolvedPath} [${resolvedPath?.class?.simpleName}]")
    }

    /**
     * Resolve a publish path against the base output directory.
     *
     * @param path
     */
    protected Path getTargetDir(String path) {
        return session.outputDir.resolve(path).normalize()
    }

    /**
     * Map every file in a value to its target path in a given directory.
     *
     * @param value
     * @param targetDir
     */
    protected Map<Path,Path> collectTargets(value, Path targetDir) {
        final result = new LinkedHashMap<Path,Path>()
        for( final file : collectFiles(new LinkedHashSet<Path>(), value) )
            result.put(file, targetDir.resolve(relativePath(file)).normalize())
        return result
    }

    /**
     * Determine whether a file should be published. Files that do not
     * originate from the work directory are published only when
     * `includeInputs` is enabled.
     *
     * @param file
     */
    protected boolean shouldPublish(Path file) {
        return includeInputs || getTaskDir(file) != null
    }

    /**
     * Get the path of a file relative to its task directory, or the
     * file name if the file is not in the work directory.
     *
     * The path is returned as a string to prevent a ProviderMismatchException
     * when the source and target use different path providers.
     *
     * @param file
     */
    protected String relativePath(Path file) {
        final sourceDir = getTaskDir(file)
        return sourceDir != null
            ? sourceDir.relativize(file).toString()
            : file.getFileName().toString()
    }

    private class PublishDsl {
        /**
         * Mapping of source file to target path, relative to the output directory.
         * It is keyed by the full source path because the same relative filename
         * can be published from multiple task directories.
         */
        private Map<Path,String> mapping = null

        void publish(Object source, String target) {
            // a no-op publish statement should still publish nothing
            // instead of falling back to the closure return value
            if( mapping == null )
                mapping = [:]
            if( source == null || target == null )
                return
            if( source instanceof Path ) {
                publish0(source, target)
            }
            else if( source instanceof Collection<Path> ) {
                if( !target.endsWith('/') )
                    throw new ScriptRuntimeException("Invalid publish target '${target}' for workflow output '${name}' -- should be a directory (end with a `/`) when publishing a collection of files")
                for( final path : source )
                    publish0(path, target)
            }
            else {
                throw new ScriptRuntimeException("Invalid publish source for workflow output '${name}' -- expected a file or collection of files, but received: ${source} [${source.class.simpleName}]")
            }
        }

        private void publish0(Path source, String target) {
            if( !shouldPublish(source) )
                return
            log.trace "Publishing ${source} to ${target}"
            final resolved = target.endsWith('/')
                ? target + relativePath(source)
                : target
            mapping[source] = resolved
        }

        Map<Path,String> build() {
            return mapping
        }
    }

    /**
     * Once all channel values have been published, publish the final
     * workflow output and index file (if enabled).
     */
    protected void onComplete() {
        // publish individual record if source is a value channel
        // NOTE: handle (invalid) empty dataflow value from legacy process, legacy collect(), etc
        final outputValue = CH.isValue(source)
            ? (publishedValues ? publishedValues.first() : null)
            : publishedValues

        // publish workflow output
        final indexPath = session.outputDir && indexOpts
            ? session.outputDir.resolve(indexOpts.path)
            : null
        session.notifyWorkflowOutput(new WorkflowOutputEvent(name, outputValue, indexPath))

        // write value to index file
        if( indexPath ) {
            final ext = indexPath.getExtension()
            indexPath.parent.mkdirs()
            if( ext == 'csv' ) {
                new CsvWriter(header: indexOpts.header, sep: indexOpts.sep).apply(publishedValues, indexPath)
            }
            else if( ext == 'json' ) {
                indexPath.text = DumpHelper.prettyPrintJson(outputValue)
            }
            else if( ext == 'yaml' || ext == 'yml' ) {
                indexPath.text = DumpHelper.prettyPrintYaml(outputValue)
            }
            else {
                throw new ScriptRuntimeException("Invalid extension '${ext}' for index file '${indexOpts.path}' -- should be CSV, JSON, or YAML")
            }
            session.notifyFilePublish(new FilePublishEvent(null, indexPath, publishOpts.labels as List))
        }

        log.trace "Completed workflow output '${name}'"
        target.bind(indexPath ?: outputValue)
    }

    /**
     * Extract files from a received value for publishing.
     *
     * @param result
     * @param value
     */
    protected Set<Path> collectFiles(Set<Path> result, value) {
        if( value instanceof Path ) {
            if( shouldPublish(value) )
                result << value
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
     * @param targets
     */
    protected Object normalizeValue(value, Map<Path,Path> targets) {
        if( value instanceof Path ) {
            return normalizePath(value, targets)
        }
        if( value instanceof Collection ) {
            return value.collect { el -> normalizeValue(el, targets) }
        }
        if( value instanceof Map ) {
            return value.collectEntries { k, v -> [k, normalizeValue(v, targets)] }
        }
        return value
    }

    /**
     * Convert a work directory path to the corresponding
     * publish destination.
     *
     * @param path
     * @param targets
     */
    private Path normalizePath(Path path, Map<Path,Path> targets) {
        // a published file is reported by its target path
        final target = targets[path]
        if( target != null )
            return target

        // an unpublished file outside the work directory is still
        // valid, so it is reported as-is
        if( getTaskDir(path) == null )
            return path

        // note: a `path` closure can omit a file in order to not publish it
        return null
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
        String path
        def /* boolean | List<String> */ header = false
        String sep = ','

        IndexOpts(Map opts) {
            this.path = opts.path as String
            if( opts.header != null )
                this.header = opts.header
            if( opts.sep )
                this.sep = opts.sep as String
        }
    }

}
