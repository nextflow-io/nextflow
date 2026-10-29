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


package nextflow.lineage.cleanup

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.util.concurrent.CompletableFuture
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ExecutionException

import com.google.common.hash.HashCode
import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.exception.ProcessException
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.lineage.LinStore
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskResolver
import nextflow.processor.TaskRun
import nextflow.trace.DefaultObserverFactory
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.TaskEvent
import nextflow.script.types.Bag
import nextflow.util.ArrayBag
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper
import nextflow.util.HashBuilder
/**
 * Track the output files of resumed tasks that were deleted by
 * eager cleanup, so that they can be hashed from their lineage
 * metadata and restored when a task needs them.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class DeletedOutputs implements TraceObserverV2 {

    private static DeletedOutputs instance

    private Session session

    private LinStore store

    private Map<Path,DeletedFile> files = new ConcurrentHashMap<>()

    private Map<Path,Producer> restoring = new ConcurrentHashMap<>()

    private Map<String,List<FileOutput>> publishedFiles

    private DeletedOutputs(Session session, LinStore store) {
        this.session = session
        this.store = store
    }

    /**
     * Get the deleted outputs for the given session, or null if resume
     * is not supported, i.e. eager cleanup and lineage are not
     * enabled or the run is not resumed.
     *
     * @param session
     */
    static synchronized DeletedOutputs get(Session session) {
        if( session == null )
            return null
        if( instance?.session !== session )
            instance = create(session)
        return instance?.store != null ? instance : null
    }

    /**
     * Get the deleted outputs for the current session.
     */
    static DeletedOutputs current() {
        return get((Session) Global.session)
    }

    private static DeletedOutputs create(Session session) {
        if( session.config?.cleanup != 'eager' || !session.resumeMode || !DefaultObserverFactory.isDefaultCache(session) )
            return new DeletedOutputs(session, null)
        final store = LinStoreFactory.getOrCreate(session)
        if( store == null )
            log.warn "Eager cleanup requires lineage to resume tasks whose outputs were deleted -- set `lineage.enabled = true` to enable it"
        return new DeletedOutputs(session, store)
    }

    /**
     * Register the output files of a task from its lineage, so that
     * the task can be resumed even if some outputs were deleted.
     *
     * @param producer
     * @return the registered output files
     */
    Set<Path> register(Producer producer) {
        final Set<Path> result = []
        try( final keys = store.getSubKeys(producer.key.toString()) ) {
            for( final key : keys.toList() ) {
                final value = store.load(key)
                if( value !instanceof FileOutput )
                    continue
                final output = (FileOutput) value
                final path = FileHelper.asPath(output.path)
                files.put(path, new DeletedFile(path, output.size, output.modifiedAt?.toInstant()?.toEpochMilli(), producer))
                result.add(path)
            }
        }
        catch( Exception e ) {
            log.debug "Unable to load lineage outputs for task hash ${producer.key} -- cause: ${e.message ?: e}"
        }
        producer.outputs = result
        return result
    }

    void unregister(Set<Path> paths) {
        for( final path : paths )
            files.remove(path)
    }

    /**
     * Determine whether a task has any inputs that were registered
     * as deleted outputs.
     *
     * @param task
     */
    boolean hasDeletedInputs(TaskRun task) {
        return task.inputFiles.any( holder -> files.containsKey(holder.sourceObj) )
    }

    /**
     * Replace any deleted files in a task input value with their
     * lineage metadata for hashing.
     *
     * @param value
     */
    Object substitute(Object value) {
        if( value instanceof FileHolder )
            return files.get(value.sourceObj) ?: value
        if( value instanceof Path )
            return files.get(value) ?: value
        if( value instanceof Bag )
            return new ArrayBag(value.collect( el -> substitute(el) ))
        if( value instanceof Set )
            return value.collect(new LinkedHashSet(), el -> substitute(el))
        if( value instanceof Collection )
            return value.collect( el -> substitute(el) )
        if( value instanceof Map )
            return value.collectEntries { k, v -> [k, substitute(v)] }
        return value
    }

    /**
     * Restore any deleted input files of a task, by re-executing the
     * tasks that produced them, and wait until they are restored.
     *
     * @param task
     */
    void restoreInputs(TaskRun task) {
        final missing = missingInputs(task)
        try {
            restoreAll(missing.values()).get()
        }
        catch( ExecutionException e ) {
            throw new ProcessException("Unable to restore input files deleted by eager cleanup for task `${task.lazyName()}`", e.cause)
        }

        // a re-executed task may not reproduce the same file names, e.g. if they contain a timestamp
        for( final entry : missing.entrySet() ) {
            if( !Files.exists(entry.key) )
                throw new ProcessException("Unable to restore input file `${entry.key.fileName}` deleted by eager cleanup for task `${task.lazyName()}` -- task `${entry.value.task.lazyName()}` did not produce it when re-executed")
        }
    }

    private Map<Path,Producer> missingInputs(TaskRun task) {
        final Map<Path,Producer> result = [:]
        for( final holder : task.inputFiles ) {
            final file = files.get(holder.sourceObj)
            if( file != null && !Files.exists(file.path) )
                result.put(file.path, file.producer)
        }
        return result
    }

    private CompletableFuture<Void> restoreAll(Collection<Producer> producers) {
        final futures = new HashSet<Producer>(producers).collect( p -> restore(p) )
        return CompletableFuture.allOf(futures as CompletableFuture[])
    }

    /**
     * Re-execute a task in a temporary work directory to restore its
     * deleted output files.
     *
     * The task is not re-executed in its original work directory, since
     * other tasks might be using the outputs that were not deleted.
     *
     * The task is launched once its own deleted input files are restored,
     * without waiting for them, so that independent restores run in parallel.
     *
     * @param producer
     */
    private CompletableFuture<Void> restore(Producer producer) {
        synchronized (producer) {
            if( producer.future != null && (!producer.future.isDone() || !producer.anyMissing()) )
                return producer.future

            final task = producer.task.makeCopy()
            task.restoreDir = producer.workDir
            final key = HashBuilder.defaultHasher().putBytes(producer.key.asBytes()).putUnencodedChars('restore').hash()
            final workDir = producer.resolver.workDirFor(key)
            FileHelper.deletePath(workDir)
            Files.createDirectories(workDir)

            log.debug "Restoring files deleted by eager cleanup > ${task.lazyName()}"
            final future = new CompletableFuture<Void>()
            producer.future = future
            restoring.put(producer.workDir, producer)
            try {
                restoreAll(missingInputs(task).values())
                    .thenRunAsync( () -> producer.resolver.launch(task, key, workDir) )
                    .exceptionally( (Throwable e) -> {
                        restoring.remove(producer.workDir)
                        future.completeExceptionally(e.cause ?: e)
                        return null
                    } )
            }
            catch( Throwable e ) {
                restoring.remove(producer.workDir)
                future.completeExceptionally(e)
            }
            return future
        }
    }

    /**
     * When a restored task completes, notify the tasks that are
     * waiting for its output files.
     *
     * @param event
     */
    @Override
    void onTaskComplete(TaskEvent event) {
        final task = event.handler.task
        if( task.restoreDir == null )
            return
        if( task.failed && task.errorAction == ErrorStrategy.RETRY )
            return
        final producer = restoring.remove(task.restoreDir)
        if( producer == null )
            return
        if( task.failed || !task.isSuccess() )
            producer.future.completeExceptionally(new ProcessException("Task `${task.lazyName()}` failed to restore its output files"))
        else
            producer.future.complete(null)
    }

    /**
     * Determine whether the deleted outputs of a task are up to date
     * with their published copies, i.e. each deleted output that was
     * published by a previous run of this session is still the source
     * of its published copy. Otherwise the task must be re-executed so
     * that its outputs are published again.
     *
     * @param producer
     */
    boolean isPublishedUpToDate(Producer producer) {
        final published = getPublishedFiles()
        for( final path : producer.outputs ) {
            if( Files.exists(path) )
                continue
            final records = published.get("lid://${producer.key}/${producer.workDir.relativize(path)}".toString())
            if( records && !records.any( record -> isSameFile(record) ) ) {
                log.debug "Published copy of deleted file ${path.toUriString()} is out of date"
                return false
            }
        }
        return true
    }

    private static boolean isSameFile(FileOutput record) {
        final target = FileHelper.asPath(record.path)
        if( !Files.exists(target) )
            return false
        final attrs = Files.readAttributes(target, BasicFileAttributes)
        return record.size == attrs.size()
            && record.modifiedAt?.toInstant()?.toEpochMilli() == attrs.lastModifiedTime().toMillis()
    }

    /**
     * Get the workflow outputs published by previous runs of this
     * session, indexed by the lineage id of the source file.
     */
    private synchronized Map<String,List<FileOutput>> getPublishedFiles() {
        if( publishedFiles != null )
            return publishedFiles

        publishedFiles = [:]
        final runs = store.getHistoryLog().getRecords().findAll( r -> r.sessionId == session.uniqueId && r.runName != session.runName )
        for( final run : runs ) {
            try( final keys = store.getSubKeys(run.runLid.replace('lid://', '')) ) {
                for( final key : keys.toList() ) {
                    final value = store.load(key)
                    if( value instanceof FileOutput && value.source )
                        publishedFiles.computeIfAbsent(value.source, k -> []).add(value)
                }
            }
            catch( Exception e ) {
                log.debug "Unable to load lineage outputs for run ${run.runLid} -- cause: ${e.message ?: e}"
            }
        }
        return publishedFiles
    }

    @PackageScope
    static class Producer {
        TaskRun task
        HashCode key
        Path workDir
        TaskResolver resolver
        Set<Path> outputs = []
        CompletableFuture<Void> future

        Producer(TaskRun task, HashCode key, Path workDir, TaskResolver resolver) {
            this.task = task
            this.key = key
            this.workDir = workDir
            this.resolver = resolver
        }

        boolean anyMissing() {
            return outputs.any( path -> !Files.exists(path) )
        }
    }

    /**
     * A deleted file which is hashed from its lineage metadata in the
     * same way as the original file.
     */
    @PackageScope
    static class DeletedFile implements CacheFunnel {
        Path path
        long size
        Long lastModified
        Producer producer

        DeletedFile(Path path, long size, Long lastModified, Producer producer) {
            this.path = path
            this.size = size
            this.lastModified = lastModified
            this.producer = producer
        }

        @Override
        Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
            hasher.putUnencodedChars(path.toAbsolutePath().toString())
            hasher.putLong(size)
            if( lastModified != null && mode != CacheHelper.HashMode.LENIENT )
                hasher.putLong(lastModified)
            return hasher
        }
    }

}
