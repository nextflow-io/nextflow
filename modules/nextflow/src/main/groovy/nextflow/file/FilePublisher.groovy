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

package nextflow.file

import static nextflow.util.CacheHelper.*

import java.nio.file.CopyOption
import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.time.temporal.ChronoUnit
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ExecutorService

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.FilesEx
import nextflow.trace.event.FilePublishEvent
import nextflow.util.HashBuilder
import nextflow.util.RetryConfig
/**
 * Publish files to the workflow output directory.
 *
 * Where {@link nextflow.processor.PublishDir} publishes the contents of a task
 * directory and derives each target from a relative filename, this class publishes
 * an explicit source to target mapping. A workflow output resolves every source file
 * to its own target, so a relative filename is never used to identify a file. Two
 * files with the same name from different task directories are therefore distinct,
 * and a source outside the work directory is publishable like any other.
 *
 * One instance serves a single workflow output for the lifetime of the run.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class FilePublisher {

    enum Mode { SYMLINK, LINK, COPY, MOVE, COPY_NO_FOLLOW, RELLINK }

    private static final List<Mode> LINK_MODES = [Mode.SYMLINK, Mode.LINK, Mode.RELLINK]

    private static final List<Mode> SYMLINK_MODES = [Mode.SYMLINK, Mode.RELLINK]

    private final Session session

    /**
     * Name of the workflow output being published, used in error messages.
     */
    private final String name

    /**
     * The publish mode. When null, it is inferred from the source and target.
     */
    private Mode mode

    /**
     * Whether to overwrite an existing target. Either a boolean or the name of a
     * hash mode, in which case the target is overwritten when its content differs.
     */
    private def /* Boolean | String */ overwrite

    /**
     * Throw an exception when a file cannot be published.
     */
    private boolean failOnError = SysEnv.getBool('NXF_PUBLISH_FAIL_ON_ERROR', true)

    /**
     * Tags to associate with the target file. Only supported by some object stores.
     */
    private def tags

    /**
     * Content type of the target file, either a MIME type or a boolean to probe it.
     * Only supported by some object stores.
     */
    private def contentType

    /**
     * Storage class of the target file. Only supported by some object stores.
     */
    private String storageClass

    /**
     * Labels to associate with the published file.
     */
    private List<String> labels

    private final RetryConfig retryConfig

    private final Map<Path,Boolean> makeCache = new ConcurrentHashMap<>()

    /**
     * Targets already published by this output, used to report two source files
     * being published to the same target.
     */
    private final Map<Path,Path> publishedTargets = new ConcurrentHashMap<>()

    @Lazy
    private ExecutorService threadPool = { session.publishDirExecutorService() }()

    @CompileDynamic
    FilePublisher(Session session, String name, Map opts) {
        this.session = session
        this.name = name
        this.retryConfig = RetryConfig.config(session.config)

        if( opts.mode )
            this.mode = parseMode(opts.mode)
        if( opts.overwrite != null )
            this.overwrite = opts.overwrite
        if( opts.failOnError != null )
            this.failOnError = Boolean.parseBoolean(opts.failOnError.toString())
        if( opts.tags != null )
            this.tags = opts.tags
        if( opts.contentType instanceof Boolean )
            this.contentType = opts.contentType
        else if( opts.contentType )
            this.contentType = opts.contentType as String
        if( opts.storageClass )
            this.storageClass = opts.storageClass as String
        if( opts.labels != null )
            this.labels = opts.labels as List<String>
    }

    protected Mode parseMode(value) {
        if( value instanceof Mode )
            return (Mode)value
        final str = value.toString()
        if( str == 'copyNoFollow' )
            return Mode.COPY_NO_FOLLOW
        try {
            return str.toUpperCase() as Mode
        }
        catch( IllegalArgumentException e ) {
            throw new ScriptRuntimeException("Invalid publish mode '${str}' for workflow output '${name}'")
        }
    }

    /**
     * Publish a set of files.
     *
     * @param mapping Source file to absolute target path
     */
    void publish(Map<Path,Path> mapping) {
        // check every target before publishing anything, so that a conflict
        // does not leave the output directory partially published
        for( final entry : mapping )
            checkTargetConflict(entry.key, entry.value.normalize())
        for( final entry : mapping )
            publish0(entry.key, entry.value.normalize())
    }

    /**
     * Publish a single file.
     *
     * @param source
     * @param target Absolute target path
     */
    void publish(Path source, Path target) {
        if( source == null || target == null )
            return
        final normalized = target.normalize()
        checkTargetConflict(source, normalized)
        publish0(source, normalized)
    }

    private void publish0(Path source, Path target) {
        final resolved = resolveMode(source, target)
        applyFileAttributes(source, target)

        // links are cheap to create, so make them in the calling thread and
        // hand the slower copy and move operations to the publish thread pool
        if( resolved in LINK_MODES )
            safePublishFile(source, target, resolved)
        else
            threadPool.submit({ safePublishFile(source, target, resolved) } as Runnable)
    }

    /**
     * Report two different source files being published to the same target, which
     * would otherwise silently publish whichever file happens to be written first.
     */
    protected void checkTargetConflict(Path source, Path target) {
        final previous = publishedTargets.putIfAbsent(target, source)
        if( previous != null && previous != source )
            throw new ScriptRuntimeException("Publish target '${target.toUriString()}' for workflow output '${name}' is used by more than one file -- offending files: ${previous.toUriString()}, ${source.toUriString()}")
    }

    /**
     * Determine the publish mode for a given source and target. Links cannot be
     * created across file systems, so they fall back to a copy.
     */
    protected Mode resolveMode(Path source, Path target) {
        final sameFileSystem = source.fileSystem == target.fileSystem
            && target.fileSystem == FileSystems.default
            && !target.toString().startsWith('/fusion/s3/')

        if( sameFileSystem )
            return mode ?: Mode.SYMLINK

        if( !mode )
            return Mode.COPY

        if( mode in LINK_MODES ) {
            log.warn1("Cannot use mode `${mode.toString().toLowerCase()}` to publish files to path: ${target.parent} -- Using mode `copy` instead", firstOnly: true)
            return Mode.COPY
        }

        return mode
    }

    protected void applyFileAttributes(Path source, Path target) {
        if( !(target instanceof TagAwareFile) )
            return
        final file = (TagAwareFile)target
        if( tags != null )
            file.setTags(resolveTags(tags))
        if( contentType )
            file.setContentType(contentType instanceof Boolean
                ? Files.probeContentType(source)
                : contentType.toString())
        if( storageClass )
            file.setStorageClass(storageClass)
    }

    protected Map<String,String> resolveTags(tags) {
        final result = tags instanceof Closure ? tags.call() : tags
        if( result instanceof Map<String,String> )
            return result
        throw new ScriptRuntimeException("Invalid publish tags for workflow output '${name}': ${tags}")
    }

    protected void safePublishFile(Path source, Path target, Mode mode) {
        try {
            retryablePublishFile(source, target, mode)
        }
        catch( Throwable e ) {
            final msg = "Failed to publish file: ${source.toUriString()}; to: ${target.toUriString()} [${mode.toString().toLowerCase()}] -- See log file for details"
            if( NF.strictMode || failOnError ) {
                log.error(msg, e)
                session.abort(e)
            }
            else {
                log.warn(msg, e)
            }
        }
    }

    protected void retryablePublishFile(Path source, Path target, Mode mode) {
        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                log.debug "Failed to publish file: ${source.toUriString()}; to: ${target.toUriString()} [${mode.toString().toLowerCase()}] -- attempt: ${event.attemptCount}; reason: ${event.lastException.message}"
            }
        }
        final retryPolicy = RetryPolicy.builder()
            .handle(Exception)
            .withBackoff(retryConfig.delay.toMillis(), retryConfig.maxDelay.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(retryConfig.maxAttempts)
            .withJitter(retryConfig.jitter)
            .onRetry(listener)
            .build()
        Failsafe
            .with(retryPolicy)
            .get({ it -> publishFile(source, target, mode) })
    }

    protected void publishFile(Path source, Path target, Mode mode) {
        makeDirs(target.parent)

        try {
            writeFile(source, target, mode)
        }
        catch( FileAlreadyExistsException e ) {
            // don't publish the source if the target already resolves to it,
            // but still report it as published
            final sameRealPath = checkIsSameRealPath(source, target, mode)

            if( !sameRealPath && shouldOverwrite(source, target, mode) ) {
                FileHelper.deletePath(target)
                writeFile(source, target, mode)
            }
        }

        session.notifyFilePublish(new FilePublishEvent(source, target, labels))
    }

    protected void writeFile(Path source, Path target, Mode mode) {
        log.trace "Publishing file: $source -[$mode]-> $target"

        if( mode == Mode.SYMLINK )
            Files.createSymbolicLink(target, source)
        else if( mode == Mode.RELLINK )
            Files.createSymbolicLink(target, target.parent.relativize(source))
        else if( mode == Mode.LINK )
            FilesEx.mklink(source, [hard: true], target)
        else if( mode == Mode.MOVE )
            FileHelper.movePath(source, target, copyOpts())
        else if( mode == Mode.COPY )
            FileHelper.copyPath(source, target, copyOpts())
        else if( mode == Mode.COPY_NO_FOLLOW )
            FileHelper.copyPath(source, target, copyOpts(LinkOption.NOFOLLOW_LINKS))
        else
            throw new IllegalArgumentException("Unknown file publish mode: ${mode}")
    }

    protected CopyOption[] copyOpts(CopyOption... opts) {
        final copyAttributes = session.config.navigate('workflow.output.copyAttributes', false)
        return copyAttributes
            ? opts + StandardCopyOption.COPY_ATTRIBUTES
            : opts
    }

    protected boolean checkIsSameRealPath(Path source, Path target, Mode mode) {
        if( mode !in SYMLINK_MODES || source.fileSystem != target.fileSystem )
            return false
        final result = realPath(source) == realPath(target)
        if( result )
            log.trace "Skipping publish since source and target real paths are the same - target=$target"
        return result
    }

    protected boolean shouldOverwrite(Path source, Path target, Mode mode) {
        if( overwrite instanceof Boolean )
            return overwrite

        // overwrite when the existing file was published with a different mode
        // (e.g. a symlink but now `copy`, or vice versa), even if its content matches.
        // placed after the explicit `overwrite` check so `overwrite false` is honored
        if( checkPublishModeMismatch(target, mode) )
            return true

        final hashMode = HashMode.of(overwrite) ?: HashMode.DEFAULT()
        final sourceHash = HashBuilder.hashPath(source, source.parent, hashMode)
        final targetHash = HashBuilder.hashPath(target, target.parent, hashMode)
        log.trace "Comparing source and target with mode=${overwrite}, source=${sourceHash}, target=${targetHash}, should overwrite=${sourceHash != targetHash}"
        return sourceHash != targetHash
    }

    protected boolean checkPublishModeMismatch(Path target, Mode mode) {
        if( target.fileSystem != FileSystems.default )
            return false
        return Files.isSymbolicLink(target) != (mode in SYMLINK_MODES)
    }

    private String realPath(Path path) {
        try {
            return path.fileSystem == FileSystems.default
                ? path.toRealPath().toString()
                : path.toUriString()
        }
        catch( NoSuchFileException e ) {
            return path.toString()
        }
        catch( Exception e ) {
            log.warn "Unable to determine real path for '$path'"
            return path.toString()
        }
    }

    protected void makeDirs(Path dir) {
        // nameCount==0 means a filesystem root e.g. an S3 bucket, which
        // always exists and cannot be created
        if( !dir || dir.nameCount == 0 || makeCache.containsKey(dir) )
            return

        try {
            Files.createDirectories(dir)
        }
        catch( FileAlreadyExistsException e ) {
            // ignore
        }
        finally {
            makeCache.put(dir, true)
        }
    }

}
