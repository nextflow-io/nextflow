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

package nextflow.processor

import static nextflow.util.CacheHelper.*

import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileSystem
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.StandardCopyOption
import java.time.temporal.ChronoUnit
import java.util.concurrent.ExecutorService

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.NF
import nextflow.Session
import nextflow.SysEnv
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.file.TagAwareFile
import nextflow.util.HashBuilder
import nextflow.util.PathTrie
/**
 * Implements the {@code publishDir} directory. It create links or copies the output
 * files of a given task to a user specified directory.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString
@EqualsAndHashCode
@CompileStatic
class PublishDir {

    enum Mode { SYMLINK, LINK, COPY, MOVE, COPY_NO_FOLLOW, RELLINK }

    private Map<Path,Boolean> makeCache = new HashMap<>()

    private Session session = Global.session as Session

    /**
     * The target path where create the links or copy the output files
     */
    Path path

    /**
     * Whether to overwrite existing files
     */
    def /* Boolean | String */ overwrite

    /**
     * The publish {@link Mode}
     */
    Mode mode

    /**
     * A glob file pattern to filter the files to be published
     */
    String pattern

    /**
     * SaveAs closure. Allows the dynamically definition of published file names
     */
    Closure saveAs

    /**
     * Enable disable publish rule
     */
    boolean enabled = true

    /**
     * Throw an exception in case publish fails
     */
    boolean failOnError = SysEnv.getBool('NXF_PUBLISH_FAIL_ON_ERROR', true)

    /**
     * Tags to be associated to the target file
     */
    private def tags

    /**
     * The content type of the file. Currently only supported by AWS S3.
     * This can be either a MIME type content type string or a Boolean value
     */
    private contentType

    /**
     * The storage class to be used for the target file.
     * Currently only supported by AWS S3.
     */
    private String storageClass

    private PublishRetryConfig retryConfig

    private PathMatcher matcher

    private FileSystem sourceFileSystem

    private Path sourceDir

    private String stageInMode

    private boolean nullPathWarn

    private TaskRun task

    @Lazy
    private ExecutorService threadPool = { session.publishDirExecutorService() }()

    protected String getTaskName() {
        return task?.getName()
    }

    void setPath( def value ) {
        final resolved = value instanceof Closure ? value.call() : value
        if( resolved instanceof String || resolved instanceof GString )
            nullPathWarn = checkNull(resolved.toString())
        this.path = FileHelper.toCanonicalPath(resolved)
    }

    void setMode( String str ) {
        this.mode = str == 'copyNoFollow' ? Mode.COPY_NO_FOLLOW : str.toUpperCase() as Mode
    }

    void setMode( Mode mode )  {
        this.mode = mode
    }

    static @PackageScope Map<String,String> resolveTags( tags ) {
        def result = tags instanceof Closure
                ? tags.call()
                : tags

        if( result instanceof Map<String,String> )
            return result

        throw new IllegalArgumentException("Invalid publishDir tags attribute: $tags")
    }

    @PackageScope boolean checkNull(String str) {
        ( str =~ /\bnull\b/  ).find()
    }

    /**
     * Object factory method
     *
     * @param params
     *      When the {@code obj} is a {@link Path} or a {@link String} object it is
     *      interpreted as the target path. Otherwise a {@link Map} object matching the class properties
     *      can be specified.
     *
     * @return An instance of {@link PublishDir} class
     */
    @CompileDynamic
    static PublishDir create( Map params ) {
        assert params

        def result = new PublishDir()
        if( params.path )
            result.path = params.path

        if( params.mode )
            result.mode = params.mode

        if( params.pattern )
            result.pattern = params.pattern

        if( params.overwrite != null )
            result.overwrite = params.overwrite

        if( params.saveAs )
            result.saveAs = (Closure) params.saveAs

        if( params.enabled != null )
            result.enabled = Boolean.parseBoolean(params.enabled.toString())

        if( params.failOnError != null )
            result.failOnError = Boolean.parseBoolean(params.failOnError.toString())

        if( params.tags != null )
            result.tags = params.tags

        if( params.contentType instanceof Boolean )
            result.contentType = params.contentType
        else if( params.contentType )
            result.contentType = params.contentType as String

        if( params.storageClass )
            result.storageClass = params.storageClass as String

        return result
    }

    protected Map getRetryOpts() {
        def result = session.config.navigate('nextflow.publish.retryPolicy') as Map
        if( result != null )
            log.warn 'The `nextflow.publish` config scope has been renamed to `workflow.output`'
        else
            result = session.config.navigate('workflow.output.retryPolicy') as Map ?: Collections.emptyMap()
        return result
    }

    protected void apply0(Set<Path> files) {
        assert path
        // setup the retry policy config to be used
        this.retryConfig = new PublishRetryConfig(getRetryOpts())

        createPublishDir()
        validatePublishMode()

        /*
         * when the publishing is using links, create them in process
         * otherwise copy and moving file can take a lot of time, thus
         * apply the operation using an external thread
         */
        final inProcess = mode == Mode.LINK || mode == Mode.SYMLINK || mode == Mode.RELLINK

        if( pattern ) {
            this.matcher = FileHelper.getPathMatcherFor("glob:${pattern}", sourceFileSystem)
        }

        /*
         * iterate over the file parameter and publish each single file
         */
        for( Path value : dedupPaths(files) ) {
            apply1(value, inProcess)
        }
    }

    /**
     * Find out only path not overlapping each other using prefix tree.
     * This is require to avoid copy multiple times the same files, when
     * the output declares both directory and files nested in the same directory.
     *
     * See also https://github.com/nextflow-io/nextflow/issues/2177
     *
     * @param files A collection of files. NOTE: MUST be local files. Remote file scheme e.g. S# is not supported
     * @return
     */
    protected List<Path> dedupPaths(Collection<Path> files) {
        if( !files )
            return Collections.emptyList()
        final trie = new PathTrie()
        for( Path it : files )
            trie.add(it)
        // convert to paths
        final result = new ArrayList()
        for( String it : trie.traverse() ) {
            result.add( FileHelper.asPath(it) )
        }
        return result
    }

    void apply( Set<Path> files, Path sourceDir ) {
        if( !files || !enabled )
            return
        this.sourceDir = sourceDir
        this.sourceFileSystem = sourceDir ? sourceDir.fileSystem : null
        apply0(files)
    }

    /**
     * Apply the publishing process to the specified {@link TaskRun} instance
     *
     * @param files Set of output files
     * @param task The task whose output need to be published
     */
    void apply( Set<Path> files, TaskRun task ) {

        if( !files || !enabled )
            return

        if( !path )
            throw new IllegalStateException("Target path for directive publishDir cannot be null")

        if( nullPathWarn )
            log.warn "Process `$task.processor.name` publishDir path contains a variable with a null value"

        this.sourceDir = task.targetDir
        this.sourceFileSystem = sourceDir.fileSystem
        this.stageInMode = task.config.stageInMode
        this.task = task

        apply0(files)
    }

    protected void apply1(Path source, boolean inProcess ) {

        def target = sourceDir ? sourceDir.relativize(source) : source.getFileName()
        if( matcher && !matcher.matches(target) ) {
            // skip not matching file
            return
        }

        if( saveAs && !(target=saveAs.call(target.toString()))) {
            // skip this file
            return
        }

        final destination = resolveDestination(target).normalize()

        // apply tags
        if( this.tags!=null && destination instanceof TagAwareFile ) {
            destination.setTags( resolveTags(this.tags) )
        }
        // apply content type
        if( contentType && destination instanceof TagAwareFile ) {
            final String type = this.contentType instanceof Boolean
                    ? Files.probeContentType(source)
                    : this.contentType.toString()
            destination.setContentType(type)
        }
        // storage class
        if( storageClass && destination instanceof TagAwareFile ) {
            destination.setStorageClass(storageClass)
        }

        if( inProcess ) {
            safeProcessFile(source, destination)
        }
        else {
            threadPool.submit({ safeProcessFile(source, destination) } as Runnable)
        }

    }

    protected Path resolveDestination(target) {

        if( target instanceof Path ) {
            if( target.isAbsolute() ) {
                return (Path)target
            }
            // note: convert to a string to avoid `ProviderMismatchException` when the
            // destination `path` is not a unix file system
            return path.resolve(target.toString())
        }

        if( target instanceof CharSequence ) {
            return path.resolve(target.toString())
        }

        throw new IllegalArgumentException("Not a valid publish target path: `$target` [${target?.class?.name}]")
    }

    protected void safeProcessFile(Path source, Path target) {
        try {
            retryableProcessFile(source, target)
        }
        catch( Throwable e ) {
            final msg =  "Failed to publish file: ${source.toUriString()}; to: ${target.toUriString()} [${mode.toString().toLowerCase()}] -- See log file for details"
            final shouldFail = NF.strictMode || failOnError
            if( shouldFail ) {
                log.error(msg,e)
                session?.abort(e)
            }
            else
                log.warn(msg,e)
        }
    }

    protected void retryableProcessFile(Path source, Path target) {
        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                log.debug "Failed to publish file: ${source.toUriString()}; to: ${target.toUriString()} [${mode.toString().toLowerCase()}] -- attempt: ${event.attemptCount}; reason: ${event.lastFailure.message}"
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
            .with( retryPolicy )
            .get({it-> processFile(source, target)})
    }

    protected void processFile( Path source, Path destination ) {

        // create target dirs if required
        makeDirs(destination.parent)

        try {
            processFileImpl(source, destination)
        }
        catch( FileAlreadyExistsException e ) {
            // don't copy source path if target is identical, but still emit the publish event
            // see also https://github.com/nextflow-io/nf-prov/issues/22
            final sameRealPath = checkIsSameRealPath(source, destination)

            // make sure destination and source does not overlap
            // see https://github.com/nextflow-io/nextflow/issues/2177
            if( !sameRealPath && checkSourcePathConflicts(destination))
                return
            
            if( !sameRealPath && shouldOverwrite(source, destination) ) {
                FileHelper.deletePath(destination)
                processFileImpl(source, destination)
            }
        }

        notifyFilePublish(destination, source)
    }

    private String real0(Path p) {
        try {
            // resolve symlink if it's file in the default (posix) file system
            return p.fileSystem == FileSystems.default
                    ? p.toRealPath().toString()
                    : p.toUriString()
        }
        catch (NoSuchFileException e) {
            return p.toString()
        }
        catch (Exception e) {
            log.warn "Unable to determine real path for '$p'"
            return p.toString()
        }
    }

    protected boolean checkIsSameRealPath(Path source, Path target) {
        if( !isSymlinkMode() || source.fileSystem!=target.fileSystem )
            return false

        final t1 = real0(target)
        final s1 = real0(source)
        final result = s1 == t1
        log.trace "Skipping publishDir since source and target real paths are the same - target=$target; real=$t1"
        return result
    }

    protected boolean checkSourcePathConflicts(Path target) {
        if( !isSymlinkMode() || sourceDir.fileSystem!=target.fileSystem )
            return false

        final t1 = real0(target)
        final s1 = real0(sourceDir)
        if( t1.startsWith(s1) ) {
            def msg = "Refusing to publish file since destination path conflicts with the task work directory!"
            def name0 = getTaskName()
            if( name0 )
                msg += "\n- offending task  : $name0"
            msg += "\n- offending file  : $target"
            if( t1 != target.toString() )
                msg += "\n- real destination: $t1"
            msg += "\n- task directory  : $s1"
            log.warn1(msg)
            return true
        }
        return false
    }

    protected boolean isSymlinkMode() {
        return !mode || mode == Mode.SYMLINK || mode == Mode.RELLINK
    }

    protected boolean shouldOverwrite(Path source, Path target) {
        if( overwrite instanceof Boolean )
            return overwrite

        final hashMode = HashMode.of(overwrite) ?: HashMode.DEFAULT()
        final sourceHash = HashBuilder.hashPath(source, source.parent, hashMode)
        final targetHash = HashBuilder.hashPath(target, target.parent, hashMode)
        log.trace "comparing source and target with mode=${overwrite}, source=${sourceHash}, target=${targetHash}, should overwrite=${sourceHash != targetHash}"
        return sourceHash != targetHash
    }

    protected void processFileImpl( Path source, Path destination ) {
        log.trace "publishing file: $source -[$mode]-> $destination"

        if( !mode || mode == Mode.SYMLINK ) {
            Files.createSymbolicLink(destination, source)
        }
        else if( mode == Mode.RELLINK ) {
            def sourceRelative = destination.getParent().relativize(source)
            Files.createSymbolicLink(destination, sourceRelative)
        }
        else if( mode == Mode.LINK ) {
            FilesEx.mklink(source, [hard:true], destination)
        }
        else if( mode == Mode.MOVE ) {
            FileHelper.movePath(source, destination, StandardCopyOption.COPY_ATTRIBUTES)
        }
        else if( mode == Mode.COPY ) {
            FileHelper.copyPath(source, destination, StandardCopyOption.COPY_ATTRIBUTES)
        }
        else if( mode == Mode.COPY_NO_FOLLOW ) {
            FileHelper.copyPath(source, destination, StandardCopyOption.COPY_ATTRIBUTES, LinkOption.NOFOLLOW_LINKS)
        }
        else {
            throw new IllegalArgumentException("Unknown file publish mode: ${mode}")
        }
    }

    protected void createPublishDir() {
        try {
            makeDirs(path)
        }
        catch( Throwable e ) {
            throw new IllegalStateException("Failed to create publish directory: ${path.toUriString()}", e)
        }
    }

    protected void makeDirs(Path dir) {
        if( !dir || makeCache.containsKey(dir) )
            return

        try {
            Files.createDirectories(dir)
        }
        catch( FileAlreadyExistsException e ) {
            // ignore
        }
        finally {
            makeCache.put(dir,true)
        }
    }

    /*
     * That valid publish mode has been selected
     * Note: link and symlinks are not allowed across different file system
     */
    @PackageScope
    void validatePublishMode() {
        if( log.isTraceEnabled() )
            log.trace "Publish path: ${path.toUriString()}; notMatchSourceFs=${sourceFileSystem && sourceFileSystem != path.fileSystem}; notMatchDefaultFs=${path.fileSystem != FileSystems.default}; isFusionFs=${path.toString().startsWith('/fusion/s3/')}"
        if( (sourceFileSystem && sourceFileSystem != path.fileSystem) || path.fileSystem != FileSystems.default || path.toString().startsWith('/fusion/s3/') ) {
            if( !mode ) {
                mode = Mode.COPY
            }
            else if( mode == Mode.SYMLINK || mode == Mode.LINK || mode == Mode.RELLINK ) {
                log.warn1("Cannot use mode `${mode.toString().toLowerCase()}` to publish files to path: $path -- Using mode `copy` instead", firstOnly:true)
                mode = Mode.COPY
            }
        }

        if( !mode ) {
            mode = stageInMode=='rellink' ? Mode.RELLINK : Mode.SYMLINK
        }
    }

    protected void notifyFilePublish(Path destination, Path source=null) {
        session.notifyFilePublish(destination, source)
    }

}
