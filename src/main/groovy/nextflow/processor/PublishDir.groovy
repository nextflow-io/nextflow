/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileSystem
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.util.concurrent.ExecutorService
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
/**
 * Implements the {@code publishDir} directory. It create links or copies the output
 * files of a given task to a user specified directory.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString
@EqualsAndHashCode
class PublishDir {

    enum Mode { SYMLINK, LINK, COPY, MOVE, COPY_NO_FOLLOW, RELLINK }

    private Map<Path,Boolean> makeCache = new HashMap<>()

    /**
     * The target path where create the links or copy the output files
     */
    Path path

    /**
     * Whenever overwrite existing files
     */
    Boolean overwrite

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

    private PathMatcher matcher

    private FileSystem sourceFileSystem

    private TaskProcessor processor

    private Path sourceDir

    private String stageInMode

    void setPath( Closure obj ) {
        setPath( obj.call() as Path )
    }

    void setPath( String str ) {
        setPath(str as Path)
    }

    void setPath( Path obj ) {
        this.path = obj.complete()
    }

    void setMode( String str ) {
        this.mode = str == 'copyNoFollow' ? Mode.COPY_NO_FOLLOW : str.toUpperCase() as Mode
    }

    void setMode( Mode mode )  {
        this.mode = mode
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
            result.saveAs = params.saveAs

        return result
    }

    /**
     * Apply the publishing process to the specified {@link TaskRun} instance
     *
     * @param task The task whose output need to be published
     */
    @CompileStatic
    void apply( List<Path> files, TaskRun task ) {

        if( !files ) {
            return
        }

        this.processor = task.processor
        this.sourceDir = task.targetDir
        this.sourceFileSystem = sourceDir.fileSystem
        this.stageInMode = task.config.stageInMode
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
        files.each { value ->
            apply(value, inProcess)
        }
    }


    @CompileStatic
    protected void apply( Path source, boolean inProcess ) {

        def target = sourceDir.relativize(source)
        if( matcher && !matcher.matches(target) ) {
            // skip not matching file
            return
        }

        if( saveAs && !(target=saveAs.call(target.toString()))) {
            // skip this file
            return
        }

        final destination = resolveDestination(target)
        if( inProcess ) {
            safeProcessFile(source, destination)
        }
        else {
            executor.submit({ safeProcessFile(source, destination) } as Runnable)
        }

    }

    @CompileStatic
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

    @CompileStatic
    protected void safeProcessFile(Path source, Path target) {
        try {
            processFile(source, target)
        }
        catch( Throwable e ) {
            log.warn "Failed to publish file: $source; to: $target [${mode.toString().toLowerCase()}] -- See log file for details", e
        }
    }

    @CompileStatic
    protected void processFile( Path source, Path destination ) {

        // create target dirs if required
        makeDirs(destination.parent)

        try {
            processFileImpl(source, destination)
        }
        catch( FileAlreadyExistsException e ) {
            if( !overwrite )
                return

            FileHelper.deletePath(destination)
            processFileImpl(source, destination)
        }
    }

    @CompileStatic
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
            FileHelper.movePath(source, destination)
        }
        else if( mode == Mode.COPY ) {
            FileHelper.copyPath(source, destination)
        }
        else if( mode == Mode.COPY_NO_FOLLOW ) {
            FileHelper.copyPath(source, destination, LinkOption.NOFOLLOW_LINKS)
        }
        else {
            throw new IllegalArgumentException("Unknown file publish mode: ${mode}")
        }
    }

    @CompileStatic
    private void createPublishDir() {
        makeDirs(this.path)
    }

    @CompileStatic
    private void makeDirs(Path dir) {
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
    @CompileStatic
    @PackageScope
    void validatePublishMode() {

        if( sourceFileSystem != path.fileSystem || path.fileSystem != FileSystems.default ) {
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

    @Memoized // <-- this guarantees that the same executor is used across different publish dir in the same session
    @CompileStatic
    static synchronized ExecutorService createExecutor(Session session) {
        final result = new ThreadPoolExecutor(0, Runtime.runtime.availableProcessors(),
                60L, TimeUnit.SECONDS,
                new LinkedBlockingQueue<Runnable>());

        session?.onShutdown {
            result.shutdown()
            result.awaitTermination(36,TimeUnit.HOURS)
        }

        return result
    }

    @PackageScope
    static ExecutorService getExecutor() {
        createExecutor(Global.session as Session)
    }
}
