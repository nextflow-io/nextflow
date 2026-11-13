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

package nextflow.processor

import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.PathMatcher

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.file.FileHelper
import nextflow.file.FilePublisher
import nextflow.util.PathTrie
/**
 * Implements the {@code publishDir} directory. It create links or copies the output
 * files of a given task to a user specified directory.
 *
 * Each file is mapped to its target path here, and published by a {@link FilePublisher}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString
@EqualsAndHashCode
@CompileStatic
class PublishDir {

    /**
     * Kept for plugins that use the publish mode of a `publishDir` directive.
     * Each value maps to the {@link FilePublisher.Mode} of the same name.
     */
    enum Mode { SYMLINK, LINK, COPY, MOVE, COPY_NO_FOLLOW, RELLINK }

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
     * Labels to be associated to the target file
     */
    private List<String> labels

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

    private PathMatcher matcher

    private FileSystem sourceFileSystem

    private Path sourceDir

    private String stageInMode

    private boolean nullPathWarn

    private TaskRun task

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
        this.mode = Mode.valueOf(FilePublisher.parseMode(str).name())
    }

    void setMode( Mode mode )  {
        this.mode = mode
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

        if( params.labels != null )
            result.labels = params.labels as List<String>

        if( params.contentType instanceof Boolean )
            result.contentType = params.contentType
        else if( params.contentType )
            result.contentType = params.contentType as String

        if( params.storageClass )
            result.storageClass = params.storageClass as String

        return result
    }

    protected void apply0(Set<Path> files) {
        assert path

        final publisher = createPublisher()
        createPublishDir(publisher)

        if( pattern ) {
            this.matcher = FileHelper.getPathMatcherFor("glob:${pattern}", sourceFileSystem)
        }

        /*
         * iterate over the file parameter and publish each single file
         */
        for( Path value : dedupPaths(files) ) {
            apply1(publisher, value)
        }
    }

    protected FilePublisher createPublisher() {
        final opts = [
            mode: mode ? FilePublisher.Mode.valueOf(mode.name()) : null,
            defaultMode: stageInMode == 'rellink' ? FilePublisher.Mode.RELLINK : null,
            overwrite: overwrite,
            failOnError: failOnError,
            tags: tags,
            contentType: contentType,
            storageClass: storageClass,
            labels: labels,
            sourceDir: sourceDir,
            taskName: getTaskName()
        ]
        return new FilePublisher(session, opts)
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

    protected void apply1(FilePublisher publisher, Path source) {

        def target = sourceDir ? sourceDir.relativize(source) : source.getFileName()
        if( matcher && !matcher.matches(target) ) {
            // skip not matching file
            return
        }

        if( saveAs && !(target=saveAs.call(target.toString()))) {
            // skip this file
            return
        }

        publisher.publish(source, resolveDestination(target))
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

    protected void createPublishDir(FilePublisher publisher) {
        try {
            publisher.makeDirs(path)
        }
        catch( Throwable e ) {
            throw new IllegalStateException("Failed to create publish directory: ${path.toUriString()}", e)
        }
    }

}
