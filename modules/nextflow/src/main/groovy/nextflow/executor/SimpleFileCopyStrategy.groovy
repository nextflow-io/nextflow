/*
 * Copyright 2020-2021, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.util.Escape
import static nextflow.util.SpuriousDeps.getS3UploaderScript

/**
 * Simple file strategy that stages input files creating symlinks
 * and copies the output files using the {@code cp} command.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SimpleFileCopyStrategy implements ScriptFileCopyStrategy {

    /**
     * Specify how output files are copied to the target path
     */
    String stageoutMode

    /**
     * Specifies how input files are copied to the work path
     */
    String stageinMode

    /**
     * File names separator
     */
    String separatorChar = '\n'

    /**
     * Target output folder
     */
    Path targetDir

    /**
     * Task work directory
     */
    Path workDir

    SimpleFileCopyStrategy() { }

    SimpleFileCopyStrategy( TaskBean bean ) {
        this.stageinMode = bean.stageInMode
        this.stageoutMode = bean.stageOutMode
        this.targetDir = bean.targetDir
        this.workDir = bean.workDir
    }

    /**
     * Copies foreign files in the pipeline work directory so that
     * those files can be accessed as regular local files.
     *
     * Foreign files are all files that are hosted in a storage remote
     * to the current work directory file system, therefore when
     * using a remote bucket as work directory local files are
     * *uploaded* to the remote bucket work directory.
     *
     * WARNING: A file upload race condition may arise if multiple
     * concurrent nextflow instances use the same pipeline work directory.
     * Therefore concurrent runs are expected to use different work
     * directory when automatic remote input file staging feature is used.
     *
     * @param files
     *      The map on input file, holds ( local file name, storage path ) pairs
     * @return
     *      A map containing remote file paths resolved as local paths 
     */
//    @Override
//    Map<String,Path> resolveForeignFiles(Map<String,Path> files) {
//        if( !files  )
//            return files
//        if( !porter )
//            porter = ((Session)Global.session).getFilePorter()
//        porter.stageForeignFiles(files, getStagingDir())
//    }

    protected Path getStagingDir() {
        final result = workDir.parent?.parent?.resolve('stage')
        if( !result )
            throw new IllegalArgumentException("Cannot resolve staging directory for task work dir: ${workDir.toUriString()}")
        return result
    }

    /**
     * Given a map of the input file parameters with respective values,
     * create the BASH script to stage them into the task working space
     *
     * @param inputFiles All the input files as a map of <stage name, store path> pairs
     * @return The BASH script to stage them
     */
    String getStageInputFilesScript(Map<String,Path> inputFiles) {
        assert inputFiles != null

        def delete = []
        def links = []
        for( Map.Entry<String,Path> entry : inputFiles ) {
            final stageName = entry.key
            final storePath = entry.value

            // delete all previous files with the same name
            delete << "rm -f ${Escape.path(stageName)}"

            // link them
            links << stageInputFile( storePath, stageName )
        }

        // return a big string containing the command
        return (delete + links).join(separatorChar)
    }


    /**
     * Stage the input file into the task working area. By default it creates a symlink
     * to the the specified path using {@code targetName} as name.
     * <p>
     *     An executor may override it to support a different staging mechanism
     *
     * @param path The {@code Path} to the file to be staged
     * @param targetName The name to be used in the task working directory
     * @return The script which will apply the staging for this file in the main script
     */
    String stageInputFile( Path path, String targetName ) {
        def cmd = ''
        def p = targetName.lastIndexOf('/')
        if( p>0 ) {
            cmd += "mkdir -p ${Escape.path(targetName.substring(0,p))} && "
        }
        cmd += stageInCommand(path.toAbsolutePath().toString(), targetName, stageinMode)

        return cmd
    }

    /**
     * Creates the script to unstage the result output files from the scratch directory
     * to the shared working directory
     */
    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {

        // collect all the expected names (pattern) for files to be un-staged
        def result = []

        // create a bash script that will copy the out file to the working directory
        log.trace "Unstaging file path: $outputFiles"
        if( outputFiles ) {
            final mode = stageoutMode ?: ( workDir==targetDir ? 'copy' : 'move' )
            def prefix = getUnstagePrefix(targetDir)
            if( prefix )
                result << prefix
            result << stageOutCommand(outputFiles, targetDir, mode)
        }

        return result.join(separatorChar)
    }

    /**
     * Prepare a bash command for staging input files for a task.
     *
     * This method should only be called when workDir is defined, and with a
     * relative target path.
     *
     * NB: 'target' and 'source' are the opposite from the ln command
     *
     * Methods:
     * * symlink: Create a symbolic link to source from target using an
     *   absolute path.
     * * rellink: Create a symbolic link to source from target using
     *   a relative path.
     * * link: Create a hard link to source from target.
     * * copy: Copy the file at source to target.
     * 
     * @param source The original file that is to be staged.
     * @param target The new path to create.
     * @param mode The method to use for staging the file. See Methods above.
     * @return The command to issue to stage the specified file.
     */ 
    protected String stageInCommand( String source, String target, String mode ) {
        if( !target || target.startsWith('/') )
            throw new IllegalArgumentException("Process input file target path must be relative: $target")

        if( mode == 'symlink' || !mode )
            return "ln -s ${Escape.path(source)} ${Escape.path(target)}"

        if( mode == 'rellink' ) {
            // GNU ln has the '-r' flag, but BSD ln doesn't, so we have to
            // manually resolve the relative path.

            def targetPath = workDir.resolve(target)
            def sourcePath = workDir.resolve(source)
            source = targetPath.getParent().relativize(sourcePath).toString()

            return "ln -s ${Escape.path(source)} ${Escape.path(target)}"
        }

        if( mode == 'link' )
            return "ln ${Escape.path(source)} ${Escape.path(target)}"

        if( mode == 'copy' || !mode )
            return "cp -fRL ${Escape.path(source)} ${Escape.path(target)}"

        throw new IllegalArgumentException("Unknown stage-in strategy: $mode")
    }

    protected String getPathScheme(Path path) {
        path?.getFileSystem()?.provider()?.getScheme()
    }

    protected String getUnstagePrefix(Path targetDir) {
        final scheme = getPathScheme(targetDir)
        if( scheme == 'file' ) {
            return "mkdir -p ${Escape.path(targetDir)}"
        }
        else
            return null
    }

    /**
     * Compose a `cp` or `mv` command given the source and target paths
     *
     * @param source
     * @param targetDir
     * @param move When {@code true} use unix {@code mv} instead of {@code cp}
     * @return A shell copy or move command string
     */

    protected String stageOutCommand( List<String> source, Path targetDir, String mode ) {
        def scheme = getPathScheme(targetDir)
        if( scheme == 'file' )
            return stageOutCommand(source, targetDir.toString(), mode)

        if( scheme == 's3' )
            return source
                    .collect {"nxf_s3_upload '$it' s3:/$targetDir || true" }
                    .join(separatorChar)

        throw new IllegalArgumentException("Unsupported target path: ${targetDir.toUriString()}")
    }

    protected String stageOutCommand( List<String> source, String target, String mode ) {
        assert mode

        String cmd
        if( mode == 'copy' )
            cmd = 'cp -fRL'
        else if( mode == 'move' )
            cmd = 'mv'
        else if( mode == 'rsync' )
            //This will not work for glob terms
            return "rsync -rRl ${normalizeGlobStarPaths( source ).collect{ Escape.path( it ) }.join(' ')} ${Escape.path(target)} || true"
        else
            throw new IllegalArgumentException("Unknown stage-out strategy: $mode")

        cmd = "sh -c 'mkdir -p \"${Escape.path(target)}/`dirname \\\"\$1\\\"`\"; $cmd \"\$1\" \"${Escape.path(target)}/`dirname \\\"\$1\\\"`\";' _ \"\$name\""

        final List<String> escape = source
                .collect {
                    String escaped = Escape.path( it, true )
                    //Bash glob behaves different than Java's Glob, if the path starts with **
                    //https://unix.stackexchange.com/questions/49913/recursive-glob
                    if( escaped.startsWith("**") && escaped.size() > 2 ) {
                        escaped = "**/*" + escaped.substring(2)
                    }
                    return escaped
                }

        return """\
            shopt -s globstar extglob || true
            IFS=\$'\\n'
            pathes=`ls -1d ${escape.join(' ')} | sort | uniq`
            set -f
            for name in \$pathes; do
                if [[ ! -e "${Escape.path(target)}/\$name" ]]; then
                    $cmd || true
                fi
            done
            set +f
            shopt -u globstar extglob || true
            unset IFS""".stripIndent(true)
    }

    /**
     * Normalize path that contains glob star wildcards (i.e. double star character) since
     * they are not supported by plain BASH
     *
     * @param files
     * @return
     */
    protected List<String> normalizeGlobStarPaths( List<String> files ) {

        def result = new ArrayList(files.size())
        for( int i=0; i<files.size(); i++ ) {
            def item = removeGlobStar(files.get(i))
            if( !result.contains(item) )
                result.add(item)
        }

        return result
    }

    /**
     * Remove a glob star specifier returning the longest parent path
     *
     * <pre>
     * /some/data/&#42;&#42;/file.txt  ->   /some/data
     * /some/data/file.txt     ->   /some/data/file.txt
     * /some&#42;&#42;                 ->   *
     * </pre>
     *
     * @param path
     * @return
     */
    protected removeGlobStar(String path) {

        def p = path?.indexOf('**')
        if( p == -1 )
            return path

        def slash = -1
        for( int i=p-1; i>=0; i-- ) {
            if( path.charAt(i) == '/' as char) {
                slash = i
                break
            }
        }

        if( slash == -1 )
            return '*'

        path.substring(0,slash)
    }

    @Override
    String getBeforeStartScript() {
        if( getPathScheme(targetDir) == 's3' ) {
            final script = getS3UploaderScript()
            if( !script ) throw new IllegalStateException("Missing required nf-amazon module")
            return script.leftTrim()
        }
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile(Path file) {
        "touch ${Escape.path(file)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr(Path file) {
        Escape.path(file)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile(String name, Path target) {
        "cp ${Escape.path(name)} ${Escape.path(target)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String exitFile(Path file) {
        "> ${Escape.path(file)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile(Path file) {
        " < ${Escape.path(file)}"
    }

    @Override
    String getEnvScript(Map environment, boolean container) {
        if( !environment )
            return null

        // create the *bash* environment script
        if( !container ) {
            return TaskProcessor.bashEnvironmentScript(environment)
        }
        else {
            final wrapper = new StringBuilder()
            wrapper << "nxf_container_env() {\n"
            wrapper << 'cat << EOF\n'
            wrapper << TaskProcessor.bashEnvironmentScript(environment, true)
            wrapper << 'EOF\n'
            wrapper << '}\n'
            return wrapper.toString()
        }
    }

}
