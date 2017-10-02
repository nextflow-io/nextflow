/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileDownloader
import nextflow.processor.TaskBean
import nextflow.util.Escape
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
     * All the input files as a map of {@code <stage name, store path>} pairs
     */
    Map<String,Path> inputFiles

    /**
     * The output file names
     */
    List<String> outputFiles

    /**
     * Path where output files need to be copied
     */
    Path targetDir

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

    private FileDownloader downloader

    SimpleFileCopyStrategy() {
        this.inputFiles = [:]
        this.outputFiles = []
        this.downloader = new FileDownloader()
    }


    SimpleFileCopyStrategy( TaskBean bean ) {
        this.inputFiles = bean.inputFiles
        this.outputFiles = bean.outputFiles
        this.targetDir = bean.targetDir
        this.stageinMode = bean.stageInMode
        this.stageoutMode = bean.stageOutMode
        this.downloader = new FileDownloader()
    }


    /**
     * Given a map of the input file parameters with respective values,
     * create the BASH script to stage them into the task working space
     *
     * @param inputs An associative array mapping each {@code FileInParam} to the corresponding file (or generic value)
     * @return The BASH script to stage them
     */
    String getStageInputFilesScript() {
        assert inputFiles != null

        final staging = downloader.downloadForeignFiles(inputFiles)

        def delete = []
        def links = []
        staging.each { stageName, storePath ->

            // delete all previous files with the same name
            delete << "rm -f ${Escape.path(stageName)}"

            // link them
            links << stageInputFile( storePath, stageName )

        }
        links << '' // just to have new-line at the end of the script

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
    String getUnstageOutputFilesScript() {

        // collect all the expected names (pattern) for files to be un-staged
        def result = []
        def normalized = normalizeGlobStarPaths(outputFiles)

        // create a bash script that will copy the out file to the working directory
        log.trace "Unstaging file path: $normalized"
        if( normalized ) {
            result << ""
            result << "mkdir -p ${Escape.path(targetDir)}"
            for( int i=0; i<normalized.size(); i++ ) {
                final path = normalized[i]
                final cmd = stageOutCommand(path, targetDir.toString(), stageoutMode) + ' || true' // <-- add true to avoid it stops on errors
                result << cmd
            }
        }

        return result.join(separatorChar)
    }

    protected String stageInCommand( String source, String target, String mode ) {

        if( mode == 'symlink' || !mode )
            return "ln -s ${Escape.path(source)} ${Escape.path(target)}"

        if( mode == 'link' )
            return "ln ${Escape.path(source)} ${Escape.path(target)}"

        if( mode == 'copy' || !mode )
            return "cp -fRL ${Escape.path(source)} ${Escape.path(target)}"

        throw new IllegalArgumentException("Unknown stage-in strategy: $mode")
    }

    /**
     * Compose a `cp` or `mv` command given the source and target paths
     *
     * @param source
     * @param target
     * @param move When {@code true} use unix {@code mv} instead of {@code cp}
     * @return A shell copy or move command string
     */
    protected String stageOutCommand( String source, String target, String mode = null) {

        def cmd
        if( mode == 'copy' || !mode )
            cmd = 'cp -fRL'
        else if( mode == 'move' )
            cmd = 'mv -f'
        else if( mode == 'rsync' )
            return "rsync -rRl ${Escape.path(source)} ${Escape.path(target)}"
        else
            throw new IllegalArgumentException("Unknown stage-out strategy: $mode")

        final p = source.lastIndexOf('/')
        if( p<=0 ) {
            return "$cmd ${Escape.path(source)} ${Escape.path(target)}"
        }

        def path  = new File(target,source.substring(0,p)).toString()
        return "mkdir -p ${Escape.path(path)} && $cmd ${Escape.path(source)} ${Escape.path(path)}"
    }

    /**
     * Normalize path that contains glob star wildcards (i.e. double star character) since
     * they are not supported by plain BASH
     *
     * @param files
     * @return
     */
    protected List<String> normalizeGlobStarPaths( List<String> files ) {

        def result = []
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
        Escape.path(file)
    }

}
