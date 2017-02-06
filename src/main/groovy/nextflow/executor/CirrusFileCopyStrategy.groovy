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
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 * Defines the script operation to handle file when running in the Cirrus cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CirrusFileCopyStrategy extends SimpleFileCopyStrategy {

    Path workDir

    CirrusFileCopyStrategy( TaskBean task ) {
        this.inputFiles = task.getInputFiles()
        this.outputFiles = task.getOutputFiles()
        this.targetDir = task.getTargetDir()
        this.workDir = task.workDir
    }

    /**
     * @return A script snipped that download from S3 the task scripts:
     * {@code .command.env}, {@code .command.sh}, {@code .command.in},
     * etc.
     */
    String getBeforeStartScript() {

        def env = workDir.resolve(TaskRun.CMD_ENV)
        def script = workDir.resolve(TaskRun.CMD_SCRIPT)
        def infile = workDir.resolve(TaskRun.CMD_INFILE)
        def wrapper = workDir.resolve(TaskRun.CMD_STUB)

        def ops = []
        ops << '# fetch scripts'
        ops << "es3 test s3:/${env} && es3 cat s3:/${env} > ${env.name}"
        ops << "es3 test s3:/${script} && es3 cat s3:/${script} > ${script.name}"
        ops << "es3 test s3:/${infile} && es3 cat s3:/${infile} > ${infile.name}"
        ops << "es3 test s3:/${wrapper} && es3 cat s3:/${wrapper} > ${wrapper.name}"
        ops << ''

        ops.join('\n')
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        def cmd = "es3 -q -v0 --no-stats sync s3:/${path} ."
        if( path.name != targetName ) {
            def p = targetName.lastIndexOf('/')
            if( p>0 ) {
                cmd += " && mkdir -p '${targetName.substring(0,p)}'"
            }
            cmd += " && mv ${path.name} ${targetName}"
        }

        return cmd
    }

    /**
     * {@inheritDoc}
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
            normalized.each {
                result << "es3 -q -v0 --no-stats sync $it s3:/${targetDir} || true" // <-- add true to avoid it stops on errors
            }
        }

        return result.join(separatorChar)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile( Path file ) {
        "es3 touch s3:/${file}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr( Path file ) {
        file.getFileName().toString()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile( String name, Path target ) {
        "es3 -q -v 0 --no-stats sync ${name} s3:/${target}"
    }

    /**
     * {@inheritDoc}
     */
    String exitFile( Path file ) {
        "${file.name} && es3 -q -v 0 --no-stats sync ${file.name} s3:/${file.parent} || true"
    }
}
