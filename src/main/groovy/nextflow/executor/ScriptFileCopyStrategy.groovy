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
/**
 * Declares file operations used in the task wrapper script
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface ScriptFileCopyStrategy {

    /**
     * @return A BASH snippet included in the wrapper script to be executed *before* the task execution itself
     */
    String getBeforeStartScript()

    /**
     * @return A BASH snippet included in the wrapper script that stages the task input files
     */
    String getStageInputFilesScript()

    /**
     * @return A BASH snippet included in the wrapper script that un-stages the task output files
     */
    String getUnstageOutputFilesScript()

    /**
     * Command to 'touch' a file
     *
     * @param file The absolute path of the file to be 'touched'
     * @return The touch command script
     */
    String touchFile( Path file )

    /**
     * @param file A file {@link Path}
     * @return the path string given a {@link Path} object
     */
    String fileStr( Path file )

    /**
     * The command to copy a path to the specified target folder
     *
     * @param name The name of the file to copy
     * @param target The target directory
     * @return The command to copy the file in the target folder
     */
    String copyFile( String name, Path target )

    /**
     * @param file The {@link Path} to the exit file
     * @return The exit file as a string
     */
    String exitFile( Path file )

    String pipeInputFile( Path file )

}
