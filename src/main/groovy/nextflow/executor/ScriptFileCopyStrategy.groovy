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
     * Resolve all foreign file to storing then in the pipeline working directory
     *
     * @param inputFiles A map of <file name, path location> of files. It can hold both local and foreign files
     * @return
     *      A map <file name, path location> of files in which all foreign files have been replaced
     *      by local copies in the working directory
     */
    Map<String,Path> resolveForeignFiles(Map<String,Path> inputFiles)

    /**
     * @param inputFiles All the input files as a map of {@code <stage name, store path>} pairs
     * @return A BASH snippet included in the wrapper script that stages the task input files
     */
    String getStageInputFilesScript(Map<String,Path> inputFiles)

    /**
     * @param outputFiles List of the output file names to unstage
     * @param targetDir The directory where output files need to be unstaged ie. stored
     * @return A BASH snippet included in the wrapper script that un-stages the task output files
     */
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir)

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

    /**
     * BASH command to redirect the specified file on the standard input
     *
     * @param file The file to be piped as an input
     * @return BASH redirection for the specified file
     */
    String pipeInputFile( Path file )

    /**
     * @param environment The task environment
     * @return The environment initialisation snippet
     */
    String getEnvScript(Map environment, String wrapName)

}
