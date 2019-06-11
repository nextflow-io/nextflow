/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import java.nio.file.Path

import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.util.Escape

/**
 * Defines the script operation to handle files when running via wr.
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * based on AwsBatchFileCopyStrategy by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WrFileCopyStrategy extends SimpleFileCopyStrategy {

    private final String outputMountLocation
    static final String inputMountLocation = ".inputs"
    private final String workDirScheme
    private final Map<String,Path> inputFiles

    WrFileCopyStrategy() { }

    WrFileCopyStrategy(TaskBean task) {
        super(task)
        this.inputFiles = task.inputFiles
        this.outputMountLocation = ".mnt" + task.workDir?.toString()
        if (this.outputMountLocation.endsWith("/")) {
            this.outputMountLocation = this.outputMountLocation.substring(0, this.outputMountLocation.length() - 1);
        }
        this.workDirScheme = getPathScheme(task.workDir)
    }

    /**
     * {@inheritDoc}
     */
    String getBeforeStartScript() {
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getEnvScript(Map environment, boolean container) {
        if(container) throw new UnsupportedOperationException()

        final result = new StringBuilder()
        final copy = environment ? new HashMap<String,String>(environment) : Collections.<String,String>emptyMap()
        final path = copy.containsKey('PATH')
        // remove any external PATH
        if( path )
            copy.remove('PATH')
        // when a remote bin directory is provided manage it properly
        // *** somehow accept a remoteBinDir opt, and assume it will be mounted
        // somewhere, and add that relative path to PATH?...
        // if( opts.remoteBinDir ) {
        //     result << "${opts.getAwsCli()} s3 cp --recursive --only-show-errors s3:/${opts.remoteBinDir} \$PWD/nextflow-bin\n"
        //     result << "chmod +x \$PWD/nextflow-bin/*\n"
        //     result << "export PATH=\$PWD/nextflow-bin:\$PATH\n"
        // }
        // finally render the environment
        final envSnippet = super.getEnvScript(copy,false)
        if( envSnippet )
            result << envSnippet
        return result.toString()
    }

    /**
     * getWorkDirScheme returns the path scheme (eg. "file" or "s3") for the
     * task's workDir.
    */
    String getWorkDirScheme() {
        return workDirScheme
    }

    /**
     * wrWorkPath returns the same paths that SimpleFileCopyStrategy would end
     * up using, unless path is in S3, in which case it is the file path in the
     * output mount location.
    */
    String wrWorkPath( Path path ) {
        return wrPath(path, outputMountLocation)
    }

    /**
     * wrInputPath returns the same paths that SimpleFileCopyStrategy would end
     * up using, unless path is in S3, in which case it is the file path in the
     * input mount location.
    */
    String wrInputPath( Path path ) {
        return wrPath(path, inputMountLocation)
    }

    /**
     * wrPath returns the same paths that SimpleFileCopyStrategy would end
     * up using, unless path is in S3, in which case it is the file path in the
     * given mount location.
    */
    private String wrPath( Path path, String mountLocation ) {
        if( getPathScheme(path) == 's3' ) {
            if (path == workDir) {
                return "$mountLocation/"
            }
            if (mountLocation == outputMountLocation) {
                String baseName = path.getFileName()
                if (baseName.charAt(0) == '/') {
                    baseName = baseName.substring(1, baseName.length());
                }
                baseName = Escape.path(baseName)
                return "$mountLocation/$baseName"
            }
            String p = path
            if (p.charAt(0) == '/') {
                p = p.substring(1, p.length());
            }
            p = Escape.path(p)
            return "$mountLocation/$p"
        }
        return path
    }

    /**
     * getInputMountPath tells you where you should mount the S3 locations
     * returned by inputBuckets(). Each different bucket should be mounted on
     * its own subdirectoy of this path, named after the bucket.
    */
    String getInputMountPath() {
        return inputMountLocation
    }

    /**
     * getOutputMountPath tells you where you should mount the S3 location
     * returned by outputBucket().
    */
    String getOutputMountPath() {
        return outputMountLocation
    }

    /**
     * inputBuckets returns all the different S3 buckets that would need to be
     * mounted at getInputMountPath() for access to the input files.
    */
    List<String> inputBuckets() {
        List<String> result = []
        inputFiles.each{
            if (getPathScheme(it.value) == 's3') {
                String dir = it.value.getRoot().toString().substring(1)
                dir = dir.substring(0, dir.length() - 1)
                if (!result.contains(dir)) {
                    result << dir
                }
            }
        }
        return result
    }

    /**
     * outputBucket returns the S3 bucket directory that would need to be
     * mounted at getOutputMountPath() in order to output files.
    */
    String outputBucket() {
        if (getPathScheme(workDir) == "s3") {
            return workDir.toString().substring(1)
        }
        return ""
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        def cmd = ''
        def p = targetName.lastIndexOf('/')
        if( p>0 ) {
            cmd += "mkdir -p ${Escape.path(targetName.substring(0,p))} && "
        }
        cmd += stageInCommand(wrInputPath(path.toAbsolutePath()), targetName, stageinMode)

        return cmd
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String stageOutCommand( String source, Path targetDir, String mode ) {
        return stageOutCommand(source, wrWorkPath(targetDir), mode)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile( Path path ) {
        "touch ${wrWorkPath(path)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr( Path path ) {
        wrWorkPath(path)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile( String name, Path target ) {
        "cp ${Escape.path(name)} ${wrWorkPath(target)}"
    }

    /**
     * {@inheritDoc}
     */
    String exitFile( Path path ) {
        "> ${wrWorkPath(path)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile( Path path ) {
        " < ${wrInputPath(path)}"
    }
}