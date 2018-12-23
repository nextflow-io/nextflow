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

package nextflow.cloud.aws.batch

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.util.Escape

/**
 * Defines the script operation to handle file when running in the Cirrus cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchFileCopyStrategy extends SimpleFileCopyStrategy {

    private AwsOptions opts

    private Map<String,String> environment

    AwsBatchFileCopyStrategy(TaskBean task, AwsOptions opts ) {
        super(task)
        this.opts = opts
        this.environment = task.environment
    }

    /**
     * @return A script snippet that download from S3 the task scripts:
     * {@code .command.env}, {@code .command.sh}, {@code .command.in},
     * etc.
     */
    String getBeforeStartScript() {
        S3Helper.getUploaderScript(opts)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getEnvScript(Map environment, String handler) {
        if( handler )
            throw new IllegalArgumentException("Parameter `wrapHandler` not supported by ${this.class.simpleName}")

        final result = new StringBuilder()
        final copy = environment ? new HashMap<String,String>(environment) : Collections.<String,String>emptyMap()
        final path = copy.containsKey('PATH')
        // remove any external PATH
        if( path )
            copy.remove('PATH')
        // when a remote bin directory is provide managed it properly
        if( opts.remoteBinDir ) {
            result << "${opts.getAwsCli()} s3 cp --recursive --only-show-errors s3:/${opts.remoteBinDir} \$PWD/nextflow-bin\n"
            result << "chmod +x \$PWD/nextflow-bin/*\n"
            result << "export PATH=\$PWD/nextflow-bin:\$PATH\n"
        }
        // finally render the environment
        final envSnippet = super.getEnvScript(copy,null)
        if( envSnippet )
            result << envSnippet
        return result.toString()
    }

    @Override
    String getStageInputFilesScript(Map<String,Path> inputFiles) {
        def result = 'downloads=()\n'
        result += super.getStageInputFilesScript(inputFiles) + '\n'
        result += 'nxf_parallel "${downloads[@]}"\n'
        return result
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        // third param should not be escaped, because it's used in the grep match rule
        "downloads+=(\"nxf_s3_download s3:/${Escape.path(path)} ${Escape.path(targetName)}\")"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {

        // collect all the expected names (pattern) for files to be un-staged
        def result = []
        def normalized = normalizeGlobStarPaths(outputFiles)

        // create a bash script that will copy the out file to the working directory
        log.trace "[AWS BATCH] Unstaging file path: $normalized"
        if( normalized ) {
            result << 'uploads=()'
            normalized.each {
                result << "uploads+=(\"nxf_s3_upload '${Escape.path(it)}' s3:/${Escape.path(targetDir)}\")" // <-- add true to avoid it stops on errors
            }
            result << 'nxf_parallel "${uploads[@]}"'
        }

        return result.join(separatorChar)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile( Path file ) {
        final aws = opts.getAwsCli()
        "echo start | $aws s3 cp --only-show-errors - s3:/${Escape.path(file)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr( Path path ) {
        Escape.path(path.getFileName())
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile( String name, Path target ) {
        "nxf_s3_upload ${Escape.path(name)} s3:/${Escape.path(target.getParent())}"
    }

    /**
     * {@inheritDoc}
     */
    String exitFile( Path path ) {
        final aws = opts.getAwsCli()
        "| $aws s3 cp --only-show-errors - s3:/${Escape.path(path)} || true"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile( Path path ) {
        " < ${Escape.path(path.getFileName())}"
    }
}