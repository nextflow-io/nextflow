/*
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

package nextflow.ga4gh.tes.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.executor.ScriptFileCopyStrategy
import nextflow.processor.TaskProcessor
import nextflow.util.Escape
/**
 * Implements file inputs/outputs staging strategy for tasks executed
 * by TES executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TesFileCopyStrategy implements ScriptFileCopyStrategy {

    /**
     * {@inheritDoc}
     */
    @Override
    String getBeforeStartScript() {
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    Map<String, Path> resolveForeignFiles(Map<String, Path> inputFile) {
        return inputFile
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {
        return null
    }

    @Override
    String touchFile(Path file) {
        return ''
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr(Path file) {
        Escape.path(file.getFileName())
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile(String name, Path target) {
        return 'true'
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String exitFile(Path file) {
        "> ${Escape.path(file.getName())}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile(Path file) {
        " < ${Escape.path(file.getName())}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getEnvScript(Map env, boolean container) {
        if(container) throw new UnsupportedOperationException()
        TaskProcessor.bashEnvironmentScript(env,false)
    }
}