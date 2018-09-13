/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.ga4gh.tes.executor

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.executor.ScriptFileCopyStrategy
import nextflow.executor.SimpleFileCopyStrategy
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
    String getEnvScript(Map env, String wrapName=null) {
        return SimpleFileCopyStrategy.getEnvScript0(env,wrapName)
    }
}