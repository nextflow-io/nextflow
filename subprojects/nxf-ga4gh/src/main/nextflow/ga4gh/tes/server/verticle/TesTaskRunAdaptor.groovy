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

package nextflow.ga4gh.tes.server.verticle

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.file.FileHelper
import nextflow.file.FileHolder
import nextflow.ga4gh.tes.server.model.TesInput
import nextflow.processor.TaskRun
import nextflow.script.InParam
import nextflow.script.OutParam
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TesTaskRunAdaptor extends TaskRun {

    /*
     * the container image name
     */
    @PackageScope String image

    /**
     * The task input files, each pair represent a < stage-path, remote-path > pair
     */
    @PackageScope Map<String,Path> inputFilesMap = [:]

    void addInputFiles(List<TesInput> inputFiles) {
        assert workDir
        if( !inputFiles ) return

        for( TesInput input : inputFiles ) {
            if( !input.path && !input.url )
                continue
            if( !input.path )
                throw new TaskServiceApiException(400, "Malformed TES request -- Input file missing path for: $input.path")
            if( !input.path.startsWith('/') )
                throw new TaskServiceApiException(400, "Malformed TES request -- Missing input path for file: $input.url")
            if( !input.url )
                throw new TaskServiceApiException(400, "Malformed TES request -- Missing input url for file: $input.path")

            final String stage = "${workDir}${input.path}"
            final Path url = FileHelper.asPath( input.url )
            inputFilesMap.put( stage, url )
        }

    }

    @Override
    String getContainer() {
        image
    }

    @Override
    protected Map<String, String> getInputEnvironment() {
        return Collections.emptyMap()
    }

    @Override
    Map<String, Path> getInputFilesMap() {
        return inputFilesMap
    }

    @Override
    List<String> getOutputFilesNames() {
        return null //super.getOutputFilesNames()
    }

    @Override
    List<String> getStagedInputs() {
        throw new UnsupportedOperationException()
    }

    @Override
    Map<InParam, List<FileHolder>> getInputFiles() {
        throw new UnsupportedOperationException()
    }

    @Override
    <T extends InParam> Map<T, Object> getInputsByType(Class<T>... types) {
        throw new UnsupportedOperationException()
    }


    @Override
    <T extends OutParam> Map<T, Object> getOutputsByType(Class<T>... types) {
        throw new UnsupportedOperationException()
    }

    @Override
    Map<OutParam, Object> getOutputs() {
        throw new UnsupportedOperationException()
    }
}
