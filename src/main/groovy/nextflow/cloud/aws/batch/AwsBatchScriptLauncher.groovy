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

package nextflow.cloud.aws.batch

import groovy.transform.CompileStatic
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 * Implements BASH launcher script for AWS Batch jobs
 */
@CompileStatic
class AwsBatchScriptLauncher extends BashWrapperBuilder {

    AwsBatchScriptLauncher(TaskBean bean, AwsOptions opts ) {
        super(bean, new AwsBatchFileCopyStrategy(bean,opts))
        // enable the copying of output file to the S3 work dir
        if( !scratch )
            scratch = true
        // include task script as an input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_SCRIPT] = bean.workDir.resolve(TaskRun.CMD_SCRIPT)
        // add the wrapper file when stats are enabled
        if( bean.statsEnabled ) {
            bean.inputFiles[TaskRun.CMD_STUB] = bean.workDir.resolve(TaskRun.CMD_STUB)
        }
        // include task stdin file
        if( bean.input != null ) {
            bean.inputFiles[TaskRun.CMD_INFILE] = bean.workDir.resolve(TaskRun.CMD_INFILE)
        }
    }

}