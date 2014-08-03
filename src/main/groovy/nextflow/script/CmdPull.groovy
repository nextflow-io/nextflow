/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.script
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException

/**
 * CLI sub-command UPDATE
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Download or updated a pipeline in the local repository")
class CmdPull implements CmdX {

    @Parameter(required=true, description = 'The name of the pipeline to pull', arity = 1)
    List<String> args

    @Parameter(names='-all', description = 'Update all installed pipelines', arity = 0)
    boolean all

    @Override
    final String getName() { 'pull' }

    @Override
    void run() {

        if( !all && !args )
            throw new AbortOperationException('Missing argument')

        (all ? new PipelineManager().list() : args.toList()) .each {
            log.info "Pulling $it ..."
            new PipelineManager(it).download()
        }

    }

}
