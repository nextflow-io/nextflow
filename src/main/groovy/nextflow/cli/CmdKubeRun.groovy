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

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.k8s.K8sDriverLauncher
/**
 * Extends `run` command to support Kubernetes deployment
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Execute a workflow in a Kubernetes cluster (experimental)")
class CmdKubeRun extends CmdRun {

    static private String POD_NAME = /[a-z0-9]([-a-z0-9]*[a-z0-9])?(\.[a-z0-9]([-a-z0-9]*[a-z0-9])?)*/

    /**
     * One or more volume claims mounts
     */
    @Parameter(names = ['-v','-volume-mount'], description = 'Volume claim mounts eg. my-pvc:/mnt/path')
    List<String> volMounts

    @Parameter(names = ['-n','-namespace'], description = 'Specify the K8s namespace to use')
    String namespace

    @Override
    String getName() { 'kuberun' }

    @Override
    protected void checkRunName() {
        if( runName && !runName.matches(POD_NAME) )
            throw new AbortOperationException("Not a valid K8s pod name -- It can only contain lower case alphanumeric characters, '-' or '.', and must start and end with an alphanumeric character")
        super.checkRunName()
        runName = runName.replace('_','-')
    }

    @Override
    void run() {
        final scriptArgs = (args?.size()>1 ? args[1..-1] : []) as List<String>
        final pipeline = stdin ? '-' : ( args ? args[0] : null )
        if( !pipeline )
            throw new AbortOperationException("No project name was specified")

        checkRunName()
        new K8sDriverLauncher(cmd: this, runName: runName).run(pipeline, scriptArgs)
    }

}
