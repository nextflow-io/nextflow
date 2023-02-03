/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cli

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.k8s.K8sDriverLauncher

/**
 * CLI `kuberun` sub-command
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class KubeRunImpl extends RunImpl {

    static private String POD_NAME = /[a-z0-9]([-a-z0-9]*[a-z0-9])?(\.[a-z0-9]([-a-z0-9]*[a-z0-9])?)*/

    interface Options extends RunImpl.Options {
        int getHeadCpus()
        String getHeadImage()
        String getHeadMemory()
        String getHeadPreScript()
        String getNamespace()
        String getPodImage()
        List<String> getRemoteConfig()
        String getRemoteProfile()
        List<String> getVolumeMounts()

        void setHeadImage(String headImage)
    }

    @Delegate
    Options options

    KubeRunImpl(Options options) {
        super(options)
        this.options = options
    }

    @Override
    protected void checkRunName() {
        if( runName && !runName.matches(POD_NAME) )
            throw new AbortOperationException("Not a valid K8s pod name -- It can only contain lower case alphanumeric characters, '-' or '.', and must start and end with an alphanumeric character")
        super.checkRunName()
        runName = runName.replace('_','-')
    }

    protected boolean background() { launcherOptions.background }

    protected hasAnsiLogFlag() { launcherOptions.hasAnsiLogFlag() }

    @Override
    void run() {
        final scriptArgs = (args?.size()>1 ? args[1..-1] : []) as List<String>
        final pipeline = stdin ? '-' : ( args ? args[0] : null )
        if( !pipeline )
            throw new AbortOperationException("No project name was specified")
        if( hasAnsiLogFlag() )
            log.warn "Ansi logging not supported by kuberun command"
        if( podImage ) {
            log.warn "-pod-image is deprecated (use -head-image instead)"
            headImage = podImage
        }
        checkRunName()
        final driver = new K8sDriverLauncher(cmd: this, runName: runName, headImage: headImage, background: background(), headCpus: headCpus, headMemory: headMemory, headPreScript: headPreScript) 
        driver.run(pipeline, scriptArgs)
        final status = driver.shutdown()
        System.exit(status)
    }

}
