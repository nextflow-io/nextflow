/*
 * Copyright 2020, Seqera Labs
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

package nextflow.container
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Implements a builder for Charliecloud containerisation
 *
 * see https://hpc.github.io/charliecloud/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Patrick Hüther <patrick.huether@gmail.com>
 */
@CompileStatic
@Slf4j
class CharliecloudBuilder extends ContainerBuilder<CharliecloudBuilder> {

    CharliecloudBuilder(String name) {
        this.image = name
    }

    @Override
    CharliecloudBuilder params(Map params) {

        if( params.containsKey('temp') )
            this.temp = params.temp

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        if( params.containsKey('runOptions') )
            addRunOptions(params.runOptions.toString())

        if( params.containsKey('readOnlyInputs') )
            this.readOnlyInputs = params.readOnlyInputs?.toString() == 'true'

        return this
    }

    CharliecloudBuilder addRunOptions(String str) {
        runOptions.add(str)
        return this
    }

    @Override
    CharliecloudBuilder build(StringBuilder result) {

        // charliecloud currently can't bind mount directories that do not exist in the container
        // charliecloud issue https://github.com/hpc/charliecloud/issues/96
        // TODO create override for makeVolumes, this needs to be done for all mounts
        result << 'ch-run --no-home -w ' + image + ' -- bash -c "mkdir -p $PWD";'

        result << 'ch-run --no-home '
        result << '--unset-env="*" '

        // get environment from container
        // this is needed to workaround the fact that charliecloud ignores ENV layers of docker images
        // charliecloud issue https://github.com/hpc/charliecloud/issues/719
        // TODO create override for appendEnv, that also sets vars defined in env scope
        result << '--set-env=' + image + '/etc/environment ' 

        if( runOptions )
            result << runOptions.join(' ') << ' '

        // mount the input folders
        result << makeVolumes(mounts)
        result << '-c "$PWD" '

        result << image
        result << ' --'

        runCommand = result.toString()

        return this
    }

    protected String composeVolumePath(String path, boolean readOnly = false) {
        return "-b ${escape(path)}:${escape(path)}"
    }
}
