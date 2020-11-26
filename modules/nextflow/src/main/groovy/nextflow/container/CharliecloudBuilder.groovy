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

import java.nio.file.Path

import nextflow.util.PathTrie
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
        assert image

        result << 'ch-run --no-home --unset-env="*" '

        appendEnv(result)

        if( runOptions )
            result << runOptions.join(' ') << ' '

        makeVolumes(mounts, result)
        
        result << '-c "$PWD" '

        result << image
        result << ' --'

        runCommand = result.toString()

        return this
    }

    protected String composeVolumePath(String path, boolean readOnly = false) {
        return "-b ${escape(path)}:${escape(path)}"
    }

    /**
    * This method override is needed because charliecloud currently can't bind mount directories that do not exist in the container
    * The workaround is to use mkdir to create the directories within the container before they are bind-mounted
    * Can probably be removed once Charliecloud issue https://github.com/hpc/charliecloud/issues/96 has been resolved
    */
    @Override
    protected CharSequence makeVolumes(List<Path> mountPaths, StringBuilder result) {

        def prependDirs = ''

        // add the work-dir to the list of container mounts
        final workDirStr = workDir?.toString()
        final allMounts = new ArrayList<Path>(mountPaths)
        if( workDir )
            allMounts << workDir

        // find the longest commons paths and mount only them
        final trie = new PathTrie()
        for( String it : allMounts ) { trie.add(it) }

        // when mounts are read-only make sure to remove the work-dir path
        final paths = trie.longest()
        if( readOnlyInputs && workDirStr && paths.contains(workDirStr) )
            paths.remove(workDirStr)

        for( String it : paths ) {
            if(!it) continue
            prependDirs += it + ' '
            result << composeVolumePath(it,readOnlyInputs)
            result << ' '
        }

        // when mounts are read-only, make sure to include the work-dir as writable
        if( readOnlyInputs && workDir ) {
            prependDirs += workDirStr + ' '
            result << composeVolumePath(workDirStr)
            result << ' '
        }

        // -- append by default the current path -- this is needed when `scratch` is set to true
        if( mountWorkDir ) {
            prependDirs += '"$PWD"'
            result << composeVolumePath('$PWD')
            result << ' '
        }

        if( prependDirs ) {
            prependDirs = 'ch-run --no-home -w ' + image + ' -- bash -c "mkdir -p ' + prependDirs + '";'
            result.insert(0, prependDirs)
        }
        return result
    }

    @Override
    protected CharSequence makeEnv( env, StringBuilder result = new StringBuilder() ) {

        if( env instanceof Map ) {
            short index = 0
            for( Map.Entry entry : env.entrySet() ) {
                if( index++ ) result << ' '
                // use bash process substitution because --set-env expects a file handle
                result << ("--set-env=<( echo \"${entry.key}=\"${entry.value}\"\" )")
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << '--set-env=<( echo "' << env << '" )'
        }
        else if( env instanceof String ) {
            result << "\${$env:+--set-env=<( echo \"$env=\"\$$env\"\" )}"
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid environment value: $env [${env.class.name}]")
        }

        return result
    }
}
