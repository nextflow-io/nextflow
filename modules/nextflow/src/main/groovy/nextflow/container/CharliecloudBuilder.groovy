/*
 * Copyright 2013-2023, Seqera Labs
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
 * @author Laurent Modolo <laurent.modolo@ens-lyon.fr>
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

        result << 'ch-run --unset-env="*" -c "$PWD" --set-env '
        if (!readOnlyInputs)
            result << '-w '

        appendEnv(result)

        if( temp )
            result << "-b $temp:/tmp "

        makeVolumes(mounts, result)

        if( runOptions )
            result << runOptions.join(' ') << ' '

        result << image
        result << ' --'

        runCommand = result.toString()

        return this
    }

    protected String getRoot(String path) {
        def rootPath = path.tokenize("/")

        if (rootPath.size() >= 1 && path[0] == '/')
            rootPath = "/${rootPath[0]}"
        else
            throw new IllegalArgumentException("Not a valid working directory value (absolute path?): ${path}")

        return rootPath
    }
    
    @Override
    protected String composeVolumePath(String path, boolean readOnlyInputs = false) {
        def mountCmd = "-b ${escape(path)}"
        if (readOnlyInputs)
            mountCmd = "-b ${getRoot(escape(path))}"
        return mountCmd
    }

    @Override
    protected StringBuilder makeEnv( env, StringBuilder result = new StringBuilder() ) {

        if( env instanceof Map ) {
            short index = 0
            for( Map.Entry entry : env.entrySet() ) {
                if( index++ ) result << ' '
                result << ("--set-env=${entry.key}=${entry.value}")
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << "--set-env=" << env
        }
        else if( env instanceof String ) {
            result << "\${$env:+--set-env=$env=\$$env}"
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid environment value: $env [${env.class.name}]")
        }

        return result
    }
}
