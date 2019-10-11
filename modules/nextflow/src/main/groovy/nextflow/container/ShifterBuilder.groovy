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

package nextflow.container
/**
 * Wrap a task execution in a Shifter container
 *
 * See https://github.com/NERSC/shifter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ShifterBuilder extends ContainerBuilder<ShifterBuilder> {

    private boolean verbose

    ShifterBuilder( String image ) {
        assert image
        this.image = image
    }

    @Override
    ShifterBuilder build(StringBuilder result) {
        assert image

        appendEnv(result)

        result << 'shifter '

        if( verbose )
            result << '--verbose '

        result << '--image ' << image

        runCommand = result.toString()
        return this
    }

    ShifterBuilder params( Map params ) {

        if( params.containsKey('verbose') )
            this.verbose = params.verbose.toString() == 'true'

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        return this
    }

    @Override
    String getRunCommand() {
        def run = super.getRunCommand()
        def result = """\
        shifterimg pull $image
        shifterimg lookup $image
        while ! shifterimg lookup $image; do
            sleep 5
            STATUS=\$(shifterimg -v pull $image | tail -n2 | head -n1 | awk \'{print \$6}\')
            [[ \$STATUS == "FAILURE" || -z \$STATUS ]] && echo "Shifter failed to pull image \'$image\'" >&2  && exit 1
        done
        """.stripIndent()
        result += run
        return result
    }

    /**
     * Get the env command line option for the give environment definition
     *
     * @param env
     * @param result
     * @return
     */
    protected CharSequence makeEnv( env, StringBuilder result = new StringBuilder() ) {

        if( env instanceof Map ) {
            short index = 0
            for( Map.Entry entry : env.entrySet() ) {
                if( index++ ) result << ' '
                result << ("${entry.key}=\"${entry.value}\"")
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << env
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid environment value: $env [${env.class.name}]")
        }

        return result
    }
}
