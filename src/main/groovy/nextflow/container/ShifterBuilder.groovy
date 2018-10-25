/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

    static private String SHIFTER_HELPERS = '''
        function shifter_img() {
          local cmd=$1
          local image=$2
          shifterimg -v $cmd $image |  awk -F: '$0~/"status":/{gsub("[\\", ]","",$2);print $2}\'
        }

        function shifter_pull() {
          local image=$1
          local STATUS=$(shifter_img lookup $image)
          if [[ $STATUS != READY && $STATUS != '' ]]; then
            STATUS=$(shifter_img pull $image)
            while [[ $STATUS != READY && $STATUS != FAILURE && $STATUS != '' ]]; do
              sleep 5
              STATUS=$(shifter_img pull $image)
            done
          fi

          [[ $STATUS == FAILURE || $STATUS == '' ]] && echo "Shifter failed to pull image \\`$image\\`" >&2  && exit 1
        }
        '''.stripIndent()


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

    StringBuilder appendHelpers( StringBuilder wrapper ) {
        wrapper << SHIFTER_HELPERS << '\n'
    }

    @Override
    String getRunCommand() {
        def run = super.getRunCommand()
        def result = "shifter_pull $image\n"
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
