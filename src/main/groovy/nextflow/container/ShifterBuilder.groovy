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
