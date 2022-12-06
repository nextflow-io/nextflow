/*
 * Copyright 2022, Pawsey Supercomputing Research Centre
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
 * Wrap a task execution in a Sarus container
 *
 * See https://sarus.readthedocs.io
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
class SarusBuilder extends ContainerBuilder<SarusBuilder> {

    private boolean tty

    private boolean verbose

    SarusBuilder( String image ) {
        assert image
        this.image = image
    }

    @Override
    SarusBuilder build(StringBuilder result) {
        assert image

        result << 'sarus '

        if( verbose )
            result << '--verbose '

        result << 'run '

        if( tty )
            result << '-t '

        // add the environment
        appendEnv(result)

        // mount the input folders
        result << makeVolumes(mounts)
        result << '-w "$PWD" '

        if( runOptions )
            result << runOptions.join(' ') << ' '

        // finally the container name
        result << image

        runCommand = result.toString()
        return this
    }

    SarusBuilder params( Map params ) {

        if( params.containsKey('verbose') )
            this.verbose = params.verbose.toString() == 'true'

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        if( params.containsKey('runOptions') )
            addRunOptions(params.runOptions.toString())

        if( params.containsKey('tty') )
            this.tty = params.tty?.toString() == 'true'

        return this
    }

    @Override
    String getRunCommand() {
        def run = super.getRunCommand()
        def result = """\
        sarus pull $image 1>&2
        """.stripIndent()
        result += run
        return result
    }

    @Override
    protected String composeVolumePath( String path, boolean readOnly = false ) {
        return "--mount=type=bind,source=${escape(path)},destination=${escape(path)}"
    }
}
