/*
 * Copyright 2013-2026, Seqera Labs
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

    SarusBuilder(String image, SarusConfig config) {
        this.image = image

        if( config.runOptions )
            addRunOptions(config.runOptions)
        this.tty = config.tty
        this.verbose = config.verbose
    }

    SarusBuilder(String image) {
        this(image, new SarusConfig([:]))
    }

    SarusBuilder params( Map params ) {
        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        return this
    }

    @Override
    SarusBuilder build(StringBuilder result) {
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
        result << '-w "$NXF_TASK_WORKDIR" '

        if( runOptions )
            result << runOptions.join(' ') << ' '

        // finally the container name
        result << image

        runCommand = result.toString()
        return this
    }

    @Override
    String getRunCommand() {
        def run = super.getRunCommand()
        def result = """\
        sarus pull $image 1>&2
        """.stripIndent(true)
        result += run
        return result
    }

    @Override
    protected String composeVolumePath( String path, boolean readOnly = false ) {
        return "--mount=type=bind,source=${escape(path)},destination=${escape(path)}"
    }
}
