/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

import groovy.util.logging.Slf4j

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import nextflow.util.Escape
/**
 * Implements a builder for Singularity containerisation
 *
 * see http://singularity.lbl.gov
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Slf4j
class SingularityBuilder extends ContainerBuilder {

    private String entryPoint

    private String runCommand

    private boolean autoMounts

    SingularityBuilder(String name) {
        this.image = name
    }

    @Override
    SingularityBuilder params(Map params) {

        if( params.containsKey('temp') )
            this.temp = params.temp

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        if( params.containsKey('engineOptions') )
            addEngineOptions(params.engineOptions.toString())

        if( params.containsKey('runOptions') )
            addRunOptions(params.runOptions.toString())

        if( params.autoMounts )
            autoMounts = params.autoMounts.toString() == 'true'

        if( params.containsKey('readOnlyInputs') )
            this.readOnlyInputs = params.readOnlyInputs?.toString() == 'true'

        return this
    }

    SingularityBuilder addRunOptions(String str) {
        runOptions.add(str)
        return this
    }

    @Override
    SingularityBuilder build(StringBuilder result) {

        result << 'set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" '

        result << 'singularity '

        if( engineOptions )
            result << engineOptions.join(' ') << ' '

        result << 'exec '

        if( autoMounts ) {
            makeVolumes(mounts, result)
            result << ' '
        }

        if( runOptions )
            result << runOptions.join(' ') << ' '

        result << image

        if( entryPoint )
            result  << ' ' << entryPoint

        runCommand = result.toString()

        return this
    }

    protected String composeVolumePath( String path, boolean readOnly = false ) {
        def result = "-B ${escape(path)}"
        if( readOnly )
            result += ":${escape(path)}:ro"
        return result
    }

    @Override
    protected CharSequence makeEnv( env, StringBuilder result = new StringBuilder() ) {
        // append the environment configuration
        if( env instanceof File ) {
            env = env.toPath()
        }
        if( env instanceof Path ) {
            result << 'export BASH_ENV="' << Escape.path(env) << '"; '
        }
        else if( env instanceof Map ) {
            short index = 0
            for( Map.Entry entry : env.entrySet() ) {
                if( index++ ) result << ' '
                result << "export ${entry.key}=\"${entry.value}\"; "
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << 'export ' << env << '; '
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid environment value: $env [${env.class.name}]")
        }

        return result
    }

    String getEnvExports() {
        def result = new StringBuilder()
        for( def entry : env ) {
            makeEnv(entry, result)
        }
        return result.toString()
    }

    @Override
    String getRunCommand() {
        return runCommand
    }

}
