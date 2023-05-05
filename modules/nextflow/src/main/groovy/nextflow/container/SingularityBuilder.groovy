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
import nextflow.SysEnv

/**
 * Implements a builder for Singularity containerisation
 *
 * see http://singularity.lbl.gov
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Slf4j
class SingularityBuilder extends ContainerBuilder<SingularityBuilder> {

    private boolean autoMounts

    private boolean newPidNamespace

    SingularityBuilder(String name) {
        this.image = name
        this.newPidNamespace = defaultNewPidNamespace()
    }

    private boolean defaultNewPidNamespace() {
        SysEnv.get("NXF_${getBinaryName().toUpperCase()}_NEW_PID_NAMESPACE", 'true').toString() == 'true'
    }

    protected String getBinaryName() { 'singularity' }

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

        if( params.newPidNamespace!=null )
            newPidNamespace = params.newPidNamespace.toString() == 'true'

        if( params.containsKey('readOnlyInputs') )
            this.readOnlyInputs = params.readOnlyInputs?.toString() == 'true'

        return this
    }

    @Override
    SingularityBuilder addRunOptions(String str) {
        super.addRunOptions(str)
    }

    @Override
    SingularityBuilder build(StringBuilder result) {

        result << 'set +u; env - PATH="$PATH" '

        appendEnv(result)

        result << getBinaryName() << ' '

        if( engineOptions )
            result << engineOptions.join(' ') << ' '

        result << 'exec '

        if( newPidNamespace )
            result << '--pid '

        if( autoMounts ) {
            makeVolumes(mounts, result)
        }

        if( runOptions )
            result << runOptions.join(' ') << ' '

        result << image

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
    protected CharSequence appendEnv(StringBuilder result) {
        makeEnv('TMP',result) .append(' ')
        makeEnv('TMPDIR',result) .append(' ')
        super.appendEnv(result)
    }

    protected String prefixEnv(String key) {
        final PREFIX = getBinaryName().toUpperCase()
        if( key.startsWith(PREFIX+'_') )
            return key
        if( key.startsWith(PREFIX+'ENV_') )
            return key
        return PREFIX+'ENV_'+key
    }

    protected String quoteValue(String env) {
        if( !env )
            return env
        final p=env.indexOf('=')
        return p==-1 ? quoteValue0(env) : env.substring(0,p) + '=' + quoteValue0(env.substring(p+1))
    }

    private String quoteValue0(String value) {
        if( !value )
            return value
        if( value.startsWith('"') && value.endsWith('"') )
            return value
        if( value.startsWith("'") && value.endsWith("'") )
            return value
        return '"'  + value + '"'
    }

    @Override
    protected StringBuilder makeEnv( env, StringBuilder result = new StringBuilder() ) {

        if( env instanceof Map ) {
            int index=0
            for( Map.Entry entry : env.entrySet() ) {
                if( index++ ) result << ' '
                result << "${prefixEnv(entry.key.toString())}=\"${entry.value}\""
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << prefixEnv(quoteValue(env))
        }
        else if( env instanceof String ) {
            result << "\${$env:+${prefixEnv(env)}=\"\$$env\"}"
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid environment value: $env [${env.getClass().name}]")
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
    String getRunCommand(String launcher) {
        if( !runCommand )
            throw new IllegalStateException("Missing `runCommand` -- make sure `build` method has been invoked")

        if( launcher ) {
            def result = getRunCommand()
            result += entryPoint ? " $entryPoint -c \"cd \$PWD; $launcher\"" : " $launcher"
            return result
        }
        return getRunCommand() + ' ' + launcher
    }

}
