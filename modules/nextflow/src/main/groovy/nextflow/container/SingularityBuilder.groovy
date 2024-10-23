/*
 * Copyright 2013-2024, Seqera Labs
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

    private boolean homeMount

    private boolean newPidNamespace

    private String runCmd0

    private Boolean ociMode

    SingularityBuilder(String name) {
        this.image = name
        this.homeMount = defaultHomeMount()
        this.autoMounts = defaultAutoMounts()
        this.newPidNamespace = defaultNewPidNamespace()
        this.runCmd0 = defaultRunCommand()
    }

    private boolean defaultHomeMount() {
        SysEnv.get("NXF_${getBinaryName().toUpperCase()}_HOME_MOUNT", 'false').toString() == 'true'
    }

    private boolean defaultNewPidNamespace() {
        SysEnv.get("NXF_${getBinaryName().toUpperCase()}_NEW_PID_NAMESPACE", 'true').toString() == 'true'
    }

    private boolean defaultAutoMounts() {
        SysEnv.get("NXF_${getBinaryName().toUpperCase()}_AUTO_MOUNTS", 'true').toString() == 'true'
    }

    private String defaultRunCommand() {
        final result = SysEnv.get("NXF_${getBinaryName().toUpperCase()}_RUN_COMMAND", 'exec')
        if( result !in ['run','exec'] )
            throw new IllegalArgumentException("Invalid singularity launch command '$result' - it should be either 'run' or 'exec'")
        return result
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

        if( params.autoMounts!=null )
            autoMounts = params.autoMounts.toString() == 'true'

        if( params.newPidNamespace!=null )
            newPidNamespace = params.newPidNamespace.toString() == 'true'

        if( params.containsKey('readOnlyInputs') )
            this.readOnlyInputs = params.readOnlyInputs?.toString() == 'true'

        // note: 'oci' flag should be ignored by Apptainer sub-class
        if( params.oci!=null && this.class==SingularityBuilder )
            ociMode = params.oci.toString() == 'true'
        else if( params.ociMode!=null && this.class==SingularityBuilder )
            ociMode = params.ociMode.toString() == 'true'

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

        result << runCmd0 << ' '

        if( !homeMount )
            result << '--no-home '

        if( newPidNamespace && !ociMode )
            result << '--pid '

        if( ociMode != null )
            result << (ociMode ? '--oci ' : '--no-oci ')

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
        // add magic variables required by singularity to run in OCI-mode
        if( ociMode ) {
            result .append('${XDG_RUNTIME_DIR:+XDG_RUNTIME_DIR="$XDG_RUNTIME_DIR"} ')
            result .append('${DBUS_SESSION_BUS_ADDRESS:+DBUS_SESSION_BUS_ADDRESS="$DBUS_SESSION_BUS_ADDRESS"} ')
        }
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
            result += entryPoint ? " $entryPoint -c \"cd \$NXF_TASK_WORKDIR; $launcher\"" : " $launcher"
            return result
        }
        return getRunCommand() + ' ' + launcher
    }

}
