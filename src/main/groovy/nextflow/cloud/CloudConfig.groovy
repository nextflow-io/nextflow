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

package nextflow.cloud
import static nextflow.Const.ROLE_WORKER

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.config.CascadingConfig
import nextflow.config.ConfigField
import nextflow.exception.AbortOperationException
import nextflow.util.Duration
/**
 * Holds the cloud configuration settings provided in the `nextflow.config` file
 * in the `cloud` scope
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class CloudConfig extends LaunchConfig {

    /**
     * Fetch the cloud config setting in the system environment
     * adding to each attribute the `NXF_CLOUD_` prefix
     */
    static class EnvConfig extends CascadingConfig<String,Object> {

        private Map env

        EnvConfig(Map env=null) {
            this.env = env != null ? env : System.getenv()
        }

        def getAttribute(String key) {
            final var = "NXF_CLOUD_${key.toUpperCase().replaceAll(/\./,'_')}".toString()
            return env .get(var)
        }
    }

    /**
     * Wrap the nextflow configuration settings in the cloud scope e.g.
     * <pre>
     *     cloud {
     *         instanceType = 'm4.xlarge'
     *         nextflow {
     *             version = '0.22.0'
     *             options = '-Dsome=var'
     *         }
     *     }
     * </pre>
     */
    static class Nextflow extends CascadingConfig<String,Object> {

        Nextflow( Map config ) {
            super(config, null)
        }

        @ConfigField
        String getPull() {
            getAttribute('pull', '')
        }

        @ConfigField
        String getVersion() {
            getAttribute('version', Const.APP_VER)
        }

        @ConfigField
        String getMode() {
            getAttribute('mode', 'ignite')
        }

        @ConfigField
        String getOptions() {
            getAttribute('options', '')
        }

        @ConfigField
        String getTrace() {
            getAttribute('trace', '')
        }
    }

    /**
     * Model the configuration of the cloud auto-scaling feature
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    static class Autoscale extends LaunchConfig {

        static final private Duration _5_MIN = Duration.of('5min')

        /**
         * Class factory method.
         *
         * @param cfg
         *      The underlying configuration object, either a {@link Boolean} value or an instance
         *      of {@link Map} holding the auto-scaling configuration settings
         * @param fallback
         *      A fallback {@link LaunchConfig} object that can be used to retrieve setting not available
         *      in the main scope
         * @return
         *      The {@link Autoscale} instance
         */
        static Autoscale create( cfg, LaunchConfig fallback=null) {

            if( cfg == null ) {
                return new Autoscale([enabled: false], fallback)
            }

            if( cfg instanceof Boolean ) {
                return new Autoscale([enabled: cfg], fallback)
            }

            if( cfg instanceof Map ) {

                if( !cfg.containsKey('enabled') )
                    cfg.enabled = true // <-- enable autoscale by default

                return new Autoscale(cfg, fallback)
            }

            throw new IllegalArgumentException("Not a valid autoscale config object: $cfg [${cfg.getClass().getName()}]")
        }


        protected Autoscale(Map config, LaunchConfig fallback) {
            super(config, fallback)
        }

        @Override
        String getRole() { ROLE_WORKER }

        @ConfigField
        boolean getEnabled() {
            getAttribute('enabled') as boolean
        }

        @ConfigField
        boolean getTerminateWhenIdle() {
            getAttribute('terminateWhenIdle') as boolean
        }

        @Deprecated
        @ConfigField
        String getTerminationPolicy() {
            getAttribute('terminationPolicy')
        }

        @ConfigField
        def Duration getIdleTimeout() {
            getAttribute('idleTimeout', _5_MIN) as Duration
        }

        @ConfigField
        Duration getStarvingTimeout()  {
            getAttribute('starvingTimeout', _5_MIN) as Duration
        }

        @ConfigField
        int getMinInstances() {
            getAttribute('minInstances', 1) as int
        }

        @ConfigField
        int getMaxInstances() {
            getAttribute('maxInstances', Integer.MAX_VALUE) as int
        }

    }


    private static final String DFLT_USER_NAME = 'ec2-user'

    static CloudConfig create( Map config ) {
        new CloudConfig( (Map)config?.cloud )
    }

    protected CloudConfig( Map config ) {
        this(config, new EnvConfig())
    }
    protected CloudConfig( Map config, Map env ) {
        this(config, new EnvConfig(env))
    }

    protected CloudConfig( Map config, CascadingConfig<String,Object> fallback ) {
        super( config ? new LinkedHashMap(config) : [:], fallback)
        // normalize auto-scale property
        if( this.config.containsKey('autoscale')) {
            def scale = this.config.get('autoscale')
            this.config.put('autoscale', Autoscale.create(scale,this))
        }
        // normalise `nextflow` config
        def nxf = (Map)this.config.get('nextflow') ?: [:]
        this.config.put('nextflow', new Nextflow(nxf))
    }

    CloudConfig build() {

        // when a `keyName` is given this is supposed to be a key provided by the cloud provider
        // thus the key will not be installed and the user not created
        if( keyName ) {
            if( !userName )
                setUserName(DFLT_USER_NAME)
        }
        else {
            createUser = true
            if( keyFile && !keyHash )
                setKeyHash( keyFile.text.trim() )

            if( !keyHash ) {
                def file = getUserPublicKeyFile()
                setKeyFile( file )
                setKeyHash( file.text.trim() )
            }

            if( !userName )
                userName = System.getProperty('user.name')
        }
        return this
    }

    CloudConfig validate(CloudDriver driver) {
        if( !imageId ) throw new IllegalStateException("Missing cloud config `imageId` attribute")
        if( !instanceType ) throw new IllegalStateException("Missing cloud config `instanceType` attribute")
        if( !clusterName ) throw new IllegalStateException("Missing cloud launch cluster name")
        if( !userName ) throw new IllegalStateException('Missing cloud launch user name')

        driver.validate(this)
        return this
    }

    private Path getUserPublicKeyFile() {
        def HOME = System.getProperty('user.home') as Path
        def rsaKey = HOME.resolve('.ssh/id_rsa.pub')

        if( rsaKey.exists() ) {
            return rsaKey
        }

        def dsaKey = HOME.resolve('.ssh/id_dsa.pub')
        if( dsaKey.exists() ) {
            return dsaKey
        }

        def message = '''
            Cannot find a SSH public key needed to access the cloud instances. You can either:
            - specify the `keyName` attribute in the configuration file with the key name
              given by the cloud provider
            - create a new key by using the command: `ssh-keygen` and store it in the default
              path i.e. $HOME/.ssh/id_rsa
            '''
                .stripIndent().leftTrim()
        throw new AbortOperationException(message)
    }

    private CloudConfig() {}


    CloudConfig setClusterName(String name) {
        setAttribute('clusterName', name)
        return this
    }

    CloudConfig setRole(String role) {
        setAttribute('role', role)
        return this
    }

    CloudConfig setCreateUser(boolean value) {
        setAttribute('createUser', value)
        return this
    }

    CloudConfig setUserName( String userName ) {
        setAttribute('userName', userName)
        return this
    }

    CloudConfig setKeyHash( String key ) {
        setAttribute('keyHash', key)
        return this
    }

    CloudConfig setKeyFile( Path file ) {
        if( file ) {
            setAttribute('keyFile', file.toString())
        }
        return this
    }

    @ConfigField(_private=true)
    Path getKeyFile() {
        getAttribute('keyFile') as Path
    }

    Path getPrivateKeyFile() {
        def file = getKeyFile()
        if( !file )
            return null
        def baseName = file.getBaseName()
        if( !baseName )
            return null

        def parent = file.getParent()
        parent ? parent.resolve(baseName) : Paths.get(baseName)
    }

    CloudConfig setDriverName( String driver ) {
        setAttribute('driver', driver)
        return this
    }

    CloudConfig setInstanceType(String instanceType) {
        setAttribute('instanceType', instanceType)
        return this
    }

    CloudConfig setImageId(String imageId) {
        setAttribute('imageId', imageId)
        return this
    }

    CloudConfig setSpotPrice(String price) {
        setAttribute('spotPrice', price)
        return this
    }

    CloudConfig setKeyName(String keyName) {
        setAttribute('keyName', keyName)
        return this
    }

    @ConfigField('driver')
    String getDriverName() {
        getAttribute('driver')
    }

    @ConfigField
    Autoscale getAutoscale() {
        (Autoscale)getAttribute('autoscale') ?: Autoscale.create(null)
    }

    String toString() {
        "imageId=$imageId, instanceType=$instanceType, spotPrice=$spotPrice"
    }

    /**
     * Render the cloud config in a human readable way
     *
     * @return A multiline string containing the config data
     */
    String prettyPrint() {
        formatEntry(this.copyPublicAttributes().sort())
    }

    @CompileDynamic
    private String formatEntry( def obj ) {
        def builder1 = new StringBuilder()

        if( obj instanceof Map ) {
            // `builder2` buffer is used only for composed values i.e. map objects that
            // need to be rendered at the bottom of the list
            def builder2 = new StringBuilder()

            ((Map)obj).sort().each { key, val ->
                if( val instanceof Map ) {
                    if( !val.isEmpty() )
                        builder2 << "- $key:\n${formatEntry(val).indent('  ')}"
                }
                else if( val instanceof CascadingConfig ) {
                    if( !val.isEmpty() )
                        builder2 << "- $key:\n${formatEntry(val.copyPublicAttributes()).indent('  ')}"
                }
                else if( key == 'keyHash' ) {
                    if( keyFile ) {
                        key = 'keyFile'
                        val = keyFile.toString()
                    }
                    else {
                        // resize the `hashKey` to just 40 chars
                        if(  val.toString().size()>40 ) { val=val.toString().substring(0,40)+'...'  }
                        val = formatEntry(val)
                    }
                    builder1 << "- $key: ${val}\n"
                }
                else {
                    builder1 << "- $key: ${formatEntry(val)}\n"
                }
            }

            return builder1.toString() + builder2
        }

        if( obj instanceof List ) {
            return obj.size()==1 ? formatEntry(obj.get(0)) : obj.collect{ formatEntry(it) } .toString()
        }

        needQuotes(obj) ? "'$obj'" : obj.toString()
    }

    private boolean needQuotes(value) {
        return !(value instanceof Boolean || value instanceof Number)
    }

}