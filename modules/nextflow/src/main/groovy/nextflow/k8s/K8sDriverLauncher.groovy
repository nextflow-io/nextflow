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

package nextflow.k8s

import java.lang.reflect.Field
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.google.common.hash.Hashing
import groovy.util.logging.Slf4j
import nextflow.cli.CmdKubeRun
import nextflow.cli.CmdRun
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.model.PodEnv
import nextflow.k8s.model.PodMountConfig
import nextflow.k8s.model.PodSpecBuilder
import nextflow.scm.AssetManager
import nextflow.scm.ProviderConfig
import nextflow.util.ConfigHelper
import nextflow.util.Escape
import org.codehaus.groovy.runtime.MethodClosure
/**
 * Configure and submit the execution of pod running the Nextflow main application
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class K8sDriverLauncher {

    /**
     * Container image to be used for the Nextflow driver pod
     */
    private String podImage

    /**
     * Nextflow execution run name
     */
    private String runName

    /**
     * Workflow project to launch
     */
    private String pipelineName

    /**
     * Command run options
     */
    private CmdKubeRun cmd

    /**
     * Kubernetes client
     */
    private K8sClient k8sClient

    /**
     * Nextflow resolved config object
     */
    private Map config

    /**
     * Kubernetes specific config settings
     */
    private K8sConfig k8sConfig

    private String paramsFile

    private boolean interactive

    /**
     * Workflow script positional parameters
     */
    private List<String> args

    /**
     * Launcher entry point. Set-up the environment and create a pod that run the Nextflow
     * application (which in turns executed each task as a pod)
     *
     * @param name Workflow project entry name
     * @param args Workflow script positional parameters
     */
    void run(String name, List<String> args) {
        this.args = args
        this.pipelineName = name
        this.interactive = name == 'login'
        this.config = makeConfig(pipelineName)
        this.k8sConfig = makeK8sConfig(config)
        this.k8sClient = makeK8sClient(k8sConfig)
        this.k8sConfig.checkStorageAndPaths(k8sClient)
        createK8sConfigMap()
        createK8sLauncherPod()
        waitPodStart()
        interactive ? launchLogin() : printK8sPodOutput()
    }

    protected void waitPodStart() {
        final name = runName
        print "Pod submitted: $name .. waiting to start"
        while( true ) {
            sleep 1000
            final state = k8sClient.podState(name)
            if( state && !state.containsKey('waiting')  ) {
                break
            }
        }
        print "\33[2K\r"
        println "Pod started: $name"
    }

    /**
     * Wait for the driver pod creation and prints the log to the
     * console standard output
     */
    protected void printK8sPodOutput() {
        k8sClient.podLog(runName, follow:true).eachLine { println it }
    }

    protected ConfigObject loadConfig( String pipelineName ) {

        // -- load local config if available
        final builder = new ConfigBuilder()
                .setShowClosures(true)
                .setOptions(cmd.launcher.options)
                .setProfile(cmd.profile)
                .setCmdRun(cmd)

        if( !interactive && !pipelineName.startsWith('/') ) {
            // -- check and parse project remote config
            final pipelineConfig = new AssetManager(pipelineName, cmd) .getConfigFile()
            builder.setUserConfigFiles(pipelineConfig)
        }

        return builder.buildConfigObject()
    }

    protected K8sConfig makeK8sConfig(Map config) {
        config.k8s instanceof Map ? new K8sConfig(config.k8s as Map) : new K8sConfig()
    }

    protected makeK8sClient( K8sConfig k8sConfig ) {
        new K8sClient(k8sConfig.getClient())
    }

    /**
     * Retrieve the workflow configuration and merge with the current local one.
     *
     * @param pipelineName Workflow project name
     * @return A {@link Map} modeling the execution configuration settings
     */
    protected Map makeConfig(String pipelineName) {

        def file = new File(pipelineName)
        if( !interactive && file.exists() ) {
            def message = "The k8s executor cannot run local ${file.directory ? 'project' : 'script'}: $pipelineName"
            message += " -- provide the absolute path of a project available in the Kubernetes cluster or the URL of a project hosted in a Git repository"
            throw new AbortOperationException(message)
        }

        def config = loadConfig(pipelineName)

        // normalize pod entries
        def k8s = config.k8s

        if( !k8s.isSet('pod') )
            k8s.pod = []
        else if( k8s.pod instanceof Map ) {
            k8s.pod = [ k8s.pod ]
        }
        else if( !(k8s.pod instanceof List) )
            throw new IllegalArgumentException("Illegal k8s.pod configuratun value: ${k8s.pod}")

        // -- use the volume claims specified in the command line
        //   to populate the pod config
        for( int i=0; i<cmd.volMounts?.size(); i++ ){
            def entry = cmd.volMounts.get(i)
            def parts = entry.tokenize(':')
            def name = parts[0]
            def path = parts[1]
            if( i==0 ) {
                k8s.storageClaimName = name
                k8s.storageMountPath = path
            }
            else {
                k8s.pod.add( [volumeClaim: name, mountPath: path] )
            }
        }

        // -- backward compatibility
        if( k8s.isSet('volumeClaims') ) {
            log.warn "Config setting k8s.volumeClaims has been deprecated -- Use k8s.storageClaimName and k8s.storageMountPath instead"
            k8s.volumeClaims.each { k,v ->
                def name = k as String
                def path = v instanceof Map ? v.mountPath : v.toString()
                if( !k8s.isSet('storageClaimName') ) {
                    k8s.storageClaimName = name
                    k8s.storageMountPath = path
                }
                else if( !cmd.volMounts ) {
                    k8s.pod.add( [volumeClaim: name, mountPath: path] )
                }
            }
            // remove it
            k8s.remove('volumeClaims')
        }

        // -- set k8s executor
        config.process.executor = 'k8s'

        // -- strip default work dir
        if( config.workDir == 'work' )
            config.remove('workDir')

        // -- check work dir
        if( cmd?.workDir )
            k8s.workDir = cmd.workDir
        else if( !k8s.isSet('workDir') && config.workDir )
            k8s.workDir = config.workDir

        // -- some cleanup
        if( !k8s.pod )
            k8s.remove('pod')

        if( !k8s.storageClaimName )
            k8s.remove('storageClaimName')
        if( !k8s.storageMountPath )
            k8s.remove('storageMountPath')

        if( !config.libDir )
            config.remove('libDir')
        
        final result = config.toMap()
        log.trace "K8s config object:\n${ConfigHelper.toCanonicalString(result).indent('  ')}"
        return result
    }


    private Field getField(CmdRun cmd, String name) {
        def clazz = cmd.class
        while( clazz != CmdRun ) {
            clazz = cmd.class.getSuperclass()
        }
        clazz.getDeclaredField(name)
    }

    private void checkUnsupportedOption(String name) {
        def field = getField(cmd,name)
        if( !field ) {
            log.warn "Unknown cli option to check: $name"
            return
        }
        field.setAccessible(true)
        if( field.get(cmd) ) {
            def param = field.getAnnotation(Parameter)
            def opt = param.names() ? param.names()[0] : "-$name"
            abort(opt)
        }
    }

    private void abort(String opt) {
        throw new AbortOperationException("Option `$opt` not supported with Kubernetes deployment")
    }

    private void unsupportedCliOptions(MethodClosure... fields) {
        unsupportedCliOptions( fields.collect { it.getMethod()}  )
    }

    private void unsupportedCliOptions(List<String> names) {
        for( String x : names ) {
            checkUnsupportedOption(x)
        }
    }

    private void addOption(List result, MethodClosure m, Closure eval=null) {
        def name = m.getMethod()
        def field = getField(cmd,name)
        field.setAccessible(true)
        def val = field.get(cmd)
        if( ( eval ? eval(val) : val ) ) {
            def param = field.getAnnotation(Parameter)
            if( param ) {
                result << "${param.names()[0]} ${Escape.wildcards(String.valueOf(val))}"
                return
            }

            param = field.getAnnotation(DynamicParameter)
            if( param && val instanceof Map ) {
                val.each { k,v ->
                    result << "${param.names()[0]}$k ${Escape.wildcards(String.valueOf(v))}"
                }
            }
        }
    }

    /**
     * @return The nextflow driver command line
     */
    protected String getLaunchCli() {
        assert cmd
        assert pipelineName

        if( interactive ) {
            return "tail -f /dev/null"
        }

        def result = []
        // -- configure NF command line
        result << "nextflow"

        if( cmd.launcher.options.trace )
            result << "-trace ${cmd.launcher.options.trace.join(',')}"
        if( cmd.launcher.options.debug )
            result << "-debug ${cmd.launcher.options.debug.join(',')}"
        if( cmd.launcher.options.jvmOpts )
            cmd.launcher.options.jvmOpts.each { k,v -> result << "-D$k=$v" }

        result << "run"
        result << pipelineName

        if( runName )
            result << '-name' << runName

        addOption(result, cmd.&cacheable, { it==false } )
        addOption(result, cmd.&resume )
        addOption(result, cmd.&poolSize )
        addOption(result, cmd.&pollInterval )
        addOption(result, cmd.&queueSize)
        addOption(result, cmd.&revision )
        addOption(result, cmd.&latest )
        addOption(result, cmd.&withTrace )
        addOption(result, cmd.&withTimeline )
        addOption(result, cmd.&withDag )
        addOption(result, cmd.&dumpHashes )
        addOption(result, cmd.&dumpChannels )
        addOption(result, cmd.&env )
        addOption(result, cmd.&process )
        addOption(result, cmd.&params )

        if( paramsFile ) {
            result << "-params-file $paramsFile"
        }

        if( cmd.process?.executor )
            abort('process.executor')

        unsupportedCliOptions(
                cmd.&libPath,
                cmd.&test,
                cmd.&executorOptions,
                cmd.&stdin,
                cmd.&withSingularity,
                cmd.&withDocker,
                cmd.&withoutDocker,
                cmd.&withMpi,
                cmd.&clusterOptions,
                cmd.&exportSysEnv
        )

        if( args )
            result.add(args)

        return result.join(' ')
    }

    /**
     * @return A {@link Map} modeling driver pod specification
     */
    protected Map makeLauncherSpec() {
        assert runName
        assert k8sClient

        // -- setup config file
        String cmd = "source /etc/nextflow/init.sh; ${getLaunchCli()}"

        // create the launcher pod
        new PodSpecBuilder()
            .withPodName(runName)
            .withImageName(podImage ?: k8sConfig.getNextflowImageName())
            .withCommand(['/bin/bash', '-c', cmd])
            .withLabels([ app: 'nextflow', runName: runName ])
            .withNamespace(k8sClient.config.namespace)
            .withServiceAccount(k8sClient.config.serviceAccount)
            .withPodOptions(k8sConfig.getPodOptions())
            .withEnv( PodEnv.value('NXF_WORK', k8sConfig.getWorkDir()) )
            .withEnv( PodEnv.value('NXF_ASSETS', k8sConfig.getProjectDir()) )
            .withEnv( PodEnv.value('NXF_EXECUTOR', 'k8s'))
            .build()

        // note: do *not* set the work directory because it may need to be created  by the init script
    }

    /**
     * Creates and executes the nextflow driver pod
     * @return A {@link nextflow.k8s.client.K8sResponseJson} response object
     */
    protected createK8sLauncherPod() {
        final spec = makeLauncherSpec()
        k8sClient.podCreate(spec, yamlDebugPath())
    }

    protected Path yamlDebugPath() {
        boolean debug = config.k8s.debug?.yaml?.toString() == 'true'
        final result = debug ? Paths.get('.nextflow.pod.yaml') : null
        if( result )
            log.info "Launcher pod spec file: $result"
        return result
    }

    protected Path getScmFile() {
        ProviderConfig.SCM_FILE.toPath()
    }

    /**
     * Creates a K8s ConfigMap to share the nextflow configuration in the K8s cluster
     */
    protected void createK8sConfigMap() {
        Map<String,String> configMap = [:]

        final launchDir = k8sConfig.getLaunchDir()
        // init file
        String initScript = ''
        initScript += "mkdir -p '$launchDir'; if [ -d '$launchDir' ]; then cd '$launchDir'; else echo 'Cannot create directory: $launchDir'; exit 1; fi; "
        initScript += '[ -f /etc/nextflow/scm ] && ln -s /etc/nextflow/scm $NXF_HOME/scm; '
        initScript += '[ -f /etc/nextflow/nextflow.config ] && cp /etc/nextflow/nextflow.config $PWD/nextflow.config; '
        configMap['init.sh'] = initScript

        // nextflow config file
        if( config ) {
            configMap['nextflow.config'] = ConfigHelper.toCanonicalString(config)
        }

        // scm config file
        final scmFile = getScmFile()
        if( scmFile.exists() ) {
            configMap['scm']  = scmFile.text
        }

        // params file
        if( cmd.paramsFile ) {
            final file = FileHelper.asPath(cmd.paramsFile)
            if( !file.exists() ) throw new NoSuchFileException("Params file does not exist: $file")
            configMap[ file.getName() ] = file.text
            paramsFile = "/etc/nextflow/$file.name"
        }

        // create the config map
        final name = makeConfigMapName(configMap)
        tryCreateConfigMap(name, configMap)
        k8sConfig.getPodOptions().getMountConfigMaps().add( new PodMountConfig(name, '/etc/nextflow') )
    }

    protected void tryCreateConfigMap(String name, Map data) {
        try {
            k8sClient.configCreate(name, data)
        }
        catch( K8sResponseException e ) {
            if( e.response.reason != 'AlreadyExists' )
                throw e
        }
    }

    protected String makeConfigMapName( Map configMap ) {
        "nf-config-${hash(configMap.values())}"
    }

    protected String hash(Collection<String> text) {
        def hasher = Hashing .murmur3_32() .newHasher()
        def itr = text.iterator()
        while( itr.hasNext() ) {
            hasher.putUnencodedChars(itr.next())
        }

        return hasher.hash().toString()
    }

    protected void launchLogin() {
        def launchDir = k8sConfig.getLaunchDir()
        def cmd = "kubectl -n ${k8sClient.config.namespace} exec -it $runName -- /bin/bash -c 'cd $launchDir; exec bash --login -i'"
        def proc = new ProcessBuilder().command('bash','-c',cmd).inheritIO().start()
        def result = proc.waitFor()
        if( result == 0 )
            k8sClient.podDelete(runName)
    }
}
