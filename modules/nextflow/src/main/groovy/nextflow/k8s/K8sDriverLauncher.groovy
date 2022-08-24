/*
 * Copyright 2020-2022, Seqera Labs
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
import nextflow.k8s.model.ResourceType
import nextflow.plugin.Plugins
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
     * Either a Pod or Job
     */
    private ResourceType resourceType = ResourceType.Pod

    /**
     * Container image to be used for the Nextflow driver pod
     */
    private String headImage

    /** 
     * Request CPUs to be used for the Nextflow driver pod
     */
    private int headCpus

    /** 
     * Request memory to be used for the Nextflow driver pod
     */
    private String headMemory

    /** 
     * Pre-script to run before nextflow
     */
    private String headPreScript

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
    private ConfigObject config

    /**
     * Name of the config map used to propagate the nextflow
     * setting in the container
     */
    private String configMapName

    /**
     * Kubernetes specific config settings
     */
    private K8sConfig k8sConfig

    private String paramsFile

    private boolean interactive

    /**
     * Runs in background mode
     */
    private boolean background

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
        if( background && interactive )
            throw new AbortOperationException("Option -bg conflicts with interactive mode")
        this.config = makeConfig(pipelineName)
        this.k8sConfig = makeK8sConfig(config.toMap())
        this.k8sClient = makeK8sClient(k8sConfig)
        this.k8sConfig.checkStorageAndPaths(k8sClient)
        createK8sConfigMap()
        createK8sLauncherPod()
        waitPodStart()
        // login into container session
        if( interactive )
            launchLogin()
        // dump pod output
        else if( !background )
            printK8sPodOutput()
        else
            log.debug "Nextflow driver launched in background mode -- pod: $runName"
        waitPodEnd()
    }

    int shutdown() {
        if( background )
            return 0
        // fetch the container exit status
        final exitCode = waitPodTermination()
        // cleanup the config map if OK
        def deleteOnSuccessByDefault = exitCode==0
        if( k8sConfig.getCleanup(deleteOnSuccessByDefault)  ) {
            deleteConfigMap()
        }
        return exitCode
    }

    protected void waitPodEnd() {
        if( background )
            return
        final currentState = k8sConfig.useJobResource() ? k8sClient.jobState(runName) : k8sClient.podState(runName)
        if (currentState && currentState?.running instanceof Map) {
            final name = runName
            println "${resourceType} running: $name ... waiting for ${resourceType.lower()} to stop running"
            try {
                while( true ) {
                    sleep 10000
                    final state = k8sConfig.useJobResource() ? k8sClient.jobState(name) : k8sClient.podState(name)
                    if ( state && !(state?.running instanceof Map) )  {
                        println "${resourceType} $name has changed from running state $state"
                        break
                    }
                }
            }
            catch( Exception e ) {
                log.warn "Caught exception while waiting for ${resourceType.lower()} to stop running"
            }
        }
    }
    
    protected boolean isWaitTimedOut(long time) {
        System.currentTimeMillis()-time > 90_000
    }

    protected int waitPodTermination() {
        log.debug "Wait for ${resourceType.lower()} termination name=$runName"
        final rnd = new Random()
        final time = System.currentTimeMillis()
        Map state = null
        try {
            while( true ) {
                sleep rnd.nextInt(500)
                state = k8sConfig.useJobResource() ? k8sClient.jobState(runName) : k8sClient.podState(runName)
                if( state?.terminated instanceof Map  )
                    return state.terminated.exitCode as int

                else if( isWaitTimedOut(time) )
                    throw new IllegalStateException("Timeout waiting for ${resourceType.lower()} terminated state=$state")
            }
        }
        catch( Exception e ) {
            log.warn "Unable to fetch ${resourceType.lower()} exit status -- ${resourceType.lower()}=$runName state=$state"
            return 127
        }
    }

    protected void deleteConfigMap() {
        try {
            k8sClient.configDelete(configMapName)
            log.debug "Deleted K8s configMap with name: $configMapName"
        }
        catch ( Exception e ) {
            log.warn "Unable to delete configMap: $configMapName", e
        }
    }

    protected void waitPodStart() {
        final name = runName
        print "${resourceType} submitted: $name .. waiting to start"
        while( true ) {
            sleep 1000
            final state = k8sConfig.useJobResource() ? k8sClient.jobState(name) : k8sClient.podState(name)
            if( state && !state.containsKey('waiting')  ) {
                break
            }
        }
        print "\33[2K\r"
        println "${resourceType} started: $name"
    }

    /**
     * Wait for the driver pod creation and prints the log to the
     * console standard output
     */
    protected void printK8sPodOutput() {
        if ( k8sConfig.useJobResource() )
            k8sClient.jobLog(runName, follow:true).eachLine { println it }
        else
            k8sClient.podLog(runName, follow:true).eachLine { println it }
    }

    protected ConfigObject loadConfig( String pipelineName ) {

        // -- load local config if available
        final builder = new ConfigBuilder()
                .setShowClosures(true)
                .setOptions(cmd.launcher.options)
                .setProfile(cmd.profile)
                .setCmdRun(cmd)

        if( !interactive && !pipelineName.startsWith('/') && !cmd.remoteProfile && !cmd.runRemoteConfig ) {
            // -- check and parse project remote config
            Plugins.init()
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
    protected ConfigObject makeConfig(String pipelineName) {

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

        log.trace "K8s config object:\n${ConfigHelper.toCanonicalString(config).indent('  ')}"
        return config
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
            log.warn "Unknown command-line option to check: $name"
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
        addOption(result, cmd.&entryName )

        if( paramsFile ) {
            result << "-params-file $paramsFile"
        }

        if ( cmd.runRemoteConfig )
            cmd.runRemoteConfig.forEach { result << "-config $it" }

        if ( cmd.remoteProfile )
            result << "-profile ${cmd.remoteProfile}"

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
        PodSpecBuilder builder = new PodSpecBuilder()
            .withPodName(runName)
            .withImageName(headImage ?: k8sConfig.getNextflowImageName())
            .withCommand(['/bin/bash', '-c', cmd])
            .withLabels([ app: 'nextflow', runName: runName ])
            .withNamespace(k8sClient.config.namespace)
            .withServiceAccount(k8sClient.config.serviceAccount)
            .withPodOptions(k8sConfig.getPodOptions())
            .withEnv( PodEnv.value('NXF_WORK', k8sConfig.getWorkDir()) )
            .withEnv( PodEnv.value('NXF_ASSETS', k8sConfig.getProjectDir()) )
            .withEnv( PodEnv.value('NXF_EXECUTOR', 'k8s'))
            .withEnv( PodEnv.value('NXF_ANSI_LOG', 'false'))
            .withMemory(headMemory?:"")
            .withCpus(headCpus)

        if ( k8sConfig.useJobResource()) {
            this.resourceType = ResourceType.Job
            return builder.buildAsJob()
        }
        else {
            return builder.build()
        }

        // note: do *not* set the work directory because it may need to be created  by the init script
    }

    /**
     * Creates and executes the nextflow driver pod
     * @return A {@link nextflow.k8s.client.K8sResponseJson} response object
     */
    protected createK8sLauncherPod() {
        final spec = makeLauncherSpec()
        if ( k8sConfig.useJobResource() ) {
            k8sClient.jobCreate(spec, yamlDebugPath())
        } else {
            k8sClient.podCreate(spec, yamlDebugPath())
        }
    }

    protected Path yamlDebugPath() {
        boolean debug = config.k8s.debug?.yaml?.toString() == 'true'
        final result = debug ? Paths.get(".nextflow.${resourceType.lower()}.yaml") : null
        if( result )
            log.info "Launcher ${resourceType.lower()} spec file: $result"
        return result
    }

    protected Path getScmFile() {
        ProviderConfig.getScmConfigPath()
    }

    protected String getPipelineName() {
        return pipelineName
    }

    protected boolean getInteractive() {
        return interactive
    }

    protected ConfigObject getConfig() {
        return config
    }

    protected K8sConfig getK8sConfig() {
        return k8sConfig
    }

    protected K8sClient getK8sClient() {
        return k8sClient
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
        if( headPreScript ) 
            initScript += "[ -f '$headPreScript' ] && '$headPreScript'; "
        configMap['init.sh'] = initScript

        // nextflow config file
        if( this.config ) {
            configMap['nextflow.config'] = ConfigHelper.toCanonicalString( this.config )
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
        configMapName = makeConfigMapName(configMap)
        tryCreateConfigMap(configMapName, configMap)
        log.debug "Created K8s configMap with name: $configMapName"
        k8sConfig.getPodOptions().getMountConfigMaps().add( new PodMountConfig(configMapName, '/etc/nextflow') )
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
        if( result == 0 ) {
            if ( k8sConfig.useJobResource() )
                k8sClient.jobDelete(runName)
            else
                k8sClient.podDelete(runName)
        }
    }
}
