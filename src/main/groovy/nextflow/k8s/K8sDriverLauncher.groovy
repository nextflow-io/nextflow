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

package nextflow.k8s
import java.nio.file.NoSuchFileException

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.google.common.hash.Hashing
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.cli.CmdRun
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.scm.AssetManager
import nextflow.scm.ProviderConfig
import nextflow.util.ConfigHelper
import org.codehaus.groovy.runtime.MethodClosure
/**
 * Configure and submit the execution of pod running the Nextlow main application
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class K8sDriverLauncher {

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
    private CmdRun cmd

    /**
     * Kubernetes client
     */
    private K8sClient client

    /**
     * Nextflow resolved config object
     */
    private Map config

    /**
     * Absolute path where the workflow temporary data is stored.
     * This must be a shared path mounted through a persistent volume claim.
     */
    private String workDir

    /**
     * Absolute path where workflow scripts and binary files are stored.
     * This must be a shared path mounted through a persistent volume claim.
     */
    private String projectDir

    /**
     * Absolute path where used data is stored
     *
     * This must be a shared path mounted through a persistent volume claim.
     */
    private String userDir

    /**
     * The current user name
     */
    private String userName = System.properties.get('user.name')

    private Map<String,String> configMounts = [:]

    private String paramsFile

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
        this.config = makeConfig(pipelineName)
        this.client = createK8sClient(config)
        checkStorageAndPaths()
        createK8sConfigMap()
        createK8sLauncherPod()
        printK8sPodOutput()
    }

    /**
     * Wait for the driver pod creation and prints the log to the
     * console standard output
     */
    protected void printK8sPodOutput() {
        final name = getPodName()
        print "Pod submitted: $name .. waiting to start"
        while( true ) {
            sleep 1000
            if( !client.podState(name).containsKey('waiting')  ) break
        }
        print "\33[2K\r"
        println "Pod started: $name"
        client.podLog(getPodName(), follow:true).eachLine { println it }
    }

    /**
     * Setup and verify required volumes and paths are properly configured
     */
    protected void checkStorageAndPaths() {
        final volumes = new VolumeClaims(config.k8s?.volumeClaims)
        if( !volumes )
            throw new AbortOperationException("Missing volume claim -- At least persistent volume claim definition needs to be provided in the nextflow configuration file")

        final defaultSharedDir = volumes.getFirstMount()

        // -- set launch dir
        userDir = config.k8s?.userDir ?: "$defaultSharedDir/$userName"
        if( !volumes.findVolumeByPath(userDir) )
            throw new AbortOperationException("Kubernetes `userDir` must be a path mounted as a persistent volume -- userDir=$userDir; volumes=${volumes.getMountPaths().join(', ')}")

        // -- set the work dir
        workDir = cmd.workDir ?: config.k8s?.workDir ?: config.workDir ?: "$userDir/work"
        if( !volumes.findVolumeByPath(workDir) )
            throw new AbortOperationException("Kubernetes workDir must be a path mounted as a persistent volume -- workDir=$workDir; volumes=${volumes.getMountPaths().join(', ')}")

        // -- set project dir
        projectDir = config.k8s?.projectDir ?: "$defaultSharedDir/projects"
        if( !volumes.findVolumeByPath(projectDir) )
            throw new AbortOperationException("Kubernetes projectDir must be a path mounted as a persistent volume -- projectDir=$projectDir; volumes=${volumes.getMountPaths().join(', ')}")


        log.debug "Kubernetes workDir=$workDir; projectDir=$projectDir; volumeClaims=$volumes"
    }

    /**
     * Retrieve the workflow configuration and merge with the current local one.
     *
     * @param pipelineName Workflow project name
     * @return A {@link Map} modeling the execution configuration settings
     */
    protected Map makeConfig(String pipelineName) {
        // -- check and parse project remote config
        final remote = new AssetManager(pipelineName, cmd)
                .checkValidRemoteRepo()
                .readRemoteConfig(cmd.profile)

        // -- load local config if available
        final local = new ConfigBuilder()
                .setOptions(cmd.launcher.options)
                .setCmdRun(cmd)
                .configObject()

        final result = (ConfigObject)remote.merge(local)

        // set k8s executor
        result.process.executor = 'k8s'
        // make sure to disable k8s auto mounts
        result.k8s.autoMountHostPaths = false
        // strip default work dir
        if( result.workDir == 'work' )
            result.remove('workDir')

        return result.toMap()
    }

    /**
     * Get a Kubernetes HTTP client to interact with the K8s cluster
     *
     * @param config The nextflow configuration with may hold the client settings
     * @return The {@link K8sClient} instance
     */
    protected K8sClient createK8sClient(Map config) {

        def k8sConfig = ( config.k8s instanceof Map && config.k8s?.client?.server
                ? ClientConfig.fromMap(config.k8s as Map)
                : ClientConfig.discover()
        )

        new K8sClient(k8sConfig)
    }

    /**
     * @return Container image name used by the driver pod
     */
    protected String getImageName() {
        final defImage = "nextflow/nextflow:${Const.APP_VER}"
        return config.navigate('k8s.nextflow.image', defImage)
    }

    /**
     * @return The driver pod name
     */
    protected String getPodName() {
        assert runName
        "nf-run-$runName".replaceAll('_','-')
    }

    private void checkUnsupportedOption(String name) {
        def field = cmd.getClass().getDeclaredField(name)
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
        def field = cmd.getClass().getDeclaredField(name)
        field.setAccessible(true)
        def val = field.get(cmd)
        if( ( eval ? eval(val) : val ) ) {
            def param = field.getAnnotation(Parameter)
            if( param ) {
                result << "${param.names()[0]} ${val}"
                return
            }

            param = field.getAnnotation(DynamicParameter)
            if( param && val instanceof Map ) {
                val.each { k,v ->
                    result << "${param.names()[0]}$k $v"
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

        def result = []
        // -- configure NF command line
        result << "nextflow"

        if( cmd.launcher.options.logFile )
            result << "-log $cmd.launcher.options.logFile"
        if( cmd.launcher.options.trace )
            result << "-trace ${cmd.launcher.options.trace.join(',')}"
        if( cmd.launcher.options.debug )
            result << "-debug ${cmd.launcher.options.debug.join(',')}"
        if( cmd.launcher.options.jvmOpts )
            cmd.launcher.options.jvmOpts.each { k,v -> result << "-D$k=$v" }

        result << "run"
        result << pipelineName

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
        addOption(result, cmd.&profile )
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

        // -- setup config file
        def cmd = ''
        cmd += '[ -f /etc/nextflow/scm ] && ln -s /etc/nextflow/scm $NXF_HOME/scm; '
        cmd +=  '[ -f /etc/nextflow/nextflow.config ] && cp /etc/nextflow/nextflow.config nextflow.config; '
        cmd += getLaunchCli()

        def params = [
                podName: getPodName(),
                imageName: getImageName(),
                command: ['/bin/bash', '-c', cmd],
                labels: [ app: 'nextflow', runName: runName ],
                workDir: userDir,
                namespace: client.config.namespace,
                volumeClaims: config.k8s.volumeClaims,
                configMounts: configMounts,
                env: [
                        NXF_WORK: workDir,
                        NXF_ASSETS: projectDir,
                        NXF_EXECUTOR: 'k8s'
                ],
        ]

        K8sHelper.createPodSpec(params)
    }

    /**
     * Creates and executes the nextflow driver pod
     * @return A {@link nextflow.k8s.client.K8sResponseJson} response object
     */
    protected createK8sLauncherPod() {
        final spec = makeLauncherSpec()
        client.podCreate(spec)
    }

    /**
     * Creates a K8s ConfigMap to share the nextflow configuration in the K8s cluster
     */
    protected void createK8sConfigMap() {
        Map<String,String> configMap = [:]
        if( config ) {
            configMap.'nextflow.config' = ConfigHelper.toCanonicalString(config)
        }

        def scm = ProviderConfig.getDefault()
        if( scm ) {
            configMap.'scm'  = ConfigHelper.toCanonicalString(scm)
        }

        if( cmd.paramsFile ) {
            final file = FileHelper.asPath(cmd.paramsFile)
            if( !file.exists() ) throw new NoSuchFileException("Params file does not exist: $file")
            configMap[ file.getName() ] = file.text
            paramsFile = "/etc/nextflow/$file.name"
        }

        if( configMap ) {
            def name = "nf-config-${hash(configMap.values())}"
            tryCreateConfigMap(name, configMap)
            configMounts[name] = '/etc/nextflow'
        }
    }

    protected void tryCreateConfigMap(String name, Map data) {
        try {
            client.configCreate(name, data)
        }
        catch( K8sResponseException e ) {
            if( e.response.reason != 'AlreadyExists' )
                throw e
        }
    }

    protected String hash(Collection<String> text) {
        def hasher = Hashing .murmur3_32() .newHasher()
        def itr = text.iterator()
        while( itr.hasNext() ) {
            hasher.putUnencodedChars(itr.next())
        }

        return hasher.hash().toString()
    }

}
