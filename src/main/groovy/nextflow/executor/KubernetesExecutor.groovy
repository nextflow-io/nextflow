/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow.executor
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger

import groovy.util.logging.Slf4j
import nextflow.container.DockerBuilder
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.MemoryUnit
import nextflow.util.PathTrie
import nextflow.util.ServiceName

/**
 * Kubernetes executor
 *
 * See http://kubernetes.io
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName('k8s')
class KubernetesExecutor extends AbstractGridExecutor {

    static final public String CMD_K8S = '.command.yml'

    final protected BashWrapperBuilder createBashWrapperBuilder(TaskRun task) {

        if( !task.config.container ) {
            throw new AbortOperationException("Missing container for process `${task.processor.name}`")
        }

        // creates the wrapper script
        final builder = new KubernetesWrapperBuilder(task)
        return builder
    }


    @Override
    protected String getHeaderToken() {
        throw new UnsupportedOperationException()
    }

    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result) {
        throw new UnsupportedOperationException()
    }

    /**
     * Define the Kubernates job execution command line
     *
     * @param task A {@link TaskRun} instance that need to be submitted for execution
     * @param scriptFile (not used)
     * @return the command line to submit the task execution
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile) {
        return ['kubectl', 'create', '-f', CMD_K8S, '-o', 'name']
    }

    /**
     * Parse the output of a `kubectl create <job>` to find out the job ID
     * @param text The `kubectl` standard output
     * @return The newly submit job ID
     */
    @Override
    def parseJobId(String text) {
        final str = text.trim()
        if( str?.startsWith('pod/') ) {
            return str.substring(4)
        }
        throw new IllegalStateException("Not a valid Kubernates job id: `$text`")
    }

    @Override
    protected List<String> getKillCommand() {
        ['kubectl', 'delete', 'pod']
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        ['kubectl', 'get', 'pods', '-a']
    }

    static private Map DECODE_STATUS = [
            'Pending': QueueStatus.PENDING,
            'Running': QueueStatus.RUNNING,
            'Completed': QueueStatus.DONE,
            'Error': QueueStatus.ERROR,
            'ContainerCreating': QueueStatus.RUNNING,
            'ContainerCannotRun': QueueStatus.ERROR,
            'ErrImagePull': QueueStatus.ERROR,
            'ImagePullBackOff': QueueStatus.ERROR,
    ]

    static private final K8S_JOB_ID = ~/^(nxf-[0-9a-f]{32})/

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {
        def result = [:]
        if( !text ) return result

        text.readLines().each { line ->
            def items = line.tokenize()
            def regex = K8S_JOB_ID.matcher(items[0])
            if( regex.matches() ) {
                def id = regex.group(1)
                def status = items[2]
                result[id] = DECODE_STATUS[ status ]
            }
        }

        return result
    }


    /**
     * Extends {@link BashWrapperBuilder} to handle Kubernetes
     * job description YAML file
     */
    static class KubernetesWrapperBuilder extends BashWrapperBuilder {

        static private volumes = new AtomicInteger()

        private String taskHash
        private cpu
        private mem

        KubernetesWrapperBuilder(TaskRun task) {
            super(task)
            taskHash = task.hash.toString()
            cpu = task.config.getCpus()
            mem = task.config.getMemory()
        }

        /* only for tests -- don not use */
        protected KubernetesWrapperBuilder(TaskBean bean) {
            super(bean)
        }

        /**
         * Container execution is managed by Kubernetes thus disable
         * the execution through the job wrapper script
         *
         * @return {@code false}
         */
        @Override
        protected boolean containerInit() {
            return false
        }

        @Override
        protected boolean fixOwnership() {
            containerConfig.fixOwnership
        }

        /**
         * Creates the job BASH wrapper file and the Kubernetes YAML descriptor
         * @return The {@link Path} to the BASH wrapper file
         */
        Path build() {
            final wrapper = super.build()
            // save the condor manifest
            workDir.resolve(CMD_K8S).text = makeYaml()
            // returns the launcher wrapper file
            return wrapper
        }

        /**
         * @return The Kubernetes job descriptor YAML string
         */
        private String makeYaml() {

            // get input files paths
            def paths = DockerBuilder.inputFilesToPaths(inputFiles)
            // add standard paths
            if( binDir ) paths << binDir
            if( workDir ) paths << workDir

            def trie = new PathTrie()
            paths.each { trie.add(it) }

            // defines the mounts
            def mounts = [:]
            trie.longest().each {
                mounts.put( "vol-${volumes.incrementAndGet()}", it )
            }

            new YamlBuilder(
                    name: "nxf-${taskHash}",
                    image: containerImage,
                    cmd: ['bash', TaskRun.CMD_RUN],
                    workDir: workDir.toString(),
                    mounts: mounts,
                    cpu: cpu,
                    mem: mem )
                    .create()
        }

    }

    /**
     * Kubernetes job YAML descriptor build helper
     */
    static class YamlBuilder {

        String name
        String image
        List<String> cmd
        String workDir
        Map<String,String> mounts
        int cpu
        MemoryUnit mem

        String create() {

"""\
apiVersion: v1
kind: Pod
metadata:
  name: $name
  labels:
    app: nextflow
spec:
  restartPolicy: Never
  containers:
  - name: $name
    image: $image
    command: ${cmd.collect { "\"$it\"" }}
    workingDir: $workDir
${getResources0()}${getMounts0(mounts)}${getVolumes0(mounts)}\
"""

        }

        /**
         * Creates the kubernetes submit YAML mounts fragment
         * @param vols A map holding the volumes to be mount
         * @return A YAML fragment representing the container mounts
         */
        private String getMounts0(Map<String,String> vols) {

            if( !vols ) return ''

            def result = new StringBuilder('    volumeMounts:\n')

            vols.each { String name, String path ->

                result.append(
"""\
    - mountPath: $path
      name: $name
"""
                )
            }

            return result.toString()
        }

        /**
         * Creates the kubernetes submit YAML volumes fragment
         * @param vols A map holding the volumes to be mount
         * @return A YAML fragment representing the container volumes
         */
        private String getVolumes0(Map<String,String> vols) {

            if( !vols ) return ''

            def result = new StringBuilder('  volumes:\n')

            vols.each { String name, String path ->

                result.append(
"""\
  - name: $name
    hostPath:
      path: $path
"""
                )

            }

            result.toString()
        }


        private String getResources0() {

            if( cpu>1 || mem ) {
                def res = ''
                if( cpu>1 ) res += "        cpu: ${cpu}\n"
                if( mem ) res += "        memory: ${mem.toMega()}Mi\n"

                def result = new StringBuilder('    resources:\n')
                result.append('      limits:\n')
                result.append(res)
                result.append('      requests:\n')
                result.append(res)
                return result.toString()
            }
            else {
                return ''
            }
        }

    }
}
