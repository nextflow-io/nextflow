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

package nextflow.executor

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.k8s.K8sExecutor
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ScriptType
import nextflow.util.ServiceDiscover
import nextflow.util.ServiceName
/**
 * Helper class to create {@link Executor} objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ExecutorFactory {

    static public String DEFAULT_EXECUTOR = System.getenv('NXF_EXECUTOR') ?: 'local'

    /*
     * Map the executor class to its 'friendly' name
     */
    final static Map<String, Class<? extends Executor>> BUILT_IN_EXECUTORS = [
            'nope': NopeExecutor,
            'local': LocalExecutor,
            'sge':  SgeExecutor,
            'oge':  SgeExecutor,
            'uge':  SgeExecutor,
            'lsf': LsfExecutor,
            'pbs': PbsExecutor,
            'pbspro': PbsProExecutor,
            'slurm': SlurmExecutor,
            'crg': CrgExecutor,
            'bsc': LsfExecutor,
            'condor': CondorExecutor,
            'k8s': K8sExecutor,
            'nqsii': NqsiiExecutor,
            'moab': MoabExecutor,
            'oar': OarExecutor
    ]

    @PackageScope Map<String, Class<? extends Executor>> executorsMap

    private Map<Class<? extends Executor>,? extends Executor> executors = new HashMap<>()

    @PackageScope Map<Class<? extends Executor>,? extends Executor> getExecutors() { executors }

    ExecutorFactory() {
        executorsMap = new HashMap(20)
        // add built-in executors
        executorsMap.putAll(BUILT_IN_EXECUTORS)
        // discover non-core executors
        for( Class<Executor> clazz : ServiceDiscover.load(Executor) ) {
            log.trace "Discovered executor class: ${clazz.toString()}"
            final name = findNameByClass(clazz)
            final current = executorsMap.get(name)
            if( current ) {
                if( current.getAnnotation(ServiceName)?.important() ) {
                    log.debug "Executor ${current.getSimpleName()} has priority - skipping ${clazz}"
                    continue
                }
                log.debug "Replacing executor ${current.getSimpleName()} with ${clazz}"
            }
            executorsMap.put(name, clazz)
        }
    }

    String getDisplayName(String key) {
        final clazz = executorsMap.get(key)
        if( !clazz ) return key
        final exec = this.executors.get(clazz)
        if( !exec ) return key
        exec.getDisplayName() ?: key
    }

    /**
     * Extract the executor name by using the annotation {@code ServiceName} or fallback to simple classname
     * if the annotation is not provided
     *
     * @param clazz
     * @return
     */
    static String findNameByClass( Class<Executor> clazz ) {
        def annotation = clazz.getAnnotation(ServiceName)
        if( annotation )
            return annotation.value()

        def name = clazz.getSimpleName().toLowerCase()
        if( name.endsWith('executor') ) {
            name = name.subSequence(0, name.size()-'executor'.length())
        }

        return name
    }

    protected Class<? extends Executor> getExecutorClass(String executorName) {
        log.debug ">> processorType: '$executorName'"
        if( !executorName )
            return LocalExecutor

        def clazz =  executorsMap[executorName.toLowerCase()]
        if( !clazz )
            throw new IllegalArgumentException("Unknown executor name: $executorName")

        if( clazz instanceof Class )
            return clazz

        if( !(clazz instanceof String ) )
            throw new IllegalArgumentException("Not a valid executor class object: $clazz")

        // if the className is empty (because the 'processorType' does not map to any class, fallback to the 'processorType' itself)
        if( !clazz ) {
            clazz = executorName
        }

        log.debug "Loading executor class: ${clazz}"
        try {
            Thread.currentThread().getContextClassLoader().loadClass(clazz as String) as Class<Executor>
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Cannot find a valid class for specified executor: '${executorName}'")
        }

    }


    protected boolean isTypeSupported( ScriptType type, executor ) {

        if( executor instanceof Executor ) {
            executor = executor.class
        }

        if( executor instanceof Class ) {
            def annotation = executor.getAnnotation(SupportedScriptTypes)
            if( !annotation )
                throw new IllegalArgumentException("Specified argument is not a valid executor class: $executor -- Missing 'SupportedScriptTypes' annotation")

            return type in annotation.value()
        }

        throw new IllegalArgumentException("Specified argument is not a valid executor class: $executor")
    }

    Executor getExecutor(String processName, ProcessConfig processConfig, BodyDef script, Session session ) {
        // -- load the executor to be used
        def name = getExecutorName(processConfig,session) ?: DEFAULT_EXECUTOR
        def clazz = getExecutorClass(name)

        if( !isTypeSupported(script.type, clazz) ) {
            log.warn "Process '$processName' cannot be executed by '$name' executor -- Using 'local' executor instead"
            name = 'local'
            clazz = LocalExecutor.class
        }

        // this code is not supposed to be executed parallel
        def result = executors.get(clazz)
        if( result )
            return result

        result = createExecutor(clazz, name, session)
        executors.put(clazz, result)
        return result
    }

    protected Executor createExecutor( Class<? extends Executor> clazz, String name, Session session) {
        def result = clazz.newInstance()
        result.session = session
        result.name = name
        result.init()
        return result
    }

    /**
     * Find out the 'executor' to be used in the process definition or in teh session configuration object
     *
     * @param taskConfig
     */
    @CompileDynamic
    protected String getExecutorName(ProcessConfig taskConfig, Session session) {
        // create the processor object
        def result = taskConfig.executor?.toString()

        if( !result ) {
            if( session.config.executor instanceof String ) {
                result = session.config.executor
            }
            else if( session.config.executor?.name instanceof String ) {
                result = session.config.executor.name
            }
        }

        log.debug "<< taskConfig executor: $result"
        return result
    }

    void signalExecutors() {
        for( Executor exec : executors.values() )
            exec.signal()
    }

}
