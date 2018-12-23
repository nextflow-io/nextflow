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

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.daemon.IgGridFactory
import nextflow.processor.TaskContext
import nextflow.processor.TaskRun
import nextflow.util.InputStreamDeserializer
import nextflow.util.RemoteSession
import org.apache.commons.lang.SerializationUtils
import org.apache.ignite.IgniteCache
import org.apache.ignite.IgniteException
/**
 * Execute a groovy closure task in a remote Ignite node
 */
@Slf4j
@CompileStatic
class IgClosureTask extends IgBaseTask<IgResultData> {

    private static final long serialVersionUID = 5515528753549263068L

    static final Map<UUID,GroovyClassLoader> classLoaderCache = new HashMap()


    /**
     * The task closure serialized as a byte array
     */
    final byte[] codeObj

    /**
     * The task delegate context serialized as a byte array
     */
    final byte[] delegateObj

    private transient StagingStrategy stagingStrategy


    IgClosureTask( TaskRun task, UUID sessionId ) {
        super(task,sessionId)
        this.codeObj = SerializationUtils.serialize(task.code.dehydrate())
        this.delegateObj = task.context.dehydrate()
    }

    void beforeExecute() {
        stagingStrategy = new IgFileStagingStrategy( task: bean, sessionId: sessionId )
        stagingStrategy.stage()
    }

    void afterExecute() {
        stagingStrategy.unstage()
    }

    @Override
    protected IgResultData execute0() throws IgniteException {
        log.debug "Running closure for task > ${bean.name}"

        def loader = getClassLoaderFor(sessionId)
        def delegate = TaskContext.rehydrate(delegateObj,loader)
        Closure closure = (Closure)InputStreamDeserializer.deserialize(codeObj,loader)
        Object result = closure.rehydrate(delegate, delegate.getScript(), delegate.getScript()).call()
        return new IgResultData(value: result, context: delegate?.getHolder())

    }

    @Override
    void cancel() {

    }


    /**
     * Create a {@link ClassLoader} object for the specified session ID
     *
     * @param sessionId
     * @param grid
     * @return
     */
    static protected ClassLoader getClassLoaderFor( UUID sessionId ) {
        assert sessionId

        (ClassLoader)classLoaderCache.getOrCreate(sessionId) {

            IgniteCache<UUID, RemoteSession> allSessions = IgGridFactory.ignite().cache( IgGridFactory.SESSIONS_CACHE )
            if( !allSessions )
                throw new IllegalStateException('Missing session cache object')

            final session = allSessions.get(sessionId)
            if( !session )
                throw new IllegalStateException("Missing session object for id: $sessionId")

            final result = new GroovyClassLoader()
            session.classpath.each { Path file ->
                log.debug "Adding to classpath: $file"
                result.addClasspath(file.toAbsolutePath().toString())
            }

            return result
        }
    }

}

