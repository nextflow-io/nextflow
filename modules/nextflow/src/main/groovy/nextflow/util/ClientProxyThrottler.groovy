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

package nextflow.util


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Implements a proxy class that forwards method call invocations to
 * a thread pool execution which throttle requests according a specified rate limit
 *
 * WARN: the caller class/method should not be compile static
 *
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class ClientProxyThrottler<T> implements GroovyInterceptable {

    /* the underlying thread pool executor */
    private ThrottlingExecutor executor

    /* the target client object that will manage the requests to be throttled */
    private T client

    /**
     * Allow the definition of an optional priority for each method name
     */
    private Map<String,Byte> priority = [:]

    /**
     * Create the proxy throttler object. Subclasses must provide the
     * target client instance as first parameter and the throttling options
     * as second parameter
     *
     * @param client The target client whose method invocation needs to be throttled
     * @param opts The throttling options
     */
    ClientProxyThrottler(T client, ThrottlingExecutor.Options opts) {
        assert client
        this.client = client
        opts.poolName = "${client.getClass().getSimpleName()}-ProxyThrottler"
        this.executor = ThrottlingExecutor.create(opts)
    }

    /**
     * Create the proxy throttler object using the specified {@link ThrottlingExecutor} instance
     *
     * @param client The target client whose method invocation needs to be throttled
     * @param executor The {@link ThrottlingExecutor} scheduling and executing the actual requests to the client
     * @param priority An optional Map that allows the specification of the execution priority for a given method name
     */
    ClientProxyThrottler(T client, ThrottlingExecutor executor, Map<String,Byte> priority = [:]) {
        assert client
        assert executor
        this.client = client
        this.executor = executor
        this.priority = priority
    }

    Object getProperty(String name) {
        name=='client' ? client : InvokerHelper.getProperty(this,name)
    }

    @Override
    Object invokeMethod(String name, Object args) {
        if( name=='getClient' && !args )
            return client

        Byte p = priority.get(name)
        if( p == null ) p = 0
        return InvokerHelper.invokeMethod(executor, 'doInvoke1', [client, name, args, p] as Object[])
    }

    /**
     * @return The client instance associate to this proxy. This allows the direct
     * access with the client bypassing the throttling logic
     */
    T getClient() { client }
}
