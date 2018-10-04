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
