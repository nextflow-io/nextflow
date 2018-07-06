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

import java.util.concurrent.Future
import java.util.concurrent.ThreadPoolExecutor

import groovy.util.logging.Slf4j
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Implements a proxy class that forwards method call invocations to
 * a thread pool execution which throttle requests according a specified rate limit
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class ClientProxyThrottler<T> implements GroovyInterceptable {

    /* the underlying thread pool executor */
    protected ThreadPoolExecutor executor

    protected T client

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
        opts.name = "${client.getClass().getSimpleName()}-ProxyThrottler"
        this.executor = ThrottlingExecutor.create(opts)
    }

    ClientProxyThrottler(T client, ThrottlingExecutor executor) {
        assert client
        assert executor
        this.client = client
        this.executor = executor
    }


    @Override
    Object invokeMethod(String name, Object args) {
        def dispatcher = name=='async' ? 'doInvoke0' : 'doInvoke1'
        return InvokerHelper.invokeMethod(executor, dispatcher, [client, name, args] as Object[])
    }

    Future async(Closure c) {
        // this is mock method
        // the real implementation is given by`invokeMethod`
        return null
    }


}
