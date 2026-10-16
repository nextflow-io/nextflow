/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.file.http

import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.util.Duration
import nextflow.util.RetryConfig

/**
 * Hold HTTP/FTP virtual file system configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class XFileSystemConfig {

    static final public String  DEFAULT_RETRY_CODES = '404,410,429,500,502,503,504'

    static final public int MAX_REDIRECT_HOPS = 5

    static final public int DEFAULT_BACK_OFF_BASE = 3

    static final public int DEFAULT_BACK_OFF_DELAY = 250

    static final public int DEFAULT_MAX_ATTEMPTS = 3

    static final public Duration DEFAULT_CONNECT_TIMEOUT = Duration.of('60s')

    static final public Duration DEFAULT_READ_TIMEOUT = Duration.of('60s')

    /**
     * The timeouts are shared with the other Nextflow HTTP clients, see {@code nextflow.scm.RepositoryProvider}
     */
    static final private String HTTP_CLIENT_ENV_PREFIX = 'NXF_HTTPCLIENT_'

    // declared after the defaults above, which the constructor reads
    static XFileSystemConfig instance = new XFileSystemConfig()

    private int maxAttempts = DEFAULT_MAX_ATTEMPTS

    private int backOffBase = DEFAULT_BACK_OFF_BASE

    private int backOffDelay = DEFAULT_BACK_OFF_DELAY

    private List<Integer> retryCodes

    private Duration connectTimeout

    private Duration readTimeout

    {
        maxAttempts = config('NXF_HTTPFS_MAX_ATTEMPTS', DEFAULT_MAX_ATTEMPTS) as Integer
        backOffBase = config('NXF_HTTPFS_BACKOFF_BASE', DEFAULT_BACK_OFF_BASE) as Integer
        backOffDelay = config('NXF_HTTPFS_DELAY', DEFAULT_BACK_OFF_DELAY) as Integer
        retryCodes = config('NXF_HTTPFS_RETRY_CODES', DEFAULT_RETRY_CODES).tokenize(',').collect( val -> val as Integer )
        connectTimeout = RetryConfig.valueOf(null, 'connectTimeout', HTTP_CLIENT_ENV_PREFIX, DEFAULT_CONNECT_TIMEOUT, Duration)
        readTimeout = RetryConfig.valueOf(null, 'readTimeout', HTTP_CLIENT_ENV_PREFIX, DEFAULT_READ_TIMEOUT, Duration)
        if( connectTimeout.toMillis() <= 0 )
            throw new IllegalArgumentException("Variable ${HTTP_CLIENT_ENV_PREFIX}CONNECT_TIMEOUT must be greater than zero - offending value: ${connectTimeout}")
    }

    static String config(String name, def defValue) {
        return SysEnv.containsKey(name) ? SysEnv.get(name) : defValue.toString()
    }

    int maxAttempts() { maxAttempts }

    int backOffBase() { backOffBase }

    int backOffDelay() { backOffDelay }

    List<Integer> retryCodes() { retryCodes }

    /**
     * @return The max time to establish a connection, in milliseconds
     */
    int connectTimeoutMillis() { toMillis(connectTimeout) }

    /**
     * @return The max time to wait for data once connected, in milliseconds, or {@code 0} to wait indefinitely
     */
    int readTimeoutMillis() { toMillis(readTimeout) }

    static private int toMillis(Duration duration) {
        return (int) Math.min(duration.toMillis(), (long) Integer.MAX_VALUE)
    }

    static XFileSystemConfig config() { return instance }
}
