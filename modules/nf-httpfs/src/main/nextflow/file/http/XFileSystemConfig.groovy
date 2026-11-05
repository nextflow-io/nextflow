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
import groovy.transform.PackageScope
import nextflow.SysEnv

/**
 * Hold HTTP/FTP virtual file system configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class XFileSystemConfig {

    private static volatile XFileSystemConfig instance = new XFileSystemConfig()

    static final public String  DEFAULT_RETRY_CODES = '404,410,429,500,502,503,504'

    static final public int MAX_REDIRECT_HOPS = 5

    static final public int DEFAULT_BACK_OFF_BASE = 3

    static final public int DEFAULT_BACK_OFF_DELAY = 250

    static final public int DEFAULT_MAX_ATTEMPTS = 3

    private int maxAttempts = DEFAULT_MAX_ATTEMPTS

    private int backOffBase = DEFAULT_BACK_OFF_BASE

    private int backOffDelay = DEFAULT_BACK_OFF_DELAY

    private List<Integer> retryCodes

    XFileSystemConfig() {
        this(SysEnv.get())
    }

    XFileSystemConfig(Map<String,String> env) {
        maxAttempts = config(env, 'NXF_HTTPFS_MAX_ATTEMPTS', DEFAULT_MAX_ATTEMPTS) as Integer
        backOffBase = config(env, 'NXF_HTTPFS_BACKOFF_BASE', DEFAULT_BACK_OFF_BASE) as Integer
        backOffDelay = config(env, 'NXF_HTTPFS_DELAY', DEFAULT_BACK_OFF_DELAY) as Integer
        retryCodes = config(env, 'NXF_HTTPFS_RETRY_CODES', DEFAULT_RETRY_CODES).tokenize(',').collect( val -> val as Integer )
    }

    private static String config(Map<String,String> env, String name, def defValue) {
        return env.containsKey(name) ? env.get(name) : defValue.toString()
    }

    static String config(String name, def defValue) {
        return config(SysEnv.get(), name, defValue)
    }

    int maxAttempts() { maxAttempts }

    int backOffBase() { backOffBase }

    int backOffDelay() { backOffDelay }

    List<Integer> retryCodes() { retryCodes }

    static XFileSystemConfig config() { return instance }

    /**
     * Replace the global config used by all providers, or reload it
     * from the environment when no config is given. Intended for tests.
     */
    @PackageScope
    static void reset(XFileSystemConfig config=null) {
        instance = config ?: new XFileSystemConfig()
    }
}
