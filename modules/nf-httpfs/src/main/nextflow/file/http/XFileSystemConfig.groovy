package nextflow.file.http

import groovy.transform.CompileStatic
import nextflow.SysEnv

/**
 * Hold HTTP/FTP virtual file system configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class XFileSystemConfig {

    static XFileSystemConfig instance = new XFileSystemConfig()

    static final public String  DEFAULT_RETRY_CODES = '404,410'

    static final public int MAX_REDIRECT_HOPS = 5

    static final public int DEFAULT_BACK_OFF_BASE = 3

    static final public int DEFAULT_BACK_OFF_DELAY = 250

    static final public int DEFAULT_MAX_ATTEMPTS = 3

    private int maxAttempts = DEFAULT_MAX_ATTEMPTS

    private int backOffBase = DEFAULT_BACK_OFF_BASE

    private int backOffDelay = DEFAULT_BACK_OFF_DELAY

    private List<Integer> retryCodes

    {
        maxAttempts = config('NXF_HTTPFS_MAX_ATTEMPTS', DEFAULT_MAX_ATTEMPTS) as Integer
        backOffBase = config('NXF_HTTPFS_BACKOFF_BASE', DEFAULT_BACK_OFF_BASE) as Integer
        backOffDelay = config('NXF_HTTPFS_DELAY', DEFAULT_BACK_OFF_DELAY) as Integer
        retryCodes = config('NXF_HTTPFS_RETRY_CODES', DEFAULT_RETRY_CODES).tokenize(',').collect( val -> val as Integer )
    }

    static String config(String name, def defValue) {
        return SysEnv.containsKey(name) ? SysEnv.get(name) : defValue.toString()
    }

    int maxAttempts() { maxAttempts }

    int backOffBase() { backOffBase }

    int backOffDelay() { backOffDelay }

    List<Integer> retryCodes() { retryCodes }

    static XFileSystemConfig config() { return instance }
}
