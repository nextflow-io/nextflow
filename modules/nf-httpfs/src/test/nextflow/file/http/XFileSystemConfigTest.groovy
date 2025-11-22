package nextflow.file.http

import spock.lang.Specification

import nextflow.SysEnv

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class XFileSystemConfigTest extends Specification {

    def 'should create with default config settings' () {
        when:
        def config = new XFileSystemConfig()
        then:
        config.retryCodes() == XFileSystemConfig.DEFAULT_RETRY_CODES.tokenize(',').collect( it -> it as int )
        config.backOffDelay() == XFileSystemConfig.DEFAULT_BACK_OFF_DELAY
        config.backOffBase() == XFileSystemConfig.DEFAULT_BACK_OFF_BASE
        config.maxAttempts() == XFileSystemConfig.DEFAULT_MAX_ATTEMPTS
    }

    def 'should create with custom config settings' () {
        given:
        SysEnv.push([NXF_HTTPFS_MAX_ATTEMPTS: '10',
                     NXF_HTTPFS_BACKOFF_BASE: '300',
                     NXF_HTTPFS_DELAY       : '400',
                     NXF_HTTPFS_RETRY_CODES : '1,2,3'])

        when:
        def config = new XFileSystemConfig()
        then:
        config.retryCodes() == [1,2,3]
        config.backOffDelay() == 400
        config.backOffBase() == 300
        config.maxAttempts() == 10

        cleanup:
        SysEnv.pop()
    }
}
