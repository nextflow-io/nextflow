package nextflow.container


import spock.lang.Specification

/**
 * Apptainer config tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class SingularityConfigTest extends Specification {

    def 'should set default values empty map'() {
        given:
        def config = new SingularityConfig([:])

        expect:
        config.envWhitelist == []
        config.pullTimeout.toMillis() == 20 * 60 * 1000 //20 min
    }

    def 'should create config with full map'(){
        given:
        def configMap = [
            autoMounts: false,
            cacheDir: 'cacheDir',
            enabled: true,
            engineOptions: '-q -v',
            envWhitelist: 'ENV_1,ENV_2',
            libraryDir: 'libraryDir',
            noHttps: false,
            ociAutoPull: false,
            pullTimeout: '50s',
            registry: 'http://registry.com',
            runOptions: '--contain --writable'
        ]
        def config = new SingularityConfig(configMap)

        expect:
        config.autoMounts == false
        config.cacheDir == 'cacheDir'
        config.enabled == true
        config.engineOptions == '-q -v'
        config.envWhitelist == ['ENV_1','ENV_2']
        config.libraryDir == 'libraryDir'
        config.noHttps == false
        config.ociAutoPull == false
        config.registry == 'http://registry.com'
        config.runOptions == '--contain --writable'
        config.pullTimeout.toMillis() == 50_000 // 50s

    }
}
