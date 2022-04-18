package nextflow.scm

import nextflow.plugin.DefaultPlugins
import nextflow.plugin.PluginSpec
import nextflow.plugin.Plugins
import nextflow.plugin.PluginsFacade
import nextflow.scm.config.ScmConfig
import spock.lang.Specification

/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class CodeCommitProviderConfigTest extends Specification{

    static final String CONFIG = '''
        providers {
              codecommit {
                user = '12732'
                password = '35454'
              }
        }
        '''

    def 'should create a ProviderConfig object' () {

        when:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj1 = new ProviderConfig('codecommit', config.providers.codecommit as ConfigObject)

        then:
        obj1.name == 'codecommit'
        obj1.server == "https://git-codecommit.[a-z1-9-]+.amazonaws.com/v1"
        obj1.auth == '12732:35454'
        obj1.domain == 'git-codecommit.[a-z1-9-]+.amazonaws.com'
        obj1.platform == 'codecommit'
        obj1.endpoint == 'https://git-codecommit.[a-z1-9-]+.amazonaws.com/v1'
    }

}
