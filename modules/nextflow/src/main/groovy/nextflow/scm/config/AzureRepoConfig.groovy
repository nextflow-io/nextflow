package nextflow.scm.config

import org.pf4j.Extension


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Extension
class AzureRepoConfig implements ScmConfig{

    @Override
    String getName() {
        'azurerepos'
    }

    @Override
    Map enrichConfiguration(Map attr) {
        attr.platform = name
        if( !attr.server ) attr.server = 'https://dev.azure.com'
        if( !attr.endpoint ) attr.endpoint = 'https://dev.azure.com'
        attr
    }
}
