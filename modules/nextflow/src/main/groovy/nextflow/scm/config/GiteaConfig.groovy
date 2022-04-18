package nextflow.scm.config

import org.pf4j.Extension


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Extension
class GiteaConfig implements ScmConfig{

    @Override
    String getName() {
        'gitea'
    }

    @Override
    Map enrichConfiguration(Map attr) {
        attr.platform = name
        if( !attr.server ) attr.server = 'https://try.gitea.io'
        if( !attr.endpoint ) attr.endpoint = attr.server.toString().stripEnd('/') + '/api/v1'
        attr
    }
}
