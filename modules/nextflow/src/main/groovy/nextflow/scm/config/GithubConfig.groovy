package nextflow.scm.config

import org.pf4j.Extension


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Extension
class GithubConfig implements ScmConfig{

    @Override
    String getName() {
        'github'
    }

    @Override
    Map enrichConfiguration(Map attr) {
        attr.platform = name
        if( !attr.server ) attr.server = 'https://github.com'
        if( !attr.endpoint ) attr.endpoint = 'https://api.github.com'
        attr
    }
}
