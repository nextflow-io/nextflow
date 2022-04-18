package nextflow.scm.config

import org.pf4j.Extension


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Extension
class GitlabConfig implements ScmConfig{

    @Override
    String getName() {
        'gitlab'
    }

    @Override
    Map enrichConfiguration(Map attr) {
        attr.platform = name
        if( !attr.server ) attr.server = 'https://gitlab.com'
        attr
    }
}
