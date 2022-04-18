package nextflow.scm.config

import org.pf4j.Extension


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Extension
class BitBucketConfig implements ScmConfig{

    @Override
    String getName() {
        'bitbucket'
    }

    @Override
    Map enrichConfiguration(Map attr) {
        attr.platform = name
        if( !attr.server ) attr.server = 'https://bitbucket.org'
        attr
    }
}
