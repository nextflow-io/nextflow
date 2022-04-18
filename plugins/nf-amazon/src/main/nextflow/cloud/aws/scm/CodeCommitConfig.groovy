package nextflow.cloud.aws.scm

import nextflow.scm.config.ScmConfig
import org.pf4j.Extension


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Extension
class CodeCommitConfig implements ScmConfig{

    @Override
    String getName() {
        'codecommit'
    }

    @Override
    Map enrichConfiguration(Map attr) {
        attr.platform = name
        attr.server = "https://git-codecommit.[a-z1-9-]+.amazonaws.com/v1"
        attr
    }
}
