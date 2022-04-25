package nextflow.cloud.aws.codecommit

import nextflow.scm.ProviderConfig
import nextflow.scm.ProviderConfigFactory


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class CodeCommitProviderConfigFactory extends ProviderConfigFactory{

    @Override
    List<ProviderConfig> allProviderConfigs() {
        [new ProviderConfig('codecommit')]
    }
}
