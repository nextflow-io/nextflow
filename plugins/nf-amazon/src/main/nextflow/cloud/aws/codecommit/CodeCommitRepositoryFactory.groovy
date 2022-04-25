package nextflow.cloud.aws.codecommit


import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryFactory
import nextflow.scm.RepositoryProvider
import org.pf4j.ExtensionPoint


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class CodeCommitRepositoryFactory extends RepositoryFactory{

    @Override
    RepositoryProvider newInstance(ProviderConfig config, String project, String url) {
        switch(config.platform) {
            case 'codecommit':
                try {
                    return new AwsCodeCommitRepositoryProvider(project, url, config)
                }catch(e){
                    e.printStackTrace()
                    throw e
                }
            default:
                return null
        }
    }
}
