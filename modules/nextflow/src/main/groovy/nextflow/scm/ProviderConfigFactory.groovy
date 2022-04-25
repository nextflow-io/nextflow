package nextflow.scm

import org.pf4j.ExtensionPoint


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class ProviderConfigFactory implements ExtensionPoint {

    private static final List<String> DEFAULT_SCMS = [
            'github',
            'gitlab',
            'gitea',
            'bitbucket',
            'azurerepos']

    protected static final ProviderConfigFactory INSTANCE = new ProviderConfigFactory()

    List<ProviderConfig> allProviderConfigs() {
        DEFAULT_SCMS.collect{
            new ProviderConfig(it)
        }
    }

}
