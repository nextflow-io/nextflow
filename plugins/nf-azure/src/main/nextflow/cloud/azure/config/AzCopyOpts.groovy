package nextflow.cloud.azure.config

/**
 * Model Azure azcopy tool config settings from nextflow config file
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */

import groovy.transform.CompileStatic


@CompileStatic
class AzCopyOpts {

    String blockSize
    String blobTier

    AzCopyOpts() {
        this(Collections.emptyMap())
    }

    AzCopyOpts(Map config, Map<String,String> env=null) {
        assert config!=null
        this.blockSize = config.blockSize
        this.blobTier = config.blobTier
    }

}
