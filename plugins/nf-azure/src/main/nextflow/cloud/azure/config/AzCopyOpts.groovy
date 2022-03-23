/*
 * Copyright 2021, Microsoft Corp
 * Copyright 2022, Seqera Labs
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.cloud.azure.config

import groovy.transform.CompileStatic

/**
 * Model Azure azcopy tool config settings from nextflow config file
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@CompileStatic
class AzCopyOpts {

    static public final String DEFAULT_BLOCK_SIZE = "4"
    static public final String DEFAULT_BLOB_TIER = "None"

    String blockSize
    String blobTier

    AzCopyOpts() {
        this.blockSize = DEFAULT_BLOCK_SIZE
        this.blobTier =  DEFAULT_BLOB_TIER
    }


    AzCopyOpts(Map config) {
        assert config!=null
        this.blockSize = config.blockSize ?: DEFAULT_BLOCK_SIZE
        this.blobTier = config.blobTier ?: DEFAULT_BLOB_TIER
    }

}
