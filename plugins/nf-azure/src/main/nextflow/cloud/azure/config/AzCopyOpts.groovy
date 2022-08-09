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

    //-----------------------------------------------------
    // Default values for azcopy copy command options
    //-----------------------------------------------------

    //Use this block size (specified in MiB) when uploading/downloading to/from Azure Storage (azcopy default: 8)
    //Can be in decimals for eg. 0.25 and maximum value is 100. https://github.com/Azure/azure-storage-azcopy/wiki/Cost-of-data-transfers
    static public final String DEFAULT_BLOCK_SIZE = "4"
    String blockSize

    //Upload block blob to Azure Storage using this blob tier. (azcopy default: "None")
    static public final String DEFAULT_BLOB_TIER = "None" // hot (None) | cool | archive
    String blobTier

    //Overwrite the conflicting files and blobs at the destination if this flag is set to true. (azcopy default: true)
    static public final String DEFAULT_OVERWRITE = "false" // true | false | prompt | ifSourceNewer
    String overwrite

    //The Azure Blob Storage service automatically computes MD5 sum for files less than 64 MB in size.
    //Content-MD5 property of the destination blob or file. (azcopy default: false)
    static public final Boolean DEFAULT_PUT_MD5 = false
    Boolean putMD5

    //Specifies how strictly MD5 hashes should be validated when downloading. (azcopy default: "FailIfDifferent")
    // NoCheck| LogOnly| FailIfDifferent| FailIfDifferentOrMissing
    static public final String DEFAULT_CHECK_MD5 = "FailIfDifferent"
    String checkMD5


    AzCopyOpts() {
        this.blockSize = DEFAULT_BLOCK_SIZE
        this.blobTier = DEFAULT_BLOB_TIER
        this.putMD5 = DEFAULT_PUT_MD5
        this.checkMD5 = DEFAULT_CHECK_MD5
        this.overwrite = DEFAULT_OVERWRITE
    }


    AzCopyOpts(Map config) {
        assert config != null

        this.blockSize = config.blockSize ?: DEFAULT_BLOCK_SIZE
        this.blobTier = config.blobTier ?: DEFAULT_BLOB_TIER
        this.putMD5 = config.putMD5 ?: DEFAULT_PUT_MD5
        this.checkMD5 = config.checkMD5 ?: DEFAULT_CHECK_MD5
        this.overwrite = config.overwrite ?: DEFAULT_OVERWRITE

    }

}
