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
    // Default values for azcopy CLI environment variables
    //-----------------------------------------------------

    //Description: Overrides where the log files are stored, to avoid filling up a disk.
    //AZCOPY_LOG_LOCATION
    // String logLocation

    //Description: Overrides how many HTTP connections work on transfers. By default, this number is determined based on the number of logical cores on the machine.
    //AZCOPY_CONCURRENCY_VALUE
    // String concurrencyValue

    //Description: Overrides the (approximate) number of files that are in progress at any one time, by controlling how many files we concurrently initiate transfers for.
    //AZCOPY_CONCURRENT_FILES
    // String concurrentFiles

    //Description: Max number of GB that AzCopy should use for buffering data between network and disk. May include decimal point, e.g. 0.5. The default is based on machine size.
    //AZCOPY_BUFFER_GB
    // String bufferGB

    //Description: Set time (in minutes) for how long AzCopy should try to upload files for each request before AzCopy times out.
    //AZCOPY_REQUEST_TRY_TIMEOUT
    // static public final String AZCOPY_REQUEST_TRY_TIMEOUT = "4"
    // String requestTryTimeout


    //-----------------------------------------------------
    // Default values for azcopy copy command options
    //-----------------------------------------------------

    //Use this block size (specified in MiB) when uploading/downloading to/from Azure Storage. Can be in decimals for eg. 0.25
    static public final String DEFAULT_BLOCK_SIZE = "4"
    String blockSize

    //Upload block blob to Azure Storage using this blob tier. (azcopy default: "None")
    static public final String DEFAULT_BLOB_TIER = "None"
    String blobTier

    // static public final String DEFAULT_METADATA = ""
    // String metadata

    // static public final String DEFAULT_BLOB_TAGS = ""
    // String blobTags

    //Define the log verbosity for the log file (azcopy default: "INFO").
    //We set the default to NONE to save space on the device
    // static public final String DEFAULT_LOG_LEVEL = "NONE"
    // String logLevel

    //Define the output verbosity (azcopy default: "default").
    //We set the default
    static public final String DEFAULT_OUTPUT_LEVEL = "quiet"
    String outputLevel

    //The azcopy default is true, which means upon `-resume` the data is uploaded again.
    //Overwrite the conflicting files and blobs at the destination if this flag is set to true. (azcopy default: true)
    static public final String DEFAULT_OVERWRITE = "false"
    String overwrite

    //The azcopy default is true, which means upon `-resume` the data is uploaded again.
    //static public final Boolean DEFAULT_RECURSIVE = false
    //String recursive

    //Check the length of a file on the destination after the transfer.
    //If there is a mismatch between source and destination, the transfer is marked as failed (azcopy default: true)
    static public final Boolean DEFAULT_CHECK_LENGTH = true
    Boolean checkLength

    //The Azure Blob Storage service automatically computes MD5 sum for files less than 256 MB in size.
    //Content-MD5 property of the destination blob or file. (azcopy default: false)
    static public final Boolean DEFAULT_PUT_MD5 = false
    Boolean putMD5

    //Specifies how strictly MD5 hashes should be validated when downloading. (azcopy default: "FailIfDifferent")
    static public final String DEFAULT_CHECK_MD5 = "FailIfDifferent"
    String checkMD5


    AzCopyOpts() {
        this.blockSize = DEFAULT_BLOCK_SIZE
        this.blobTier = DEFAULT_BLOB_TIER
        this.putMD5 = DEFAULT_PUT_MD5
        this.checkMD5 = DEFAULT_CHECK_MD5
        this.overwrite = DEFAULT_OVERWRITE
        this.outputLevel = DEFAULT_OUTPUT_LEVEL
    }


    AzCopyOpts(Map config) {
        assert config != null

        this.blockSize = config.blockSize ?: DEFAULT_BLOCK_SIZE
        this.blobTier = config.blobTier ?: DEFAULT_BLOB_TIER
        this.putMD5 = config.putMD5 ?: DEFAULT_PUT_MD5
        this.checkMD5 = config.checkMD5 ?: DEFAULT_CHECK_MD5
        this.overwrite = config.overwrite ?: DEFAULT_OVERWRITE
        this.outputLevel = config.outputLevel ?: DEFAULT_OUTPUT_LEVEL

    }

}
