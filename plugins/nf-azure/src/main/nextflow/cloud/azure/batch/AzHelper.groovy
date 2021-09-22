/*
 * Copyright 2021, Microsoft Corp
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
package nextflow.cloud.azure.batch

import java.nio.file.Path
import java.time.OffsetDateTime

import com.azure.storage.blob.BlobContainerClient
import com.azure.storage.blob.sas.BlobContainerSasPermission
import com.azure.storage.blob.sas.BlobSasPermission
import com.azure.storage.blob.sas.BlobServiceSasSignatureValues
import groovy.transform.CompileStatic
import nextflow.cloud.azure.nio.AzPath
import nextflow.util.Duration
/**
 * Azure helper functions
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzHelper {

    static private AzPath az0(Path path){
        if( path !instanceof AzPath )
            throw new IllegalArgumentException("Not a valid Azure path: $path [${path?.getClass()?.getName()}]")
        return (AzPath)path
    }

    static String toHttpUrl(Path path, String sas=null) {
        def url = az0(path).blobClient().getBlobUrl()
        url = URLDecoder.decode(url, 'UTF-8').stripEnd('/')
        return !sas ? url : "${url}?${sas}"
    }

    static String toContainerUrl(Path path, String sas) {
        def url = az0(path).containerClient().getBlobContainerUrl()
        url = URLDecoder.decode(url, 'UTF-8').stripEnd('/')
        return !sas ? url : "${url}?${sas}"
    }

    static String generateContainerSas(Path path, Duration duration) {
        generateSas(az0(path).containerClient(), duration)
    }

    static BlobContainerSasPermission CONTAINER_PERMS = new BlobContainerSasPermission()
            .setAddPermission(true)
            .setCreatePermission(true)
            .setDeletePermission(true)
            .setListPermission(true)
            .setMovePermission(true)
            .setReadPermission(true)
            .setTagsPermission(true)
            .setWritePermission(true)

    static BlobSasPermission BLOB_PERMS = new BlobSasPermission()
            .setAddPermission(true)
            .setCreatePermission(true)
            .setDeletePermission(true)
            .setListPermission(true)
            .setMovePermission(true)
            .setReadPermission(true)
            .setTagsPermission(true)
            .setWritePermission(true)


    static String generateSas(BlobContainerClient client, Duration duration) {
        final now = OffsetDateTime .now()

        final signature = new BlobServiceSasSignatureValues()
                .setPermissions(BLOB_PERMS)
                .setPermissions(CONTAINER_PERMS)
                .setStartTime(now)
                .setExpiryTime( now.plusSeconds(duration.toSeconds()) )

        return client .generateSas(signature)
    }

}
