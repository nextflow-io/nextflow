package nextflow.cloud.azure.batch

import java.nio.file.Path
import java.time.OffsetDateTime

import com.azure.storage.blob.BlobClient
import com.azure.storage.blob.nio.AzurePath
import com.azure.storage.blob.sas.BlobSasPermission
import com.azure.storage.blob.sas.BlobServiceSasSignatureValues
import groovy.transform.CompileStatic
import nextflow.util.Duration
/**
 * Azure helper functions
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzHelper {

    static private AzurePath az0(Path path){
        if( path !instanceof AzurePath )
            throw new IllegalArgumentException("Not a valid Azure path: $path [${path?.getClass()?.getName()}]")
        return (AzurePath)path
    }

    static String toHttpUrl(Path path, String sas=null) {
        final client = az0(path).toBlobClient()
        return !sas ? client.getBlobUrl() : "${client.getBlobUrl()}?${sas}"
    }

    static String generateSas(Path path, Duration duration, String perms) {
        generateSas(az0(path).toBlobClient(), duration, perms)
    }

    static String generateSas(BlobClient client, Duration duration, String perms) {
        final offset = OffsetDateTime
                .now()
                .plusSeconds(duration.seconds)

        return client
                .generateSas(new BlobServiceSasSignatureValues(offset, BlobSasPermission.parse(perms)))
    }

}
