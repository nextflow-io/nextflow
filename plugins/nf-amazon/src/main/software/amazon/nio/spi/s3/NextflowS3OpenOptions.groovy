package software.amazon.nio.spi.s3

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.util.AwsHelper
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration
import software.amazon.awssdk.services.s3.model.GetObjectRequest
import software.amazon.awssdk.services.s3.model.ObjectCannedACL
import software.amazon.awssdk.services.s3.model.PutObjectRequest
import software.amazon.awssdk.services.s3.model.RequestPayer
import software.amazon.awssdk.services.s3.model.ServerSideEncryption

import java.nio.file.Path

@Slf4j
@CompileStatic
class NextflowS3OpenOptions extends S3OpenOption {
    private Boolean isRequesterPays
    private ObjectCannedACL cannedACL
    private String kmsKeyId
    private ServerSideEncryption storageEncryption
    private ClientOverrideConfiguration clientOverride

    NextflowS3OpenOptions(){}

    protected NextflowS3OpenOptions(Boolean isRequesterPays, ObjectCannedACL cannedACL, String kmsKeyId,
                                    ServerSideEncryption storageEncryption, ClientOverrideConfiguration clientOverride) {
        this.isRequesterPays = isRequesterPays
        this.cannedACL = cannedACL
        this.kmsKeyId = kmsKeyId
        this.storageEncryption = storageEncryption
        this.clientOverride = clientOverride
    }

    @Override
    S3OpenOption copy() {
        return new NextflowS3OpenOptions(this.isRequesterPays, this.cannedACL, this.kmsKeyId, this.storageEncryption, this.clientOverride)
    }

    @Override
    protected void apply(GetObjectRequest.Builder getObjectRequest) {
        if( isRequesterPays )
            getObjectRequest.requestPayer(RequestPayer.REQUESTER)
    }

    @Override
    protected void apply(PutObjectRequest.Builder putObjectRequest, Path file) {
        if( cannedACL )
            putObjectRequest.acl(cannedACL)
        if( kmsKeyId )
            putObjectRequest.ssekmsKeyId(kmsKeyId)
        if( storageEncryption )
            putObjectRequest.serverSideEncryption(storageEncryption)
        if( isRequesterPays )
            putObjectRequest.requestPayer(RequestPayer.REQUESTER)

    }

    void setCannedAcl(String acl) {
        if( acl == null )
            return
        this.cannedAcl = AwsHelper.parseS3Acl(acl)
    }

    void setKmsKeyId(String kmsKeyId) {
        if( kmsKeyId == null )
            return
        this.kmsKeyId = kmsKeyId
    }

    void setStorageEncryption(String alg) {
        if( alg == null )
            return;
        this.storageEncryption = ServerSideEncryption.fromValue(alg)
    }

    void setRequesterPays(String requesterPaysEnabled) {
        if( requesterPaysEnabled == null )
            return;
        this.isRequesterPays = Boolean.valueOf(requesterPaysEnabled)
    }

    void setClientOverride(ClientOverrideConfiguration clientOverride){
        if( clientOverride == null)
            return;
        this.clientOverride = clientOverride
    }
}
