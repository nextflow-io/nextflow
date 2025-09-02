/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package nextflow.cloud.aws.nio;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.*;
import java.util.function.Consumer;
import java.util.function.Supplier;

import nextflow.cloud.aws.nio.util.S3SyncClientConfiguration;
import nextflow.cloud.aws.nio.util.S3TransferManagerSync;
import nextflow.cloud.aws.AwsClientFactory;
import nextflow.cloud.aws.nio.util.S3AsyncClientConfiguration;
import nextflow.cloud.aws.util.AwsHelper;
import nextflow.util.ThreadPoolManager;
import nextflow.util.Threads;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.core.ResponseInputStream;
import software.amazon.awssdk.core.exception.SdkException;
import software.amazon.awssdk.core.sync.RequestBody;
import software.amazon.awssdk.services.s3.model.*;
import software.amazon.awssdk.services.s3.paginators.ListObjectsV2Iterable;
import software.amazon.awssdk.transfer.s3.S3TransferManager;
import software.amazon.awssdk.transfer.s3.model.*;

/**
 * Client Amazon S3
 * @see software.amazon.awssdk.services.s3.S3Client
 */
public class S3Client {

	private static final Logger log = LoggerFactory.getLogger(S3Client.class);

	private software.amazon.awssdk.services.s3.S3Client client;

	// Semaphore to limit concurrent client connections when using virtual threads.
	private Semaphore semaphore;

	private ObjectCannedACL cannedAcl;

	private String kmsKeyId;

	private ServerSideEncryption storageEncryption;

	private S3TransferManagerSync transferManager;

	private ExecutorService transferPool;

	private Boolean isRequesterPaysEnabled = false;

	private String callerAccount;

	private AwsClientFactory factory;

	private Properties props;

	private boolean global;

	public S3Client(AwsClientFactory factory, Properties props, boolean global) {
		S3SyncClientConfiguration clientConfig = S3SyncClientConfiguration.create(props);
		this.factory = factory;
		this.props = props;
		this.global = global;
		this.client = factory.getS3Client(clientConfig, global);
		this.semaphore = Threads.useVirtual() ? new Semaphore(clientConfig.getMaxConnections()) : null;
		this.callerAccount = fetchCallerAccount();
	}

	/**
	 * Perform an action that requires the S3 semaphore to limit concurrent connections.
	 *
	 * @param action
	 */
	private <T> T runWithPermit(Supplier<T> action) {
		try {
			if (semaphore != null) semaphore.acquire();
			try {
				return action.get();
			} finally {
				if (semaphore != null) semaphore.release();
			}
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new RuntimeException("Interrupted while acquiring S3 client semaphore", e);
		}
	}

	/**
	 * AmazonS3Client#getS3AccountOwner() is not available in SDK v2.
	 * The STSClient#getCallerIdentity returns the account, but it does not include the canonical ID required for ACLs.
	 *
	 * This function and the fetchCallerAccount() emulate the old behavior retrieving the canonicalId can only be
	 * retrieved if the user owns a bucket.
	 */
	public String getCallerAccount() {
		return callerAccount;
	}

	private String fetchCallerAccount(){
		try {
			List<Bucket> buckets = runWithPermit(() -> client.listBuckets(ListBucketsRequest.builder().maxBuckets(1).build()).buckets());
			if (buckets == null || buckets.isEmpty())
				return null;
			return getBucketAcl(buckets.get(0).name()).owner().id();
		}catch (Throwable e){
			log.debug("Unable to fetch caller account - {} ", e.getMessage());
			return null;
		}
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#listBuckets()
	 */
	public List<Bucket> listBuckets() {
		return runWithPermit(() -> client.listBuckets().buckets());
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#listObjects(ListObjectsRequest)
	 */
	public ListObjectsResponse listObjects(ListObjectsRequest request) {
		return runWithPermit(() -> client.listObjects(request));
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#getObject
	 */
	public ResponseInputStream<GetObjectResponse> getObject(String bucketName, String key) {
		GetObjectRequest.Builder reqBuilder = GetObjectRequest.builder().bucket(bucketName).key(key);
		if( this.isRequesterPaysEnabled )
			reqBuilder.requestPayer(RequestPayer.REQUESTER);
		return runWithPermit(() -> client.getObject(reqBuilder.build()));
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#putObject
	 */
	public PutObjectResponse putObject(String bucket, String key, File file) {
		PutObjectRequest.Builder builder = PutObjectRequest.builder().bucket(bucket).key(key);
		if( cannedAcl != null ) {
			log.trace("Setting canned ACL={}; bucket={}; key={}", cannedAcl, bucket, key);
			builder.acl(cannedAcl);
		}
		return runWithPermit(() -> client.putObject(builder.build(), file.toPath()));
	}

	private PutObjectRequest preparePutObjectRequest(PutObjectRequest.Builder reqBuilder, List<Tag> tags, String contentType, String storageClass) {
		if( cannedAcl != null ) {
			reqBuilder.acl(cannedAcl);
		}
		if( tags != null && tags.size()>0 ) {
			reqBuilder.tagging(Tagging.builder().tagSet(tags).build());
		}
		if( kmsKeyId != null ) {
			reqBuilder.ssekmsKeyId(kmsKeyId);
		}
		if( storageEncryption!=null ) {
			reqBuilder.serverSideEncryption(storageEncryption);
		}
		if( contentType!=null ) {
			reqBuilder.contentType(contentType);
		}
		if( storageClass!=null ) {
			reqBuilder.storageClass(storageClass);
		}
		return reqBuilder.build();
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#putObject
	 */
	public PutObjectResponse putObject(String bucket, String keyName, InputStream inputStream, List<Tag> tags, String contentType, long contentLength) {
		PutObjectRequest.Builder reqBuilder = PutObjectRequest.builder()
			.bucket(bucket)
			.key(keyName);
		if( cannedAcl != null ) {
			reqBuilder.acl(cannedAcl);
		}
		if( tags != null && tags.size()>0 ) {
			reqBuilder.tagging(Tagging.builder().tagSet(tags).build());
		}
		if( kmsKeyId != null ) {
			reqBuilder.ssekmsKeyId(kmsKeyId);
		}
		if( storageEncryption!=null ) {
			reqBuilder.serverSideEncryption(storageEncryption);
		}
		if( contentType!=null ) {
			reqBuilder.contentType(contentType);
		}
		PutObjectRequest req = reqBuilder.build();
		if( log.isTraceEnabled() ) {
			log.trace("S3 PutObject request {}", req);
		}
		return runWithPermit(() -> client.putObject(req, RequestBody.fromInputStream(inputStream, contentLength)));
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#deleteObject
	 */
	public void deleteObject(String bucket, String key) {
		client.deleteObject(DeleteObjectRequest.builder().bucket(bucket).key(key).build());
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#copyObject(CopyObjectRequest)
	 */
	public void copyObject(CopyObjectRequest.Builder reqBuilder, List<Tag> tags, String contentType, String storageClass) {
		if( tags !=null && !tags.isEmpty()) {
			log.debug("Setting tags: {}", tags);
			reqBuilder.taggingDirective(TaggingDirective.REPLACE);
			reqBuilder.tagging(Tagging.builder().tagSet(tags).build());
		}
		if( cannedAcl != null ) {
			reqBuilder.acl(cannedAcl);
		}
		if( storageEncryption != null ) {
			reqBuilder.serverSideEncryption(storageEncryption);
		}
		if( kmsKeyId !=null ) {
			reqBuilder.ssekmsKeyId(kmsKeyId);
		}
		if( contentType!=null ) {
			reqBuilder.metadataDirective(MetadataDirective.REPLACE);
			reqBuilder.contentType(contentType);
		}
		if( storageClass!=null ) {
			reqBuilder.storageClass(storageClass);
		}
		CopyObjectRequest req = reqBuilder.build();
		if( log.isTraceEnabled() ) {
			log.trace("S3 CopyObject request {}", req);
		}

		runWithPermit(() -> client.copyObject(req) );
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#getBucketAcl
	 */
	public AccessControlPolicy getBucketAcl(String bucket) {
		GetBucketAclResponse response = runWithPermit(() -> client.getBucketAcl(GetBucketAclRequest.builder().bucket(bucket).build()));
		return AccessControlPolicy.builder().grants(response.grants()).owner(response.owner()).build();
	}

	public void setCannedAcl(String acl) {
		if( acl==null )
			return;
		this.cannedAcl = AwsHelper.parseS3Acl(acl);
		log.debug("Setting S3 canned ACL={} [{}]", this.cannedAcl, acl);
	}

	public void setKmsKeyId(String kmsKeyId) {
		if( kmsKeyId==null )
			return;
		this.kmsKeyId = kmsKeyId;
		log.debug("Setting S3 SSE kms Id={}", kmsKeyId);
	}

	public void setStorageEncryption(String alg) {
		if( alg == null )
			return;
		this.storageEncryption = ServerSideEncryption.fromValue(alg);
		log.debug("Setting S3 SSE storage encryption algorithm={}", alg);
	}

	public void setRequesterPaysEnabled(String requesterPaysEnabled) {
		if( requesterPaysEnabled == null )
			return;
		this.isRequesterPaysEnabled = Boolean.valueOf(requesterPaysEnabled);
		log.debug("Setting S3 requester pays enabled={}", isRequesterPaysEnabled);
	}

	public ObjectCannedACL getCannedAcl() {
		return cannedAcl;
	}

	public software.amazon.awssdk.services.s3.S3Client getClient() {
		return client;
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#getObjectAcl
	 */
	public AccessControlPolicy getObjectAcl(String bucketName, String key) {
		GetObjectAclResponse response = runWithPermit(() -> client.getObjectAcl(GetObjectAclRequest.builder().bucket(bucketName).key(key).build()));
		return AccessControlPolicy.builder().grants(response.grants()).owner(response.owner()).build();
	}
	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#headObject
	 */
	public HeadObjectResponse getObjectMetadata(String bucketName, String key) {
		return runWithPermit(() -> client.headObject(HeadObjectRequest.builder().bucket(bucketName).key(key).build()));
	}

	public List<Tag> getObjectTags(String bucketName, String key) {
		return runWithPermit(() -> client.getObjectTagging(GetObjectTaggingRequest.builder().bucket(bucketName).key(key).build()).tagSet());
	}

	public String getObjectKmsKeyId(String bucketName, String key) {
		return getObjectMetadata(bucketName, key).ssekmsKeyId();
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#listObjectsV2Paginator
	 */
	public ListObjectsV2Iterable listObjectsV2Paginator(ListObjectsV2Request request) {
		return runWithPermit(() -> client.listObjectsV2Paginator(request));
	}

	// ===== transfer manager section =====

	synchronized S3TransferManagerSync transferManager() {
		if( transferManager==null ) {
			transferPool = ThreadPoolManager.create("S3TransferManager");
			var delegate = S3TransferManager.builder()
					.s3Client(factory.getS3AsyncClient(S3AsyncClientConfiguration.create(props), global))
					.executor(transferPool)
					.build();
			transferManager = new S3TransferManagerSync(delegate);
		}
		return transferManager;
	}

	public void downloadFile(S3Path source, File target) throws IOException {
		DownloadFileRequest downloadFileRequest = DownloadFileRequest.builder()
			.getObjectRequest(b -> b.bucket(source.getBucket()).key(source.getKey()))
			.destination(target)
			.build();

		IOException exception = runWithPermit(() -> {
			try {
				transferManager().downloadFile(downloadFileRequest);
			} catch (InterruptedException e) {
				log.debug("S3 download file: s3://{}/{} cancelled", source.getBucket(), source.getKey());
				Thread.currentThread().interrupt();
			} catch (ExecutionException e) {
				String msg = String.format("Exception thrown downloading S3 object s3://{}/{}", source.getBucket(), source.getKey());
				return new IOException(msg, e.getCause());
			}
			return null;
		});

		if( exception != null )
			throw exception;
	}

	public void downloadDirectory(S3Path source, File targetFile) throws IOException {
		DownloadDirectoryRequest downloadDirRequest = DownloadDirectoryRequest.builder()
				.bucket(source.getBucket())
				.listObjectsV2RequestTransformer(builder -> builder.prefix(source.getKey()))
				.destination(targetFile.toPath())
				.build();

		IOException exception = runWithPermit(() -> {
			try {
				CompletedDirectoryDownload completed = transferManager().downloadDirectory(downloadDirRequest);
				if( !completed.failedTransfers().isEmpty() ) {
					log.debug("S3 download directory: s3://{}/{} failed transfers", source.getBucket(), source.getKey());
					return new IOException("Some transfers in S3 download directory: s3://"+ source.getBucket() +"/"+ source.getKey() +" has failed - Transfers: " +  completed.failedTransfers() );
				}
			} catch (InterruptedException e){
				log.debug("S3 download directory: s3://{}/{} interrupted", source.getBucket(), source.getKey());
				Thread.currentThread().interrupt();
			} catch (ExecutionException e) {
				String msg = String.format("Exception thrown downloading S3 object s3://{}/{}", source.getBucket(), source.getKey());
				return new IOException(msg, e.getCause());
			}
			return null;
		});

		if( exception != null )
			throw exception;
	}

	public void uploadFile(File source, S3Path target) throws IOException{
		PutObjectRequest.Builder req = PutObjectRequest.builder().bucket(target.getBucket()).key(target.getKey());
		preparePutObjectRequest(req, target.getTagsList(), target.getContentType(), target.getStorageClass());

		UploadFileRequest uploadFileRequest = UploadFileRequest.builder()
				.putObjectRequest(req.build())
				.source(source)
				.build();

		IOException exception = runWithPermit(() -> {
			try {
				transferManager().uploadFile(uploadFileRequest);
			} catch (InterruptedException e) {
				log.debug("S3 upload file: s3://{}/{} interrupted", target.getBucket(), target.getKey());
				Thread.currentThread().interrupt();
			} catch (ExecutionException e) {
				String msg = String.format("Exception thrown uploading S3 object s3://{}/{}", target.getBucket(), target.getKey());
				return new IOException(msg, e.getCause());
			}
			return null;
		});

		if( exception != null )
			throw exception;
	}

	private Consumer<UploadFileRequest.Builder> transformUploadRequest(List<Tag> tags) {
		return builder -> builder.putObjectRequest(updateBuilder(builder.build().putObjectRequest().toBuilder(), tags).build());
	}

	private PutObjectRequest.Builder updateBuilder(PutObjectRequest.Builder porBuilder, List<Tag> tags) {

		if( cannedAcl != null )
			porBuilder.acl(cannedAcl);
		if( storageEncryption != null )
			porBuilder.serverSideEncryption(storageEncryption);
		if( kmsKeyId != null )
			porBuilder.ssekmsKeyId(kmsKeyId);
		if( tags != null && !tags.isEmpty() )
			porBuilder.tagging(Tagging.builder().tagSet(tags).build());
		return porBuilder;
	}

	public void uploadDirectory(File source, S3Path target) throws IOException {
		UploadDirectoryRequest request = UploadDirectoryRequest.builder()
				.bucket(target.getBucket())
				.s3Prefix(target.getKey())
				.source(source.toPath())
				.uploadFileRequestTransformer(transformUploadRequest(target.getTagsList()))
				.build();

		IOException exception = runWithPermit(() -> {
			try {
				CompletedDirectoryUpload completed = transferManager().uploadDirectory(request);
				if( !completed.failedTransfers().isEmpty() ) {
					log.debug("S3 upload directory: s3://{}/{} failed transfers", target.getBucket(), target.getKey());
					return new IOException("Some transfers in S3 upload directory: s3://"+ target.getBucket() +"/"+ target.getKey() +" has failed - Transfers: " +  completed.failedTransfers() );
				}
			} catch (InterruptedException e) {
				log.debug("S3 upload directory: s3://{}/{} interrupted", target.getBucket(), target.getKey());
				Thread.currentThread().interrupt();
			} catch (ExecutionException e) {
				String msg = String.format("Exception thrown uploading S3 object s3://{}/{}", target.getBucket(), target.getKey());
				return new IOException(msg, e.getCause());
			}
			return null;
		});

		if( exception != null )
			throw exception;
	}

    public void copyFile(CopyObjectRequest.Builder reqBuilder, List<Tag> tags, String contentType, String storageClass) throws IOException {
		if( tags !=null && !tags.isEmpty()) {
			log.debug("Setting tags: {}", tags);
            reqBuilder.taggingDirective(TaggingDirective.REPLACE);
            reqBuilder.tagging(Tagging.builder().tagSet(tags).build());
		}
		if( cannedAcl != null ) {
			reqBuilder.acl(cannedAcl);
		}
		if( storageEncryption != null ) {
			reqBuilder.serverSideEncryption(storageEncryption);
		}
		if( kmsKeyId !=null ) {
            reqBuilder.ssekmsKeyId(kmsKeyId);
		}
		if( contentType!=null ) {
			reqBuilder.metadataDirective(MetadataDirective.REPLACE);
            reqBuilder.contentType(contentType);
		}
		if( storageClass!=null ) {
			reqBuilder.storageClass(storageClass);
		}
        CopyObjectRequest req = reqBuilder.build();
		if( log.isTraceEnabled() ) {
			log.trace("S3 CopyObject request {}", req);
		}
		CopyRequest copyRequest = CopyRequest.builder().copyObjectRequest(req).build();

		IOException exception = runWithPermit(() -> {
			try {
				transferManager().copy(copyRequest);
			} catch (InterruptedException e) {
				log.debug("S3 copy s3://{}/{} to s3://{}/{} interrupted", req.sourceBucket(), req.sourceKey(), req.destinationBucket(), req.destinationKey());
				Thread.currentThread().interrupt();
			} catch (ExecutionException e) {
				String msg = String.format("Exception thrown copying S3 object form s3://{}/{} to s3://{}/{}", req.sourceBucket(), req.sourceKey(), req.destinationBucket(), req.destinationKey());
				return new IOException(msg, e.getCause());
			}
			return null;
		});

		if( exception != null )
			throw exception;
	}
}
