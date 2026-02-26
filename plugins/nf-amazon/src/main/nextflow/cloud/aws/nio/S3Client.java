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

import nextflow.cloud.aws.nio.util.S3SyncClientConfiguration;
import software.amazon.awssdk.core.ResponseInputStream;
import software.amazon.awssdk.core.exception.SdkException;
import software.amazon.awssdk.core.sync.RequestBody;
import software.amazon.awssdk.services.s3.model.*;
import software.amazon.awssdk.services.s3.paginators.ListObjectsV2Iterable;
import software.amazon.awssdk.transfer.s3.S3TransferManager;
import software.amazon.awssdk.transfer.s3.model.*;
import nextflow.cloud.aws.AwsClientFactory;
import nextflow.cloud.aws.nio.util.S3AsyncClientConfiguration;
import nextflow.cloud.aws.nio.util.S3MultipartOptions;
import nextflow.cloud.aws.util.AwsHelper;
import nextflow.util.ThreadPoolManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static nextflow.cloud.aws.nio.util.S3UploadHelper.*;

/**
 * Client Amazon S3
 * @see software.amazon.awssdk.services.s3.S3Client
 */
public class S3Client {

	private static final Logger log = LoggerFactory.getLogger(S3Client.class);

	private software.amazon.awssdk.services.s3.S3Client client;

	private ObjectCannedACL cannedAcl;

	private String kmsKeyId;

	private ServerSideEncryption storageEncryption;

	private S3TransferManager transferManager;

	private ExecutorService transferPool;

	private Integer transferManagerThreads = 10;

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
		this.callerAccount = fetchCallerAccount();
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
			List<Bucket> buckets = client.listBuckets(ListBucketsRequest.builder().maxBuckets(1).build()).buckets();
			if (buckets == null || buckets.isEmpty())
				return null;
			return getBucketAcl(buckets.get(0).name()).owner().id();
		}catch (Throwable e){
			log.debug("Exception fetching caller account", e);
			return null;
		}
	}


	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#listBuckets()
	 */
	public List<Bucket> listBuckets() {
		return client.listBuckets().buckets();
	}
	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#listObjects(ListObjectsRequest)
	 */
	public ListObjectsResponse listObjects(ListObjectsRequest request) {
		return client.listObjects(request);
	}
	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#getObject
	 */
	public ResponseInputStream<GetObjectResponse> getObject(String bucketName, String key) {
		GetObjectRequest.Builder reqBuilder = GetObjectRequest.builder().bucket(bucketName).key(key);
		if( this.isRequesterPaysEnabled )
			reqBuilder.requestPayer(RequestPayer.REQUESTER);
		return client.getObject(reqBuilder.build());
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
		return client.putObject(builder.build(), file.toPath());
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
		return client.putObject(req, RequestBody.fromInputStream(inputStream, contentLength));
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

		client.copyObject(req);
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#getBucketAcl
	 */
	public AccessControlPolicy getBucketAcl(String bucket) {
		GetBucketAclResponse response = client.getBucketAcl(GetBucketAclRequest.builder().bucket(bucket).build());
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

	public void setTransferManagerThreads(String value) {
		if( value==null )
			return;

		try {
			this.transferManagerThreads = Integer.valueOf(value);
			log.debug("Setting S3 upload max threads={}", transferManagerThreads);
		}
		catch( NumberFormatException e ) {
			log.warn("Not a valid AWS S3 upload max threads: `{}` -- Using default", value);
		}
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
		GetObjectAclResponse response = client.getObjectAcl(GetObjectAclRequest.builder().bucket(bucketName).key(key).build());
		return AccessControlPolicy.builder().grants(response.grants()).owner(response.owner()).build();
	}
	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#headObject
	 */
	public HeadObjectResponse getObjectMetadata(String bucketName, String key) {
		return client.headObject(HeadObjectRequest.builder().bucket(bucketName).key(key).build());
	}

	public List<Tag> getObjectTags(String bucketName, String key) {
		return client.getObjectTagging(GetObjectTaggingRequest.builder().bucket(bucketName).key(key).build()).tagSet();
	}

	public String getObjectKmsKeyId(String bucketName, String key) {
		return getObjectMetadata(bucketName, key).ssekmsKeyId();
	}

	/**
	 * @see software.amazon.awssdk.services.s3.S3Client#listObjectsV2Paginator
	 */
	public ListObjectsV2Iterable listObjectsV2Paginator(ListObjectsV2Request request) {
		return client.listObjectsV2Paginator(request);
	}

	public void multipartCopyObject(S3Path s3Source, S3Path s3Target, Long objectSize, S3MultipartOptions opts, List<Tag> tags, String contentType, String storageClass ) {

		final String sourceBucketName = s3Source.getBucket();
		final String sourceObjectKey = s3Source.getKey();
		final String sourceS3Path = "s3://"+sourceBucketName+'/'+sourceObjectKey;
		final String targetBucketName = s3Target.getBucket();
		final String targetObjectKey = s3Target.getKey();

		// Step 2: Initialize
		CreateMultipartUploadRequest.Builder reqBuilder = CreateMultipartUploadRequest.builder()
				.bucket(targetBucketName)
				.key(targetObjectKey);

		if( cannedAcl!=null ) {
			reqBuilder.acl(cannedAcl);
		}
		if( storageEncryption!=null ) {
			reqBuilder.serverSideEncryption(storageEncryption);
		}
		if( kmsKeyId != null ) {
			reqBuilder.ssekmsKeyId(kmsKeyId);
		}

		if( tags != null && tags.size()>0 ) {
			reqBuilder.tagging( Tagging.builder().tagSet(tags).build() );
		}

		if( contentType!=null ) {
			reqBuilder.contentType(contentType);
		}

		if( storageClass!=null ) {
			reqBuilder.storageClass(StorageClass.fromValue(storageClass));
		}

		CreateMultipartUploadResponse initResult = client.createMultipartUpload(reqBuilder.build());


		// Step 3: Save upload Id.
		String uploadId = initResult.uploadId();

		// Multipart upload and copy allows max 10_000 parts
		// each part can be up to 5 GB
		// Max file size is 5 TB
		// See https://docs.aws.amazon.com/AmazonS3/latest/userguide/qfacts.html
		final int defChunkSize = opts.getChunkSize();
		final long partSize = computePartSize(objectSize, defChunkSize);
		ExecutorService executor = S3OutputStream.getOrCreateExecutor(opts.getMaxThreads());
		List<Callable<CompletedPart>> copyPartRequests = new ArrayList<>();
		checkPartSize(partSize);

		// Step 4. create copy part requests
		long bytePosition = 0;
		for (int i = 1; bytePosition < objectSize; i++)
		{
			checkPartIndex(i, sourceS3Path, objectSize, partSize);

			long lastPosition = bytePosition + partSize -1;
			if( lastPosition >= objectSize )
				lastPosition = objectSize - 1;

			UploadPartCopyRequest copyRequest = UploadPartCopyRequest.builder()
					.sourceBucket(sourceBucketName)
					.sourceKey(sourceObjectKey)
					.destinationBucket(targetBucketName)
					.destinationKey(targetObjectKey)
					.uploadId(uploadId)
					.partNumber(i)
					.copySourceRange("bytes=" + bytePosition + "-" + lastPosition) // e.g., "bytes=0-5242879"
					.build();

			copyPartRequests.add( copyPart(client, copyRequest, opts) );
			bytePosition += partSize;
		}

		log.trace("Starting multipart copy from: {} to {} -- uploadId={}; objectSize={}; chunkSize={}; numOfChunks={}", s3Source, s3Target, uploadId, objectSize, partSize, copyPartRequests.size() );


		List<CompletedPart> completedParts = new ArrayList<>();
		try {
			// Step 5. Start parallel parts copy
			List<Future<CompletedPart>> futures = executor.invokeAll(copyPartRequests);
			// Step 6. Fetch all results
			for (Future<CompletedPart> future : futures) {
				completedParts.add(future.get());
			}
		} catch( Exception e ) {
			throw new IllegalStateException("Multipart copy reported an unexpected error -- uploadId=" + uploadId, e);
		}

		// Step 7. Complete copy operation
		CompletedMultipartUpload completedUpload = CompletedMultipartUpload.builder()
		.parts(completedParts)
		.build();

		CompleteMultipartUploadRequest completeRequest = CompleteMultipartUploadRequest.builder()
		.bucket(targetBucketName)
		.key(targetObjectKey)
		.uploadId(uploadId)
		.multipartUpload(completedUpload)
		.build();

		log.trace("Completing multipart copy uploadId={}", uploadId);
		client.completeMultipartUpload(completeRequest);
	}

	static Callable<CompletedPart> copyPart( final software.amazon.awssdk.services.s3.S3Client client, final UploadPartCopyRequest request, final S3MultipartOptions opts ) {
		return new Callable<CompletedPart>() {
			@Override
			public CompletedPart call() throws Exception {
				return copyPart0(client,request,opts);
			}
		};
	}


	static CompletedPart copyPart0(software.amazon.awssdk.services.s3.S3Client client, UploadPartCopyRequest request, S3MultipartOptions opts) throws IOException, InterruptedException {

		final String objectId = request.uploadId();
		final int partNumber = request.partNumber();
		final String range = request.copySourceRange();

		int attempt=0;
		CompletedPart result=null;
		while( result == null ) {
			attempt++;
			try {
				log.trace("Copying multipart {} with length {} attempt {} for {} ", partNumber, range, attempt, objectId);
				UploadPartCopyResponse response = client.uploadPartCopy(request);
				result = CompletedPart.builder()
						.partNumber(partNumber)
						.eTag(response.copyPartResult().eTag())
						.build();
			}
			catch (SdkException e) {
				if( attempt >= opts.getMaxAttempts() )
					throw new IOException("Failed to upload multipart data to Amazon S3", e);

				log.debug("Failed to upload part {} attempt {} for {} -- Caused by: {}", partNumber, attempt, objectId, e.getMessage());
				Thread.sleep(opts.getRetrySleepWithAttempt(attempt));
			}
		}

		return result;
	}

	// ===== transfer manager section =====

	synchronized S3TransferManager transferManager() {
		if( transferManager==null ) {
			log.debug("Creating S3 transfer manager pool - max-treads={};", transferManagerThreads);
			transferPool = ThreadPoolManager.create("S3TransferManager", transferManagerThreads);
			transferManager = S3TransferManager.builder()
					.s3Client(factory.getS3AsyncClient(S3AsyncClientConfiguration.create(props), global))
					.executor(transferPool)
					.build();
		}
		return transferManager;
	}

	public void downloadFile(S3Path source, File target) throws IOException {
		DownloadFileRequest downloadFileRequest = DownloadFileRequest.builder()
			.getObjectRequest(b -> b.bucket(source.getBucket()).key(source.getKey()))
			.destination(target)
			.build();

		FileDownload downloadFile = transferManager().downloadFile(downloadFileRequest);
		try{
			downloadFile.completionFuture().get();
		} catch (InterruptedException e){
			log.debug("S3 download file: s3://{}/{} cancelled", source.getBucket(), source.getKey());
			Thread.currentThread().interrupt();
		} catch (ExecutionException e) {
			String msg = String.format("Exception thrown downloading S3 object s3://{}/{}", source.getBucket(), source.getKey());
			throw new IOException(msg, e.getCause());
		}

	}

	public void downloadDirectory(S3Path source, File targetFile) throws IOException {
		DownloadDirectoryRequest downloadDirRequest = DownloadDirectoryRequest.builder()
				.bucket(source.getBucket())
				.listObjectsV2RequestTransformer(builder -> builder.prefix(source.getKey()))
				.destination(targetFile.toPath())
				.build();

		DirectoryDownload downloadDirectory = transferManager().downloadDirectory(downloadDirRequest);
		try{
			downloadDirectory.completionFuture().get();
		} catch (InterruptedException e){
			log.debug("S3 download directory: s3://{}/{} interrupted", source.getBucket(), source.getKey());
			Thread.currentThread().interrupt();
		} catch (ExecutionException e) {
			String msg = String.format("Exception thrown downloading S3 object s3://{}/{}", source.getBucket(), source.getKey());
			throw new IOException(msg, e.getCause());
		}
	}

	public void uploadFile(File source, S3Path target) throws IOException{
		PutObjectRequest.Builder req = PutObjectRequest.builder().bucket(target.getBucket()).key(target.getKey());
		preparePutObjectRequest(req, target.getTagsList(), target.getContentType(), target.getStorageClass());
		// initiate transfer
		FileUpload upload = transferManager().uploadFile(UploadFileRequest.builder().putObjectRequest(req.build()).source(source).build());
        try{
		    upload.completionFuture().get();
        } catch (InterruptedException e){
            log.debug("S3 upload file: s3://{}/{} interrupted", target.getBucket(), target.getKey());
            Thread.currentThread().interrupt();
        } catch (ExecutionException e) {
            String msg = String.format("Exception thrown uploading S3 object s3://{}/{}", target.getBucket(), target.getKey());
            throw new IOException(msg, e.getCause());
        }
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

		// initiate transfer
		DirectoryUpload upload = transferManager().uploadDirectory(request);
		try {
			CompletedDirectoryUpload completed = upload.completionFuture().get();
			if (!completed.failedTransfers().isEmpty()){
				log.debug("S3 upload directory: s3://{}/{} failed transfers", target.getBucket(), target.getKey());
				throw new IOException("Some transfers in S3 upload directory: s3://"+ target.getBucket() +"/"+ target.getKey() +" has failed - Transfers: " +  completed.failedTransfers() );
			}
		} catch (InterruptedException e){
			log.debug("S3 upload directory: s3://{}/{} interrupted", target.getBucket(), target.getKey());
			Thread.currentThread().interrupt();
		} catch (ExecutionException e) {
			String msg = String.format("Exception thrown uploading S3 object s3://{}/{}", target.getBucket(), target.getKey());
			throw new IOException(msg, e.getCause());
		}
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
		Copy copy = transferManager().copy(CopyRequest.builder().copyObjectRequest(req).build());
        try {
            copy.completionFuture().get();
        } catch (InterruptedException e){
            log.debug("S3 copy s3://{}/{} to s3://{}/{} interrupted", req.sourceBucket(), req.sourceKey(), req.destinationBucket(), req.destinationKey());
			Thread.currentThread().interrupt();
        } catch (ExecutionException e) {
            String msg = String.format("Exception thrown copying S3 object form s3://{}/{} to s3://{}/{}", req.sourceBucket(), req.sourceKey(), req.destinationBucket(), req.destinationKey());
            throw new IOException(msg, e.getCause());
        }

	}
}
