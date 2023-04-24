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
import java.nio.file.FileVisitOption;
import java.nio.file.FileVisitResult;
import java.nio.file.FileVisitor;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.amazonaws.AmazonClientException;
import com.amazonaws.ClientConfiguration;
import com.amazonaws.auth.AWSCredentials;
import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.regions.Region;
import com.amazonaws.regions.RegionUtils;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.Headers;
import com.amazonaws.services.s3.model.AccessControlList;
import com.amazonaws.services.s3.model.AmazonS3Exception;
import com.amazonaws.services.s3.model.Bucket;
import com.amazonaws.services.s3.model.CannedAccessControlList;
import com.amazonaws.services.s3.model.CompleteMultipartUploadRequest;
import com.amazonaws.services.s3.model.CopyObjectRequest;
import com.amazonaws.services.s3.model.CopyPartRequest;
import com.amazonaws.services.s3.model.CopyPartResult;
import com.amazonaws.services.s3.model.GetObjectTaggingRequest;
import com.amazonaws.services.s3.model.GlacierJobParameters;
import com.amazonaws.services.s3.model.InitiateMultipartUploadRequest;
import com.amazonaws.services.s3.model.InitiateMultipartUploadResult;
import com.amazonaws.services.s3.model.ListObjectsRequest;
import com.amazonaws.services.s3.model.ObjectListing;
import com.amazonaws.services.s3.model.ObjectMetadata;
import com.amazonaws.services.s3.model.ObjectTagging;
import com.amazonaws.services.s3.model.Owner;
import com.amazonaws.services.s3.model.PartETag;
import com.amazonaws.services.s3.model.PutObjectRequest;
import com.amazonaws.services.s3.model.PutObjectResult;
import com.amazonaws.services.s3.model.RestoreObjectRequest;
import com.amazonaws.services.s3.model.S3Object;
import com.amazonaws.services.s3.model.SSEAlgorithm;
import com.amazonaws.services.s3.model.SSEAwsKeyManagementParams;
import com.amazonaws.services.s3.model.StorageClass;
import com.amazonaws.services.s3.model.Tag;
import com.amazonaws.services.s3.transfer.Download;
import com.amazonaws.services.s3.transfer.MultipleFileUpload;
import com.amazonaws.services.s3.transfer.ObjectCannedAclProvider;
import com.amazonaws.services.s3.transfer.ObjectMetadataProvider;
import com.amazonaws.services.s3.transfer.ObjectTaggingProvider;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.TransferManagerBuilder;
import com.amazonaws.services.s3.transfer.Upload;
import com.amazonaws.services.s3.transfer.UploadContext;
import nextflow.cloud.aws.nio.util.S3MultipartOptions;
import nextflow.cloud.aws.util.AwsHelper;
import nextflow.util.Duration;
import nextflow.util.ThreadPoolHelper;
import nextflow.util.ThreadPoolManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static nextflow.cloud.aws.nio.util.S3UploadHelper.*;

/**
 * Client Amazon S3
 * @see com.amazonaws.services.s3.AmazonS3Client
 */
public class S3Client {

	private static final Logger log = LoggerFactory.getLogger(S3Client.class);
	
	private AmazonS3 client;

	private CannedAccessControlList cannedAcl;

	private String kmsKeyId;

	private SSEAlgorithm storageEncryption;

	private TransferManager transferManager;

	private ExecutorService transferPool;

	private Long uploadChunkSize = Long.valueOf(S3MultipartOptions.DEFAULT_CHUNK_SIZE);

	private Integer uploadMaxThreads = 10;

	private boolean glacierAutoRetrieval;

	private int glacierExpirationDays = 7;

	private String glacierRetrievalTier;

	public S3Client(AmazonS3 client) {
		this.client = client;
	}

	public S3Client(ClientConfiguration config, AWSCredentials creds, String region) {
		this.client = AmazonS3ClientBuilder
				.standard()
				.withCredentials(new AWSStaticCredentialsProvider(creds))
				.withClientConfiguration(config)
				.withRegion(region)
				.build();
	}

	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#listBuckets()
	 */
	public List<Bucket> listBuckets() {
		return client.listBuckets();
	}
	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#listObjects(ListObjectsRequest)
	 */
	public ObjectListing listObjects(ListObjectsRequest request) {
		return client.listObjects(request);
	}
	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#getObject(String, String)
	 */
	public S3Object getObject(String bucketName, String key) {
		return client.getObject(bucketName, key);
	}
	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#putObject(String, String, File)
	 */
	public PutObjectResult putObject(String bucket, String key, File file) {
		PutObjectRequest req = new PutObjectRequest(bucket, key, file);
		if( cannedAcl != null ) {
			log.trace("Setting canned ACL={}; bucket={}; key={}", cannedAcl, bucket, key);
			req.withCannedAcl(cannedAcl);
		}
		return client.putObject(req);
	}

	private PutObjectRequest preparePutObjectRequest(PutObjectRequest req, ObjectMetadata metadata, List<Tag> tags, String contentType, String storageClass) {
		req.withMetadata(metadata);
		if( cannedAcl != null ) {
			req.withCannedAcl(cannedAcl);
		}
		if( tags != null && tags.size()>0 ) {
			req.setTagging(new ObjectTagging(tags));
		}
		if( kmsKeyId != null ) {
			req.withSSEAwsKeyManagementParams( new SSEAwsKeyManagementParams(kmsKeyId) );
		}
		if( storageEncryption!=null ) {
			metadata.setSSEAlgorithm(storageEncryption.toString());
		}
		if( contentType!=null ) {
			metadata.setContentType(contentType);
		}
		if( storageClass!=null ) {
			req.setStorageClass(storageClass);
		}
		return req;
	}

	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#putObject(String, String, java.io.InputStream, ObjectMetadata)
	 */
	public PutObjectResult putObject(String bucket, String keyName, InputStream inputStream, ObjectMetadata metadata, List<Tag> tags, String contentType) {
		PutObjectRequest req = new PutObjectRequest(bucket, keyName, inputStream, metadata);
		if( cannedAcl != null ) {
			req.withCannedAcl(cannedAcl);
		}
		if( tags != null && tags.size()>0 ) {
			req.setTagging(new ObjectTagging(tags));
		}
		if( kmsKeyId != null ) {
			req.withSSEAwsKeyManagementParams( new SSEAwsKeyManagementParams(kmsKeyId) );
		}
		if( storageEncryption!=null ) {
			metadata.setSSEAlgorithm(storageEncryption.toString());
		}
		if( contentType!=null ) {
			metadata.setContentType(contentType);
		}
		if( log.isTraceEnabled() ) {
			log.trace("S3 PutObject request {}", req);
		}
		return client.putObject(req);
	}
	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#deleteObject(String, String)
	 */
	public void deleteObject(String bucket, String key) {
		client.deleteObject(bucket, key);
	}

	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#copyObject(CopyObjectRequest)
	 */
	public void copyObject(CopyObjectRequest req, List<Tag> tags, String contentType, String storageClass) {
		if( tags !=null && tags.size()>0 ) {
			req.setNewObjectTagging(new ObjectTagging(tags));
		}
		if( cannedAcl != null ) {
			req.withCannedAccessControlList(cannedAcl);
		}
		// getNewObjectMetada returns null if no object metadata has been specified.
		ObjectMetadata meta = req.getNewObjectMetadata() != null ? req.getNewObjectMetadata() : new ObjectMetadata();
		if( storageEncryption != null ) {
			meta.setSSEAlgorithm(storageEncryption.toString());
			req.setNewObjectMetadata(meta);
		}
		if( kmsKeyId !=null ) {
			req.withSSEAwsKeyManagementParams(new SSEAwsKeyManagementParams(kmsKeyId));
		}
		if( contentType!=null ) {
			meta.setContentType(contentType);
			req.setNewObjectMetadata(meta);
		}
		if( storageClass!=null ) {
			req.setStorageClass(storageClass);
		}
		if( log.isTraceEnabled() ) {
			log.trace("S3 CopyObject request {}", req);
		}

		client.copyObject(req);
	}

	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#getBucketAcl(String)
	 */
	public AccessControlList getBucketAcl(String bucket) {
		return client.getBucketAcl(bucket);
	}
	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#getS3AccountOwner()
	 */
	public Owner getS3AccountOwner() {
		return client.getS3AccountOwner();
	}
	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#setEndpoint(String)
	 */
	public void setEndpoint(String endpoint) {
		client.setEndpoint(endpoint);
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
		this.storageEncryption = SSEAlgorithm.fromString(alg);
		log.debug("Setting S3 SSE storage encryption algorithm={}", alg);
	}

	public void setUploadChunkSize(String value) {
		if( value==null )
			return;

		try {
			this.uploadChunkSize = Long.valueOf(value);
			log.debug("Setting S3 upload chunk size={}", uploadChunkSize);
		}
		catch( NumberFormatException e ) {
			log.warn("Not a valid AWS S3 upload chunk size: `{}` -- Using default", value);
		}
	}

	public void setUploadMaxThreads(String value) {
		if( value==null )
			return;

		try {
			this.uploadMaxThreads = Integer.valueOf(value);
			log.debug("Setting S3 upload max threads={}", uploadMaxThreads);
		}
		catch( NumberFormatException e ) {
			log.warn("Not a valid AWS S3 upload max threads: `{}` -- Using default", value);
		}
	}

	public CannedAccessControlList getCannedAcl() {
		return cannedAcl;
	}

	public void setGlacierAutoRetrieval(boolean value) {
		this.glacierAutoRetrieval = value;
		log.debug("Setting S3 glacierAutoRetrieval={}", glacierAutoRetrieval);
	}

	public void setGlacierAutoRetrieval(String value) {
		if( value==null )
			return;
	  	setGlacierAutoRetrieval(Boolean.parseBoolean(value));
	}

	public void setGlacierExpirationDays(int days) {
		this.glacierExpirationDays = days;
		log.debug("Setting S3 glacierExpirationDays={}", glacierExpirationDays);
	}

	public void setGlacierExpirationDays(String days) {
		if( days==null )
			return;
		try {
			setGlacierExpirationDays(Integer.parseInt(days));
		}
		catch( NumberFormatException e ) {
			log.warn("Not a valid AWS S3 glacierExpirationDays: `{}` -- Using default", days);
		}
	}

	public void setGlacierRetrievalTier(String tier) {
		this.glacierRetrievalTier = tier;
		log.debug("Setting S3 glacierRetrievalTier={}", glacierRetrievalTier);
	}

	public AmazonS3 getClient() {
		return client;
	}

	public void setRegion(String regionName) {
		Region region = RegionUtils.getRegion(regionName);
		if( region == null )
			throw new IllegalArgumentException("Not a valid S3 region name: " + regionName);
		client.setRegion(region);
	}


	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#getObjectAcl(String, String)
	 */
	public AccessControlList getObjectAcl(String bucketName, String key) {
		return client.getObjectAcl(bucketName, key);
	}
	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#getObjectMetadata(String, String)
	 */
	public ObjectMetadata getObjectMetadata(String bucketName, String key) {
		return client.getObjectMetadata(bucketName, key);
	}

	public List<Tag> getObjectTags(String bucketName, String key) {
		return client.getObjectTagging(new GetObjectTaggingRequest(bucketName,key)).getTagSet();
	}

	/**
     * @see com.amazonaws.services.s3.AmazonS3Client#listNextBatchOfObjects(com.amazonaws.services.s3.model.ObjectListing)
     */
    public ObjectListing listNextBatchOfObjects(ObjectListing objectListing) {
        return client.listNextBatchOfObjects(objectListing);
    }


	public void multipartCopyObject(S3Path s3Source, S3Path s3Target, Long objectSize, S3MultipartOptions opts, List<Tag> tags, String contentType, String storageClass ) {

		final String sourceBucketName = s3Source.getBucket();
		final String sourceObjectKey = s3Source.getKey();
		final String sourceS3Path = "s3://"+sourceBucketName+'/'+sourceObjectKey;
		final String targetBucketName = s3Target.getBucket();
		final String targetObjectKey = s3Target.getKey();
	  	final ObjectMetadata meta = new ObjectMetadata();

		// Step 2: Initialize
		InitiateMultipartUploadRequest initiateRequest = new InitiateMultipartUploadRequest(targetBucketName, targetObjectKey);
		if( cannedAcl!=null ) {
			initiateRequest.withCannedACL(cannedAcl);
		}
		if( storageEncryption!=null ) {
			meta.setSSEAlgorithm(storageEncryption.toString());
			initiateRequest.withObjectMetadata(meta);
		}
		if( kmsKeyId != null ) {
			initiateRequest.setSSEAwsKeyManagementParams( new SSEAwsKeyManagementParams(kmsKeyId) );
		}

		if( tags != null && tags.size()>0 ) {
			initiateRequest.setTagging( new ObjectTagging(tags));
		}

		if( contentType!=null ) {
			meta.setContentType(contentType);
			initiateRequest.withObjectMetadata(meta);
		}

		if( storageClass!=null ) {
			initiateRequest.setStorageClass(StorageClass.fromValue(storageClass));
		}

		InitiateMultipartUploadResult initResult = client.initiateMultipartUpload(initiateRequest);


		// Step 3: Save upload Id.
		String uploadId = initResult.getUploadId();

		// Multipart upload and copy allows max 10_000 parts
		// each part can be up to 5 GB
		// Max file size is 5 TB
		// See https://docs.aws.amazon.com/AmazonS3/latest/userguide/qfacts.html
		final int defChunkSize = opts.getChunkSize();
		final long partSize = computePartSize(objectSize, defChunkSize);
		ExecutorService executor = S3OutputStream.getOrCreateExecutor(opts.getMaxThreads());
		List<Callable<CopyPartResult>> copyPartRequests = new ArrayList<>();
		checkPartSize(partSize);

		// Step 4. create copy part requests
		long bytePosition = 0;
		for (int i = 1; bytePosition < objectSize; i++)
		{
			checkPartIndex(i, sourceS3Path, objectSize, partSize);

			long lastPosition = bytePosition + partSize -1;
			if( lastPosition >= objectSize )
				lastPosition = objectSize - 1;

			CopyPartRequest copyRequest = new CopyPartRequest()
					.withDestinationBucketName(targetBucketName)
					.withDestinationKey(targetObjectKey)
					.withSourceBucketName(sourceBucketName)
					.withSourceKey(sourceObjectKey)
					.withUploadId(uploadId)
					.withFirstByte(bytePosition)
					.withLastByte(lastPosition)
					.withPartNumber(i);

			copyPartRequests.add( copyPart(client, copyRequest, opts) );
			bytePosition += partSize;
		}

		log.trace("Starting multipart copy from: {} to {} -- uploadId={}; objectSize={}; chunkSize={}; numOfChunks={}", s3Source, s3Target, uploadId, objectSize, partSize, copyPartRequests.size() );

		List<PartETag> etags = new ArrayList<>();
		List<Future<CopyPartResult>> responses;
		try {
			// Step 5. Start parallel parts copy
			responses = executor.invokeAll(copyPartRequests);

			// Step 6. Fetch all results
			for (Future<CopyPartResult> response : responses) {
				CopyPartResult result = response.get();
				etags.add(new PartETag(result.getPartNumber(), result.getETag()));
			}
		}
		catch( Exception e ) {
			throw new IllegalStateException("Multipart copy reported an unexpected error -- uploadId=" + uploadId, e);
		}

		// Step 7. Complete copy operation
		CompleteMultipartUploadRequest completeRequest = new
				CompleteMultipartUploadRequest(
				targetBucketName,
				targetObjectKey,
				initResult.getUploadId(),
				etags);

		log.trace("Completing multipart copy uploadId={}", uploadId);
		client.completeMultipartUpload(completeRequest);
	}

	static Callable<CopyPartResult> copyPart( final AmazonS3 client, final CopyPartRequest request, final S3MultipartOptions opts ) {
		return new Callable<CopyPartResult>() {
			@Override
			public CopyPartResult call() throws Exception {
				return copyPart0(client,request,opts);
			}
		};
	}


	static CopyPartResult copyPart0(AmazonS3 client, CopyPartRequest request, S3MultipartOptions opts) throws IOException, InterruptedException {

		final String objectId = request.getUploadId();
		final int partNumber = request.getPartNumber();
		final long len = request.getLastByte() - request.getFirstByte();

		int attempt=0;
		CopyPartResult result=null;
		while( result == null ) {
			attempt++;
			try {
				log.trace("Copying multipart {} with length {} attempt {} for {} ", partNumber, len, attempt, objectId);
				result = client.copyPart(request);
			}
			catch (AmazonClientException e) {
				if( attempt >= opts.getMaxAttempts() )
					throw new IOException("Failed to upload multipart data to Amazon S3", e);

				log.debug("Failed to upload part {} attempt {} for {} -- Caused by: {}", partNumber, attempt, objectId, e.getMessage());
				Thread.sleep(opts.getRetrySleepWithAttempt(attempt));
			}
		}

		return result;
	}

	// ===== transfer manager section =====

	synchronized TransferManager transferManager() {
		if( transferManager==null ) {
			log.debug("Creating S3 transfer manager pool - chunk-size={}; max-treads={};", uploadChunkSize, uploadMaxThreads);
			transferPool = ThreadPoolManager.create("S3TransferManager", uploadMaxThreads);
			transferManager = TransferManagerBuilder.standard()
					.withS3Client(getClient())
					.withMinimumUploadPartSize(uploadChunkSize)
					.withExecutorFactory(() -> transferPool)
					.build();
		}
		return transferManager;
	}

	public void downloadFile(S3Path source, File target) {
		Download download = transferManager()
				.download(source.getBucket(), source.getKey(), target);
		try {
			download.waitForCompletion();
		}
		catch (InterruptedException e) {
			log.debug("S3 download file: s3://{}/{} interrupted",source.getBucket(), source.getKey());
			Thread.currentThread().interrupt();
		}
		catch (AmazonS3Exception e) {
			handleAmazonException(source, target, e);
		}
	}

	private void handleAmazonException(S3Path source, File target, AmazonS3Exception e) {
		// the following message is returned when accessing a Glacier stored file
		// "The operation is not valid for the object's storage class"
		final boolean isGlacierError = e.getMessage().contains("storage class")
				&& e.getErrorCode().equals("InvalidObjectState");

		if( isGlacierError && glacierAutoRetrieval ) {
			log.info("S3 download s3://{}/{} failed due to invalid storage class -- Retrieving from Glacier", source.getBucket(), source.getKey());
			restoreFromGlacier(source.getBucket(), source.getKey());
			downloadFile(source, target);
		}
		else {
			throw e;
		}
	}

	protected void restoreFromGlacier(String bucketName, String key) {
		final int sleepMillis = 30_000;
		final long _5_mins = 5 * 60 * 1_000;

		try {
			RestoreObjectRequest request = new RestoreObjectRequest(bucketName, key);

			String storageClass = client.getObjectMetadata(bucketName, key).getStorageClass();
			if( storageClass!=null && !storageClass.equals("INTELLIGENT_TIERING") )
				request.setExpirationInDays(glacierExpirationDays);

			if( glacierRetrievalTier != null )
				request.setGlacierJobParameters(
					new GlacierJobParameters()
						.withTier(glacierRetrievalTier)
				);

			client.restoreObjectV2(request);
		}
		catch (AmazonS3Exception e) {
			if( e.getMessage().contains("RestoreAlreadyInProgress") ) {
				log.debug("S3 Glacier restore already initiated for object s3://{}/{}", bucketName, key);
			}
			else {
				throw e;
			}
		}

		try {
			boolean ongoingRestore = true;
			long begin = System.currentTimeMillis();
			while( ongoingRestore ) {
				final long now = System.currentTimeMillis();
				if( now-begin>_5_mins ) {
					log.info("S3 Glacier restore ongoing for object s3://{}/{}", bucketName, key);
					begin = now;
				}
				Thread.sleep(sleepMillis);
				ongoingRestore = client.getObjectMetadata(bucketName, key).getOngoingRestore();
			}
		}
		catch (InterruptedException e) {
			log.debug("S3 Glacier restore s3://{}/{} interrupted", bucketName, key);
			Thread.currentThread().interrupt();
		}
	}

	public void downloadDirectory(S3Path source, File targetFile) throws IOException {
		//
		// the download directory method provided by the TransferManager replicates
		// the source files directory structure in the target path
		// see https://github.com/aws/aws-sdk-java/issues/1321
		//
		// just traverse to source path a copy all files
		// 
		final Path target = targetFile.toPath();
		final List<Download> allDownloads = new ArrayList<>();

		FileVisitor<Path> visitor = new SimpleFileVisitor<Path>() {

			public FileVisitResult preVisitDirectory(Path current, BasicFileAttributes attr) throws IOException {
				// get the *delta* path against the source path
				Path rel = source.relativize(current);
				String delta = rel != null ? rel.toString() : null;
				Path newFolder = delta != null ? target.resolve(delta) : target;
				if(log.isTraceEnabled())
					log.trace("Copy DIR: " + current + " -> " + newFolder);
				// this `copy` creates the new folder, but does not copy the contained files
				Files.createDirectory(newFolder);
				return FileVisitResult.CONTINUE;
			}

			@Override
			public FileVisitResult visitFile(Path current, BasicFileAttributes attr) {
				// get the *delta* path against the source path
				Path rel = source.relativize(current);
				String delta = rel != null ? rel.toString() : null;
				Path newFile = delta != null ? target.resolve(delta) : target;
				if( log.isTraceEnabled())
					log.trace("Copy file: " + current + " -> "+newFile.toUri());

				String sourceKey = ((S3Path) current).getKey();
				Download it = transferManager() .download(source.getBucket(), sourceKey, newFile.toFile());
				allDownloads.add(it);

				return FileVisitResult.CONTINUE;
			}

		};

		Files.walkFileTree(source, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, visitor);

		try {
			while(allDownloads.size()>0) {
				allDownloads.get(0).waitForCompletion();
				allDownloads.remove(0);
			}
		}
		catch (InterruptedException e) {
			log.debug("S3 download directory: s3://{}/{} interrupted", source.getBucket(), source.getKey());
			Thread.currentThread().interrupt();
		}
	}

	public void uploadFile(File source, S3Path target) {
		PutObjectRequest req = new PutObjectRequest(target.getBucket(), target.getKey(), source);
		ObjectMetadata metadata = new ObjectMetadata();
		preparePutObjectRequest(req, metadata, target.getTagsList(), target.getContentType(), target.getStorageClass());
		// initiate transfer
		Upload upload = transferManager() .upload(req);
		// await for completion
		try {
			upload.waitForCompletion();
		}
		catch (InterruptedException e) {
			log.debug("S3 upload file: s3://{}/{} interrupted", target.getBucket(), target.getKey());
			Thread.currentThread().interrupt();
		}
	}

	/**
	 * This class is used by the upload directory operation to acquire the mecessary meta info
	 */
	private class MetadataProvider implements ObjectMetadataProvider, ObjectTaggingProvider, ObjectCannedAclProvider {

		@Override
		public CannedAccessControlList provideObjectCannedAcl(File file) {
			return cannedAcl;
		}

		@Override
		public void provideObjectMetadata(File file, ObjectMetadata metadata) {
			if( storageEncryption!=null ) {
				metadata.setSSEAlgorithm(storageEncryption.toString());
			}
			if( kmsKeyId!=null ) {
				// metadata.setHeader(Headers.SERVER_SIDE_ENCRYPTION, SSEAlgorithm.KMS.getAlgorithm());
				metadata.setHeader(Headers.SERVER_SIDE_ENCRYPTION_AWS_KMS_KEYID, kmsKeyId);
			}
		}

		@Override
		public ObjectTagging provideObjectTags(UploadContext context) {
			List<Tag> tags = uploadTags.get();
			if( tags==null || tags.size()==0 )
				return null;
			return new ObjectTagging(new ArrayList<>(tags));
		}
	}

	final private MetadataProvider metaProvider = new MetadataProvider();

	final private ThreadLocal<List<Tag>> uploadTags = new ThreadLocal<>();

	public void uploadDirectory(File source, S3Path target) {
		// set the tags to be used in a thread local
		uploadTags.set( target.getTagsList() );
		// initiate transfer
		MultipleFileUpload upload = transferManager()
				.uploadDirectory(target.getBucket(), target.getKey(), source, true, metaProvider, metaProvider, metaProvider);
		// the tags are fetched by the previous operation
		// the thread local can be cleared
		uploadTags.remove();
		// await for completion
		try {
			upload.waitForCompletion();
		}
		catch (InterruptedException e) {
			log.debug("S3 upload file: s3://{}/{} interrupted", target.getBucket(), target.getKey());
			Thread.currentThread().interrupt();
		}
	}

	String getObjectKmsKeyId(String bucketName, String key) {
		return getObjectMetadata(bucketName,key).getSSEAwsKmsKeyId();
	}

	protected void showdownTransferPool(boolean hard) {
		log.debug("Initiating transfer manager shutdown (hard={})", hard);
		if( hard ) {
			transferPool.shutdownNow();
		}
		else {
			// await pool completion
			transferPool.shutdown();
			final String waitMsg = "[AWS S3] Waiting files transfer to complete (%d files)";
			final String exitMsg = "[AWS S3] Exiting before FileTransfer thread pool complete -- Some files maybe lost";
			ThreadPoolHelper.await(transferPool, Duration.of("1h"), waitMsg, exitMsg);
		}
	}
}
