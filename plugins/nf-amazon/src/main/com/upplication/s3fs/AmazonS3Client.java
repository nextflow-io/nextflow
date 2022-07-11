/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Javier Arnáiz @arnaix
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package com.upplication.s3fs;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.amazonaws.AmazonClientException;
import com.amazonaws.ClientConfiguration;
import com.amazonaws.auth.AWSCredentials;
import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.regions.Region;
import com.amazonaws.regions.RegionUtils;
import com.amazonaws.regions.Regions;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.AccessControlList;
import com.amazonaws.services.s3.model.Bucket;
import com.amazonaws.services.s3.model.CannedAccessControlList;
import com.amazonaws.services.s3.model.CompleteMultipartUploadRequest;
import com.amazonaws.services.s3.model.CopyObjectRequest;
import com.amazonaws.services.s3.model.CopyObjectResult;
import com.amazonaws.services.s3.model.CopyPartRequest;
import com.amazonaws.services.s3.model.CopyPartResult;
import com.amazonaws.services.s3.model.GetObjectTaggingRequest;
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
import com.amazonaws.services.s3.model.S3Object;
import com.amazonaws.services.s3.model.SSEAlgorithm;
import com.amazonaws.services.s3.model.SSEAwsKeyManagementParams;
import com.amazonaws.services.s3.model.Tag;
import com.amazonaws.services.s3.transfer.Download;
import com.amazonaws.services.s3.transfer.MultipleFileDownload;
import com.amazonaws.services.s3.transfer.MultipleFileUpload;
import com.amazonaws.services.s3.transfer.ObjectCannedAclProvider;
import com.amazonaws.services.s3.transfer.ObjectMetadataProvider;
import com.amazonaws.services.s3.transfer.ObjectTaggingProvider;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.TransferManagerBuilder;
import com.amazonaws.services.s3.transfer.Upload;
import com.amazonaws.services.s3.transfer.UploadContext;
import com.amazonaws.services.s3.transfer.internal.TransferManagerUtils;
import com.upplication.s3fs.util.S3MultipartOptions;
import nextflow.extension.FilesEx;
import nextflow.util.Duration;
import nextflow.util.ThreadPoolHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Client Amazon S3
 * @see com.amazonaws.services.s3.AmazonS3Client
 */
public class AmazonS3Client implements ObjectMetadataProvider, ObjectTaggingProvider, ObjectCannedAclProvider {

	private static final Logger log = LoggerFactory.getLogger(AmazonS3Client.class);
	
	private AmazonS3 client;

	private CannedAccessControlList cannedAcl;

	private String kmsKeyId;

	private SSEAlgorithm storageEncryption;

	private TransferManager transferManager;

	private ExecutorService transferPool;

	private static final List<AmazonS3Client> allSessions = new ArrayList<>();

	{
		allSessions.add(this);
	}

	public AmazonS3Client(AmazonS3 client){
		this.client = client;
	}

	public AmazonS3Client(ClientConfiguration config) {
		this.client = AmazonS3ClientBuilder.standard().withClientConfiguration(config).build();
	}

	public AmazonS3Client(ClientConfiguration config, AWSCredentials creds ) {
		this(config, creds, Regions.DEFAULT_REGION.getName());
	}

	public AmazonS3Client(ClientConfiguration config, AWSCredentials creds, String region) {
		this.client = AmazonS3ClientBuilder
				.standard()
				.withCredentials(new AWSStaticCredentialsProvider(creds))
				.withClientConfiguration(config)
				.withRegion( region )
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

	private PutObjectRequest preparePutObjectRequest(PutObjectRequest req, ObjectMetadata metadata, List<Tag> tags) {
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
		return req;
	}

	/**
	 * @see com.amazonaws.services.s3.AmazonS3Client#putObject(String, String, java.io.InputStream, ObjectMetadata)
	 */
	public PutObjectResult putObject(String bucket, String keyName, InputStream inputStream, ObjectMetadata metadata, List<Tag> tags) {
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
	public CopyObjectResult copyObject(CopyObjectRequest req, List<Tag> tags) {
		if( tags !=null && tags.size()>0 ) {
			req.setNewObjectTagging(new ObjectTagging(tags));
		}
		if( cannedAcl != null ) {
			req.withCannedAccessControlList(cannedAcl);
		}
		if( storageEncryption != null ) {
			ObjectMetadata meta = req.getNewObjectMetadata();
			meta.setSSEAlgorithm(storageEncryption.toString());
			req.setNewObjectMetadata(meta);
		}
		if( kmsKeyId !=null ) {
			req.withSSEAwsKeyManagementParams(new SSEAwsKeyManagementParams(kmsKeyId));
		}
		if( log.isTraceEnabled() ) {
			log.trace("S3 CopyObject request {}", req);
		}
		return client.copyObject(req);
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
		this.cannedAcl = CannedAccessControlList.valueOf(acl);
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

	public CannedAccessControlList getCannedAcl() {
		return cannedAcl;
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

	public void multipartCopyObject(S3Path s3Source, S3Path s3Target, Long objectSize, S3MultipartOptions opts, List<Tag> tags ) {

		final String sourceBucketName = s3Source.getBucket();
		final String sourceObjectKey = s3Source.getKey();
		final String targetBucketName = s3Target.getBucket();
		final String targetObjectKey = s3Target.getKey();

		// Step 2: Initialize
		InitiateMultipartUploadRequest initiateRequest = new InitiateMultipartUploadRequest(targetBucketName, targetObjectKey);
		if( cannedAcl!=null ) {
			initiateRequest.withCannedACL(cannedAcl);
		}
		if( storageEncryption!=null ) {
			ObjectMetadata meta = new ObjectMetadata();
			meta.setSSEAlgorithm(storageEncryption.toString());
			initiateRequest.withObjectMetadata(meta);
		}
		if( kmsKeyId != null ) {
			initiateRequest.setSSEAwsKeyManagementParams( new SSEAwsKeyManagementParams(kmsKeyId) );
		}

		if( tags != null && tags.size()>0 ) {
			initiateRequest.setTagging( new ObjectTagging(tags));
		}
		
		InitiateMultipartUploadResult initResult = client.initiateMultipartUpload(initiateRequest);


		// Step 3: Save upload Id.
		String uploadId = initResult.getUploadId();

		final int partSize = opts.getChunkSize(objectSize);
		ExecutorService executor = S3OutputStream.getOrCreateExecutor(opts.getMaxThreads());
		List<Callable<CopyPartResult>> copyPartRequests = new ArrayList<>();

		// Step 4. create copy part requests
		long bytePosition = 0;
		for (int i = 1; bytePosition < objectSize; i++)
		{
			long lastPosition = bytePosition + partSize -1 >= objectSize ? objectSize - 1 : bytePosition + partSize - 1;

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

	synchronized TransferManager transferManager() {
		if( transferManager==null ) {
			transferPool = TransferManagerUtils.createDefaultExecutorService();
			transferManager = TransferManagerBuilder.standard()
					.withS3Client(getClient())
					.withExecutorFactory(() -> transferPool)
					.build();
		}
		return transferManager;
	}

	public void downloadFile(S3Path source, File target) {
		log.debug("S3 download file from={} to={}", FilesEx.toUriString(source), target);
		Download download = transferManager()
				.download(source.getBucket(), source.getKey(), target);
		try {
			download.waitForCompletion();
		}
		catch (InterruptedException e) {
			log.debug("S3 download file: s3://{}/{} interrupted",source.getBucket(), source.getKey());
			Thread.currentThread().interrupt();
		}
	}

	public void downloadDirectory(S3Path source, File target) {
		log.debug("S3 download directory from={} to={}", FilesEx.toUriString(source), target);
		MultipleFileDownload download = transferManager()
				.downloadDirectory(source.getBucket(), source.getKey(), target);
		try {
			download.waitForCompletion();
		}
		catch (InterruptedException e) {
			log.debug("S3 download directory: s3://{}/{} interrupted",source.getBucket(), source.getKey());
			Thread.currentThread().interrupt();
		}
	}

	public void uploadFile(File source, S3Path target) {
		log.debug("S3 upload file from={} to={}", source, FilesEx.toUriString(target));
		PutObjectRequest req = new PutObjectRequest(target.getBucket(), target.getKey(), source);
		ObjectMetadata metadata = new ObjectMetadata();
		preparePutObjectRequest(req,metadata, target.getTagsList());
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

	private Map<String,List<Tag>> uploadTags = new ConcurrentHashMap<>();

	private List<Tag> uploadTags(String bucket, String key) {
		final String fullKey = bucket + "/" + key;
		return uploadTags.get(fullKey);
	}

	private void uploadTags(String bucket, String key, List<Tag> tags) {
		final String fullKey = bucket + "/" + key;
		if( tags!=null ) {
			uploadTags.put(fullKey, tags);
		}
		else
			uploadTags.remove(fullKey);
	}

	public void uploadDirectory(File source, S3Path target) {
		log.debug("S3 upload file from={} to={}", source, FilesEx.toUriString(target));
		uploadTags(target.getBucket(), target.getKey(), target.getTagsList());
		// initiate transfer
		MultipleFileUpload upload = transferManager()
				.uploadDirectory(target.getBucket(), target.getKey(), source, true, this, this, this);
		// await for completion
		try {
			upload.waitForCompletion();
		}
		catch (InterruptedException e) {
			log.debug("S3 upload file: s3://{}/{} interrupted", target.getBucket(), target.getKey());
			Thread.currentThread().interrupt();
		}
	}

	void showdown0(boolean hard) {
		if( transferManager==null )
			return;
		if( hard ) {
			transferManager.shutdownNow();
		}
		else {
			// await pool completion
			transferPool.shutdown();
			final String waitMsg = "[AWS S3] Waiting files transfer to complete (%d files)";
			final String exitMsg = "[AWS S3] Exiting before FileTransfer thread pool complete -- Some files maybe lost";
			ThreadPoolHelper.await(transferPool, Duration.of("1h"), waitMsg, exitMsg);
		}
	}

	static public void shutdown(boolean hard) {
		for( AmazonS3Client it : allSessions ) {
			try {
				it.showdown0(hard);
			}
			catch (Exception e) {
				log.debug("Unexpected error during S3 session shutdown - cause: " + e.getMessage(), e);
			}
		}
	}

	@Override
	public CannedAccessControlList provideObjectCannedAcl(File file) {
		return cannedAcl;
	}

	@Override
	public void provideObjectMetadata(File file, ObjectMetadata metadata) {
		if( storageEncryption!=null ) {
			metadata.setSSEAlgorithm(storageEncryption.toString());
		}
	}

	@Override
	public ObjectTagging provideObjectTags(UploadContext context) {
		List<Tag> tags = uploadTags(context.getBucket(), context.getKey());
		if( tags==null || tags.size()==0 )
			return null;
		return new ObjectTagging(new ArrayList<>(tags));
	}
}
