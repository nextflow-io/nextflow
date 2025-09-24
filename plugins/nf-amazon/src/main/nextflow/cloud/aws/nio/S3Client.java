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
import java.io.InterruptedIOException;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.*;
import java.util.function.Consumer;

import nextflow.cloud.aws.nio.util.ExtendedS3TransferManager;
import nextflow.cloud.aws.nio.util.S3SyncClientConfiguration;
import nextflow.extension.FilesEx;
import software.amazon.awssdk.core.ResponseInputStream;
import software.amazon.awssdk.core.sync.RequestBody;
import software.amazon.awssdk.services.s3.model.*;
import software.amazon.awssdk.services.s3.paginators.ListObjectsV2Iterable;
import software.amazon.awssdk.transfer.s3.S3TransferManager;
import software.amazon.awssdk.transfer.s3.model.*;
import nextflow.cloud.aws.AwsClientFactory;
import nextflow.cloud.aws.nio.util.S3AsyncClientConfiguration;
import nextflow.cloud.aws.util.AwsHelper;
import nextflow.util.ThreadPoolManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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

	private ExtendedS3TransferManager transferManager;

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
			log.debug("Unable to fetch caller account - {} ", e.getMessage());
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
		if( tags != null && !tags.isEmpty()) {
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
		if( tags != null && !tags.isEmpty()) {
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

	// ===== transfer manager section =====

	synchronized ExtendedS3TransferManager transferManager() {
		if( transferManager == null ) {
			log.debug("Creating S3 transfer manager pool - max-treads={};", transferManagerThreads);
			transferPool = ThreadPoolManager.create("S3TransferManager", transferManagerThreads);
			var delegate = S3TransferManager.builder()
					.s3Client(factory.getS3AsyncClient(S3AsyncClientConfiguration.create(props), global))
					.executor(transferPool)
					.build();
            transferManager = new ExtendedS3TransferManager(delegate, props);

		}
		return transferManager;
	}

    public void downloadFile(S3Path source, File target, long size) throws IOException {
		try{
            DownloadFileRequest downloadFileRequest = DownloadFileRequest.builder()
                .getObjectRequest(b -> b.bucket(source.getBucket()).key(source.getKey()))
                .destination(target)
                .build();
			transferManager().downloadFile(downloadFileRequest,size).completionFuture().get();
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new InterruptedIOException(String.format("S3 download file: s3://%s/%s cancelled", source.getBucket(), source.getKey()));
		} catch (ExecutionException e) {
			String msg = String.format("Exception thrown downloading S3 object s3://%s/%s", source.getBucket(), source.getKey());
			throw new IOException(msg, e.getCause());
		}
	}

    private static void createDirectory(Path dir) throws IOException {
        try {
            Files.createDirectory(dir);
        } catch (FileAlreadyExistsException e) {
            log.trace("File already exists: " + dir);
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
        final List<OngoingFileDownload> allDownloads = new ArrayList<>();

        FileVisitor<Path> visitor = new SimpleFileVisitor<Path>() {

            public FileVisitResult preVisitDirectory(Path current, BasicFileAttributes attr) throws IOException {
                // get the *delta* path against the source path
                final Path rel = source.relativize(current);
                final String delta = rel != null ? rel.toString() : null;
                final Path newFolder = delta != null ? target.resolve(delta) : target;
                if(log.isTraceEnabled())
                    log.trace("Download DIR: " + current + " -> " + newFolder);
                // this `copy` creates the new folder, but does not copy the contained files
                createDirectory(newFolder);
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult visitFile(Path current, BasicFileAttributes attr) throws IOException {
                // get the *delta* path against the source path
                final Path rel = source.relativize(current);
                final String delta = rel != null ? rel.toString() : null;
                final Path newFile = delta != null ? target.resolve(delta) : target;
                if( log.isTraceEnabled())
                    log.trace("Download file: " + current + " -> "+ FilesEx.toUriString(newFile));
                try {
                S3Path s3Path = (S3Path)current;
                    DownloadFileRequest downloadFileRequest = DownloadFileRequest.builder()
                        .getObjectRequest(b -> b.bucket(s3Path.getBucket()).key(s3Path.getKey()))
                        .destination(newFile)
                        .build();
                    FileDownload it = transferManager().downloadFile(downloadFileRequest, attr.size());
                    allDownloads.add(new OngoingFileDownload(s3Path.getBucket(), s3Path.getKey(), it));
                }catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new InterruptedIOException(String.format("S3 download directory: s3://%s/%s interrupted", source.getBucket(), source.getKey()));
                }
                return FileVisitResult.CONTINUE;
            }
        };

        Files.walkFileTree(source, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, visitor);

        try {
            Throwable cause = null;
            while(!allDownloads.isEmpty()) {
                OngoingFileDownload current = allDownloads.get(0);
                try{
                    current.download.completionFuture().get();
                } catch (ExecutionException e) {
                    cause = e.getCause();
                    log.debug("Exception thrown downloading S3 object s3://{}/{}", current.bucket, current.key, cause);
                }
                allDownloads.remove(0);
            }
            if (cause != null)
                throw new IOException(String.format("Some transfers from S3 download directory: s3://%s/%s has failed", source.getBucket(), source.getKey() ), cause);
        }
        catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new InterruptedIOException(String.format("Interrupted while download directory s3://%s/%s", source.getBucket(), source.getKey()));
        }
    }

	public void uploadFile(File source, S3Path target) throws IOException {
		PutObjectRequest.Builder req = PutObjectRequest.builder().bucket(target.getBucket()).key(target.getKey());
		preparePutObjectRequest(req, target.getTagsList(), target.getContentType(), target.getStorageClass());
		// initiate transfer
		FileUpload upload = transferManager().uploadFile(UploadFileRequest.builder().putObjectRequest(req.build()).source(source).build());
        try{
		    upload.completionFuture().get();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new InterruptedIOException(String.format("S3 upload file: s3://%s/%s interrupted", target.getBucket(), target.getKey()));
        } catch (ExecutionException e) {
            String msg = String.format("Exception thrown uploading S3 object s3://%s/%s", target.getBucket(), target.getKey());
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
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new InterruptedIOException(String.format("S3 upload directory: s3://%s/%s interrupted", target.getBucket(), target.getKey()));
		} catch (ExecutionException e) {
			String msg = String.format("Exception thrown uploading S3 object s3://%s/%s", target.getBucket(), target.getKey());
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
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new InterruptedIOException(String.format("S3 copy s3://%s/%s to s3://%s/%s interrupted", req.sourceBucket(), req.sourceKey(), req.destinationBucket(), req.destinationKey()));
        } catch (ExecutionException e) {
            String msg = String.format("Exception thrown copying S3 object form s3://%s/%s to s3://%s/%s", req.sourceBucket(), req.sourceKey(), req.destinationBucket(), req.destinationKey());
            throw new IOException(msg, e.getCause());
        }

	}

    static class OngoingFileDownload {
        String bucket;
        String key;
        FileDownload download;

        public OngoingFileDownload(String bucket, String key, FileDownload download) {
            this.bucket = bucket;
            this.key = key;
            this.download = download;
        }

    }
}
