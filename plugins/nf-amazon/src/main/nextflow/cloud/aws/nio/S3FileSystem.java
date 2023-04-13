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

import java.io.IOException;
import java.net.URI;
import java.nio.file.FileStore;
import java.nio.file.FileSystem;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.WatchService;
import java.nio.file.attribute.UserPrincipalLookupService;
import java.nio.file.spi.FileSystemProvider;
import java.util.Properties;
import java.util.Set;

import com.amazonaws.services.s3.model.Bucket;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;

public class S3FileSystem extends FileSystem {
	
	private final S3FileSystemProvider provider;
	private final S3Client client;
	private final String endpoint;
	private final String bucketName;

	private final Properties properties;

	public S3FileSystem(S3FileSystemProvider provider, S3Client client, URI uri, Properties props) {
		this.provider = provider;
		this.client = client;
		this.endpoint = uri.getHost();
		this.bucketName = S3Path.bucketName(uri);
		this.properties = props;
	}

	@Override
	public FileSystemProvider provider() {
		return provider;
	}

	public Properties properties() {
		return properties;
	}

	@Override
	public void close() {
		this.provider.fileSystems.remove(bucketName);
	}

	@Override
	public boolean isOpen() {
		return this.provider.fileSystems.containsKey(bucketName);
	}

	@Override
	public boolean isReadOnly() {
		return false;
	}

	@Override
	public String getSeparator() {
		return S3Path.PATH_SEPARATOR;
	}

	@Override
	public Iterable<Path> getRootDirectories() {
		ImmutableList.Builder<Path> builder = ImmutableList.builder();

		for (Bucket bucket : client.listBuckets()) {
			builder.add(new S3Path(this, bucket.getName()));
		}

		return builder.build();
	}

	@Override
	public Iterable<FileStore> getFileStores() {
		return ImmutableList.of();
	}

	@Override
	public Set<String> supportedFileAttributeViews() {
		return ImmutableSet.of("basic");
	}

	@Override
	public Path getPath(String first, String... more) {
		if (more.length == 0) {
			return new S3Path(this, first);
		}

		return new S3Path(this, first, more);
	}

	@Override
	public PathMatcher getPathMatcher(String syntaxAndPattern) {
		throw new UnsupportedOperationException();
	}

	@Override
	public UserPrincipalLookupService getUserPrincipalLookupService() {
		throw new UnsupportedOperationException();
	}

	@Override
	public WatchService newWatchService() throws IOException {
		throw new UnsupportedOperationException();
	}

	public S3Client getClient() {
		return client;
	}

	/**
	 * get the endpoint associated with this fileSystem.
	 * 
	 * @see <a href="http://docs.aws.amazon.com/general/latest/gr/rande.html">http://docs.aws.amazon.com/general/latest/gr/rande.html</a>
	 * @return string
	 */
	public String getEndpoint() {
		return endpoint;
	}

	public String getBucketName() {
		return bucketName;
	}
}
