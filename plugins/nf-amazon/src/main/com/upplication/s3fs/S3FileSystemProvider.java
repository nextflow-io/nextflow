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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.ByteBuffer;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.AccessDeniedException;
import java.nio.file.AccessMode;
import java.nio.file.CopyOption;
import java.nio.file.DirectoryNotEmptyException;
import java.nio.file.DirectoryStream;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.FileStore;
import java.nio.file.FileSystem;
import java.nio.file.FileSystemAlreadyExistsException;
import java.nio.file.FileSystemNotFoundException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.NoSuchFileException;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileAttribute;
import java.nio.file.attribute.FileAttributeView;
import java.nio.file.attribute.FileTime;
import java.nio.file.spi.FileSystemProvider;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

import com.amazonaws.ClientConfiguration;
import com.amazonaws.Protocol;
import com.amazonaws.auth.AWSCredentials;
import com.amazonaws.auth.BasicAWSCredentials;
import com.amazonaws.auth.BasicSessionCredentials;
import com.amazonaws.services.s3.S3ClientOptions;
import com.amazonaws.services.s3.model.AccessControlList;
import com.amazonaws.services.s3.model.AmazonS3Exception;
import com.amazonaws.services.s3.model.CopyObjectRequest;
import com.amazonaws.services.s3.model.Grant;
import com.amazonaws.services.s3.model.ObjectMetadata;
import com.amazonaws.services.s3.model.Owner;
import com.amazonaws.services.s3.model.Permission;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import com.upplication.s3fs.ng.DownloadOpts;
import com.upplication.s3fs.ng.S3ParallelDownload;
import com.upplication.s3fs.util.IOUtils;
import com.upplication.s3fs.util.S3MultipartOptions;
import com.upplication.s3fs.util.S3ObjectSummaryLookup;
import com.upplication.s3fs.util.S3UploadRequest;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import static com.google.common.collect.Sets.difference;
import static java.lang.String.format;

/**
 * Spec:
 * 
 * URI: s3://[endpoint]/{bucket}/{key} If endpoint is missing, it's assumed to
 * be the default S3 endpoint (s3.amazonaws.com)
 * 
 * FileSystem roots: /{bucket}/
 * 
 * Treatment of S3 objects: - If a key ends in "/" it's considered a directory
 * *and* a regular file. Otherwise, it's just a regular file. - It is legal for
 * a key "xyz" and "xyz/" to exist at the same time. The latter is treated as a
 * directory. - If a file "a/b/c" exists but there's no "a" or "a/b/", these are
 * considered "implicit" directories. They can be listed, traversed and deleted.
 * 
 * Deviations from FileSystem provider API: - Deleting a file or directory
 * always succeeds, regardless of whether the file/directory existed before the
 * operation was issued i.e. Files.delete() and Files.deleteIfExists() are
 * equivalent.
 * 
 * 
 * Future versions of this provider might allow for a strict mode that mimics
 * the semantics of the FileSystem provider API on a best effort basis, at an
 * increased processing cost.
 * 
 * 
 */
public class S3FileSystemProvider extends FileSystemProvider {

	private static Logger log = LoggerFactory.getLogger(S3FileSystemProvider.class);

	public static final String ACCESS_KEY = "access_key";
	public static final String SECRET_KEY = "secret_key";

	public static final String SESSION_TOKEN = "session_token";

	final AtomicReference<S3FileSystem> fileSystem = new AtomicReference<>();

    private final S3ObjectSummaryLookup s3ObjectSummaryLookup = new S3ObjectSummaryLookup();

	private Properties props;

	private S3ParallelDownload downloader;

	@Override
	public String getScheme() {
		return "s3";
	}

	@Override
	public FileSystem newFileSystem(URI uri, Map<String, ?> env)
			throws IOException {
		Preconditions.checkNotNull(uri, "uri is null");
		Preconditions.checkArgument(uri.getScheme().equals("s3"),
				"uri scheme must be 's3': '%s'", uri);
		// first try to load amazon props
		props = loadAmazonProperties();
		Object accessKey = props.getProperty(ACCESS_KEY);
		Object secretKey = props.getProperty(SECRET_KEY);
		Object sessionToken = props.getProperty(SESSION_TOKEN);
		// but can overload by envs vars
		if (env.get(ACCESS_KEY) != null){
			accessKey = env.get(ACCESS_KEY);
		}
		if (env.get(SECRET_KEY) != null){
			secretKey = env.get(SECRET_KEY);
		}
		if (env.get(SESSION_TOKEN) != null){
			sessionToken = env.get(SESSION_TOKEN);
		}

		// allows the env variables to override the ones in the property file
		props.putAll(env);

		Preconditions.checkArgument((accessKey == null && secretKey == null)
				|| (accessKey != null && secretKey != null),
				"%s and %s (and optionally %s) should be provided or should be omitted",
				ACCESS_KEY, SECRET_KEY, SESSION_TOKEN);

		S3FileSystem result = createFileSystem(uri, accessKey, secretKey, sessionToken);

		// if this instance already has a S3FileSystem, throw exception
		// otherwise set
		if (!fileSystem.compareAndSet(null, result)) {
			throw new FileSystemAlreadyExistsException(
					"S3 filesystem already exists. Use getFileSystem() instead");
		}

		// create s3 downloader
		DownloadOpts opts = DownloadOpts.from(props, System.getenv());
		if( opts.parallelEnabled() ) {
			log.debug("Using S3 multi-part downloader");
			this.downloader = S3ParallelDownload.create(result.getClient().getClient(), opts);
		}
		else {
			log.debug("Using S3 serial downloader");
		}

		return result;
	}

	@Override
	public FileSystem getFileSystem(URI uri) {
		FileSystem fileSystem = this.fileSystem.get();

		if (fileSystem == null) {
			throw new FileSystemNotFoundException(
					String.format("S3 filesystem not yet created. Use newFileSystem() instead"));
		}

		return fileSystem;
	}

	/**
	 * Deviation from spec: throws FileSystemNotFoundException if FileSystem
	 * hasn't yet been initialized. Call newFileSystem() first.
	 * Need credentials. Maybe set credentials after? how?
	 */
	@Override
	public Path getPath(URI uri) {
		Preconditions.checkArgument(uri.getScheme().equals(getScheme()),
				"URI scheme must be %s", getScheme());

		if (uri.getHost() != null && !uri.getHost().isEmpty() &&
				!uri.getHost().equals(fileSystem.get().getEndpoint())) {
			throw new IllegalArgumentException(format(
					"only empty URI host or URI host that matching the current fileSystem: %s",
					fileSystem.get().getEndpoint())); // TODO
		}
		/**
		 * TODO: set as a list. one s3FileSystem by region
		 */
		return getFileSystem(uri).getPath(uri.getPath());
	}

    @Override
    public DirectoryStream<Path> newDirectoryStream(Path dir,
                                                    DirectoryStream.Filter<? super Path> filter) throws IOException {

        Preconditions.checkArgument(dir instanceof S3Path,
                "path must be an instance of %s", S3Path.class.getName());
        final S3Path s3Path = (S3Path) dir;

        return new DirectoryStream<Path>() {
            @Override
            public void close() throws IOException {
                // nothing to do here
            }

            @Override
            public Iterator<Path> iterator() {
                return new S3Iterator(s3Path.getFileSystem(), s3Path.getBucket(), s3Path.getKey() + "/");
            }
        };
    }

	@Override
	public InputStream newInputStream(Path path, OpenOption... options)
			throws IOException {
		Preconditions.checkArgument(options.length == 0,
				"OpenOptions not yet supported: %s",
				ImmutableList.copyOf(options)); // TODO

		Preconditions.checkArgument(path instanceof S3Path,
				"path must be an instance of %s", S3Path.class.getName());
		S3Path s3Path = (S3Path) path;

		Preconditions.checkArgument(!s3Path.getKey().equals(""),
				"cannot create InputStream for root directory: %s", s3Path);

		InputStream result;
		try {
			if( downloader==null ) {
				result = s3Path.getFileSystem().getClient()
						.getObject(s3Path.getBucket(), s3Path.getKey())
						.getObjectContent();

				if (result == null)
					throw new IOException(String.format("The specified path is a directory: %s", path));
			}
			else {
				log.debug("S3 parallel download with direct buffer: s3://{}/{}",s3Path.getBucket(), s3Path.getKey());
				result = downloader.download(s3Path.getBucket(), s3Path.getKey());
			}
		}
		catch (AmazonS3Exception e) {
			if (e.getStatusCode() == 404)
				throw new NoSuchFileException(path.toString());
			// otherwise throws a generic IO exception
			throw new IOException(String.format("Cannot access file: %s", path),e);
		}

		return result;
	}

	@Override
	public OutputStream newOutputStream(final Path path, final OpenOption... options) throws IOException {
		Preconditions.checkArgument(path instanceof S3Path, "path must be an instance of %s", S3Path.class.getName());
		S3Path s3Path = (S3Path)path;

		// validate options
		if (options.length > 0) {
			Set<OpenOption> opts = new LinkedHashSet<>(Arrays.asList(options));

			// cannot handle APPEND here -> use newByteChannel() implementation
			if (opts.contains(StandardOpenOption.APPEND)) {
				return super.newOutputStream(path, options);
			}

			if (opts.contains(StandardOpenOption.READ)) {
				throw new IllegalArgumentException("READ not allowed");
			}

			boolean create = opts.remove(StandardOpenOption.CREATE);
			boolean createNew = opts.remove(StandardOpenOption.CREATE_NEW);
			boolean truncateExisting = opts.remove(StandardOpenOption.TRUNCATE_EXISTING);

			// remove irrelevant/ignored options
			opts.remove(StandardOpenOption.WRITE);
			opts.remove(StandardOpenOption.SPARSE);

			if (!opts.isEmpty()) {
				throw new UnsupportedOperationException(opts.iterator().next() + " not supported");
			}

			if (!(create && truncateExisting)) {
				if (exists(s3Path)) {
					if (createNew || !truncateExisting) {
						throw new FileAlreadyExistsException(path.toString());
					}
				} else {
					if (!createNew && !create) {
						throw new NoSuchFileException(path.toString());
					}
				}
			}
		}

		return createUploaderOutputStream(s3Path);
	}

	private S3OutputStream createUploaderOutputStream( S3Path fileToUpload ) {
		AmazonS3Client s3 = fileToUpload.getFileSystem().getClient();

		S3UploadRequest req = props != null ? new S3UploadRequest(props) : new S3UploadRequest();
		req.setObjectId(fileToUpload.toS3ObjectId());

		S3OutputStream stream = new S3OutputStream(s3.getClient(), req);
		stream.setCannedAcl(s3.getCannedAcl());
		return stream;
	}

	protected boolean isAES256Enabled() {
		String encryption = props.getProperty("storage_encryption");
		if ( "AES256".equals(encryption) ) {
			return true;
		}
		if( encryption!=null ) {
			log.warn("Not a valid S3 server-side encryption type: `{}` -- Currently only AES256 is supported",encryption);
		}
		return false;
	}

	@Override
	public SeekableByteChannel newByteChannel(Path path,
			Set<? extends OpenOption> options, FileAttribute<?>... attrs)
			throws IOException {
		Preconditions.checkArgument(path instanceof S3Path,
				"path must be an instance of %s", S3Path.class.getName());
		final S3Path s3Path = (S3Path) path;
		// we resolve to a file inside the temp folder with the s3path name
        final Path tempFile = createTempDir().resolve(path.getFileName().toString());

		try {
			InputStream is = s3Path.getFileSystem().getClient()
					.getObject(s3Path.getBucket(), s3Path.getKey())
					.getObjectContent();

			if (is == null)
				throw new IOException(String.format("The specified path is a directory: %s", path));

			Files.write(tempFile, IOUtils.toByteArray(is));
		}
		catch (AmazonS3Exception e) {
			if (e.getStatusCode() != 404)
				throw new IOException(String.format("Cannot access file: %s", path),e);
		}

        // and we can use the File SeekableByteChannel implementation
		final SeekableByteChannel seekable = Files .newByteChannel(tempFile, options);

		return new SeekableByteChannel() {
			@Override
			public boolean isOpen() {
				return seekable.isOpen();
			}

			@Override
			public void close() throws IOException {

                if (!seekable.isOpen()) {
                    return;
                }
				seekable.close();
				// upload the content where the seekable ends (close)
                if (Files.exists(tempFile)) {
                    ObjectMetadata metadata = new ObjectMetadata();
                    metadata.setContentLength(Files.size(tempFile));
                    if( isAES256Enabled() )
                    	metadata.setSSEAlgorithm(ObjectMetadata.AES_256_SERVER_SIDE_ENCRYPTION);
                    // FIXME: #20 ServiceLoader cant load com.upplication.s3fs.util.FileTypeDetector when this library is used inside a ear :(
					metadata.setContentType(Files.probeContentType(tempFile));

                    try (InputStream stream = Files.newInputStream(tempFile)) {
                        /*
                         FIXME: if the stream is {@link InputStream#markSupported()} i can reuse the same stream
                         and evict the close and open methods of probeContentType. By this way:
                         metadata.setContentType(new Tika().detect(stream, tempFile.getFileName().toString()));
                        */
                        s3Path.getFileSystem()
                                .getClient()
                                .putObject(s3Path.getBucket(), s3Path.getKey(),
                                        stream,
                                        metadata);
                    }
                }
                else {
                    // delete: check option delete_on_close
                    s3Path.getFileSystem().
                        getClient().deleteObject(s3Path.getBucket(), s3Path.getKey());
                }
				// and delete the temp dir
                Files.deleteIfExists(tempFile);
                Files.deleteIfExists(tempFile.getParent());
			}

			@Override
			public int write(ByteBuffer src) throws IOException {
				return seekable.write(src);
			}

			@Override
			public SeekableByteChannel truncate(long size) throws IOException {
				return seekable.truncate(size);
			}

			@Override
			public long size() throws IOException {
				return seekable.size();
			}

			@Override
			public int read(ByteBuffer dst) throws IOException {
				return seekable.read(dst);
			}

			@Override
			public SeekableByteChannel position(long newPosition)
					throws IOException {
				return seekable.position(newPosition);
			}

			@Override
			public long position() throws IOException {
				return seekable.position();
			}
		};
	}

	/**
	 * Deviations from spec: Does not perform atomic check-and-create. Since a
	 * directory is just an S3 object, all directories in the hierarchy are
	 * created or it already existed.
	 */
	@Override
	public void createDirectory(Path dir, FileAttribute<?>... attrs)
			throws IOException {
		
		// FIXME: throw exception if the same key already exists at amazon s3
		
		S3Path s3Path = (S3Path) dir;

		Preconditions.checkArgument(attrs.length == 0,
				"attrs not yet supported: %s", ImmutableList.copyOf(attrs)); // TODO

		ObjectMetadata metadata = new ObjectMetadata();
		metadata.setContentLength(0);
		if( isAES256Enabled() )
			metadata.setSSEAlgorithm(ObjectMetadata.AES_256_SERVER_SIDE_ENCRYPTION);

		String keyName = s3Path.getKey()
				+ (s3Path.getKey().endsWith("/") ? "" : "/");

		s3Path.getFileSystem()
				.getClient()
				.putObject(s3Path.getBucket(), keyName,
						new ByteArrayInputStream(new byte[0]), metadata);
	}

	@Override
	public void delete(Path path) throws IOException {
		Preconditions.checkArgument(path instanceof S3Path,
				"path must be an instance of %s", S3Path.class.getName());

		S3Path s3Path = (S3Path) path;

        if (Files.notExists(path)){
            throw new NoSuchFileException("the path: " + path + " not exists");
        }

        if (Files.isDirectory(path)){
            try (DirectoryStream<Path> stream = Files.newDirectoryStream(path)){
                if (stream.iterator().hasNext()){
                    throw new DirectoryNotEmptyException("the path: " + path + " is a directory and is not empty");
                }
            }
        }

		// we delete the two objects (sometimes exists the key '/' and sometimes not)
		s3Path.getFileSystem().getClient()
			.deleteObject(s3Path.getBucket(), s3Path.getKey());
		s3Path.getFileSystem().getClient()
			.deleteObject(s3Path.getBucket(), s3Path.getKey() + "/");
	}

	@Override
	public void copy(Path source, Path target, CopyOption... options)
			throws IOException {
		Preconditions.checkArgument(source instanceof S3Path,
				"source must be an instance of %s", S3Path.class.getName());
		Preconditions.checkArgument(target instanceof S3Path,
				"target must be an instance of %s", S3Path.class.getName());

		if (isSameFile(source, target)) {
			return;
		}

		S3Path s3Source = (S3Path) source;
		S3Path s3Target = (S3Path) target;
		/*
		 * Preconditions.checkArgument(!s3Source.isDirectory(),
		 * "copying directories is not yet supported: %s", source); // TODO
		 * Preconditions.checkArgument(!s3Target.isDirectory(),
		 * "copying directories is not yet supported: %s", target); // TODO
		 */
		ImmutableSet<CopyOption> actualOptions = ImmutableSet.copyOf(options);
		verifySupportedOptions(EnumSet.of(StandardCopyOption.REPLACE_EXISTING),
				actualOptions);

		if (!actualOptions.contains(StandardCopyOption.REPLACE_EXISTING)) {
			if (exists(s3Target)) {
				throw new FileAlreadyExistsException(format(
						"target already exists: %s", target));
			}
		}

		AmazonS3Client client = s3Source.getFileSystem() .getClient();

        final ObjectMetadata sourceObjMetadata = s3Source.getFileSystem().getClient().getObjectMetadata(s3Source.getBucket(), s3Source.getKey());
		final S3MultipartOptions opts = props != null ? new S3MultipartOptions<>(props) : new S3MultipartOptions();
		final int chunkSize = opts.getChunkSize();
		final long length = sourceObjMetadata.getContentLength();

		if( length <= chunkSize ) {

			CopyObjectRequest copyObjRequest = new CopyObjectRequest(s3Source.getBucket(), s3Source.getKey(),s3Target.getBucket(), s3Target.getKey());

			ObjectMetadata targetObjectMetadata = null;
			if( isAES256Enabled() ) {
				targetObjectMetadata = new ObjectMetadata();
				targetObjectMetadata.setSSEAlgorithm(ObjectMetadata.AES_256_SERVER_SIDE_ENCRYPTION);
				copyObjRequest.setNewObjectMetadata(targetObjectMetadata);
			}
			log.trace("Copy file via copy object - source: source={}, target={}, metadata={}", s3Source, s3Target, targetObjectMetadata);
			client.copyObject(copyObjRequest);
		}
		else {
			ObjectMetadata targetObjectMetadata = null;
			if( isAES256Enabled() ) {
				targetObjectMetadata = new ObjectMetadata();
				targetObjectMetadata.setSSEAlgorithm(ObjectMetadata.AES_256_SERVER_SIDE_ENCRYPTION);
			}
			log.trace("Copy file via multipart upload - source: source={}, target={}, metadata={}", s3Source, s3Target, targetObjectMetadata);
			client.multipartCopyObject(s3Source, s3Target, length, opts, targetObjectMetadata);
		}
	}


	@Override
	public void move(Path source, Path target, CopyOption... options)
			throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isSameFile(Path path1, Path path2) throws IOException {
		return path1.isAbsolute() && path2.isAbsolute() && path1.equals(path2);
	}

	@Override
	public boolean isHidden(Path path) throws IOException {
		return false;
	}

	@Override
	public FileStore getFileStore(Path path) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void checkAccess(Path path, AccessMode... modes) throws IOException {
		S3Path s3Path = (S3Path) path;
		Preconditions.checkArgument(s3Path.isAbsolute(),
				"path must be absolute: %s", s3Path);

		AmazonS3Client client = s3Path.getFileSystem().getClient();

		if( modes==null || modes.length==0 ) {
			// when no modes are given, the method is invoked
			// by `Files.exists` method, therefore just use summary lookup
			s3ObjectSummaryLookup.lookup((S3Path)path);
			return;
		}

		// get ACL and check if the file exists as a side-effect
		AccessControlList acl = getAccessControl(s3Path);

		for (AccessMode accessMode : modes) {
			switch (accessMode) {
			case EXECUTE:
				throw new AccessDeniedException(s3Path.toString(), null,
						"file is not executable");
			case READ:
				if (!hasPermissions(acl, client.getS3AccountOwner(),
						EnumSet.of(Permission.FullControl, Permission.Read))) {
					throw new AccessDeniedException(s3Path.toString(), null,
							"file is not readable");
				}
				break;
			case WRITE:
				if (!hasPermissions(acl, client.getS3AccountOwner(),
						EnumSet.of(Permission.FullControl, Permission.Write))) {
					throw new AccessDeniedException(s3Path.toString(), null,
							format("bucket '%s' is not writable",
									s3Path.getBucket()));
				}
				break;
			}
		}
	}

    /**
     * check if the param acl has the same owner than the parameter owner and
     * have almost one of the permission set in the parameter permissions
     * @param acl
     * @param owner
     * @param permissions almost one
     * @return
     */
	private boolean hasPermissions(AccessControlList acl, Owner owner,
			EnumSet<Permission> permissions) {
		boolean result = false;
		for (Grant grant : acl.getGrants()) {
			if (grant.getGrantee().getIdentifier().equals(owner.getId())
					&& permissions.contains(grant.getPermission())) {
				result = true;
				break;
			}
		}
		return result;
	}

	@Override
	public <V extends FileAttributeView> V getFileAttributeView(Path path,
			Class<V> type, LinkOption... options) {
		throw new UnsupportedOperationException();
	}

	@Override
	public <A extends BasicFileAttributes> A readAttributes(Path path,
			Class<A> type, LinkOption... options) throws IOException {
		Preconditions.checkArgument(path instanceof S3Path,
				"path must be an instance of %s", S3Path.class.getName());
		S3Path s3Path = (S3Path) path;

		if (type == BasicFileAttributes.class) {

			S3ObjectSummary objectSummary = s3ObjectSummaryLookup.lookup(s3Path);

			// parse the data to BasicFileAttributes.
			FileTime lastModifiedTime = null;
			if( objectSummary.getLastModified() != null ) {
				lastModifiedTime = FileTime.from(objectSummary.getLastModified().getTime(), TimeUnit.MILLISECONDS);
			}

			long size =  objectSummary.getSize();
			boolean directory = false;
			boolean regularFile = false;
			String key = objectSummary.getKey();
            // check if is a directory and exists the key of this directory at amazon s3
			if (objectSummary.getKey().equals(s3Path.getKey() + "/") && objectSummary.getKey().endsWith("/")) {
				directory = true;
			}
			// is a directory but not exists at amazon s3
			else if ((!objectSummary.getKey().equals(s3Path.getKey()) || "".equals(s3Path.getKey())) && objectSummary.getKey().startsWith(s3Path.getKey())){
				directory = true;
				// no metadata, we fake one
				size = 0;
                // delete extra part
                key = s3Path.getKey() + "/";
			}
			// is a file:
			else {
                regularFile = true;
			}

			return type.cast(new S3FileAttributes(key, lastModifiedTime, size, directory, regularFile));
		}

		// not support attribute class
		throw new UnsupportedOperationException(format("only %s supported", BasicFileAttributes.class));
	}

	@Override
	public Map<String, Object> readAttributes(Path path, String attributes,
			LinkOption... options) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setAttribute(Path path, String attribute, Object value,
			LinkOption... options) throws IOException {
		throw new UnsupportedOperationException();
	}

	protected ClientConfiguration createClientConfig(Properties props) {
		ClientConfiguration config = new ClientConfiguration();

		if( props == null )
			return config;

		if( props.containsKey("connection_timeout") ) {
			log.trace("AWS client config - connection_timeout: {}", props.getProperty("connection_timeout"));
			config.setConnectionTimeout(Integer.parseInt(props.getProperty("connection_timeout")));
		}

		if( props.containsKey("max_connections")) {
			log.trace("AWS client config - max_connections: {}", props.getProperty("max_connections"));
			config.setMaxConnections(Integer.parseInt(props.getProperty("max_connections")));
		}

		if( props.containsKey("max_error_retry")) {
			log.trace("AWS client config - max_error_retry: {}", props.getProperty("max_error_retry"));
			config.setMaxErrorRetry(Integer.parseInt(props.getProperty("max_error_retry")));
		}

		if( props.containsKey("protocol")) {
			log.trace("AWS client config - protocol: {}", props.getProperty("protocol"));
			config.setProtocol(Protocol.valueOf(props.getProperty("protocol").toUpperCase()));
		}

		if( props.containsKey("proxy_domain")) {
			log.trace("AWS client config - proxy_domain: {}", props.getProperty("proxy_domain"));
			config.setProxyDomain(props.getProperty("proxy_domain"));
		}

		if( props.containsKey("proxy_host")) {
			log.trace("AWS client config - proxy_host: {}", props.getProperty("proxy_host"));
			config.setProxyHost(props.getProperty("proxy_host"));
		}

		if( props.containsKey("proxy_port")) {
			log.trace("AWS client config - proxy_port: {}", props.getProperty("proxy_port"));
			config.setProxyPort(Integer.parseInt(props.getProperty("proxy_port")));
		}

		if( props.containsKey("proxy_username")) {
			log.trace("AWS client config - proxy_username: {}", props.getProperty("proxy_username"));
			config.setProxyUsername(props.getProperty("proxy_username"));
		}

		if( props.containsKey("proxy_password")) {
			log.trace("AWS client config - proxy_password: {}", props.getProperty("proxy_password"));
			config.setProxyPassword(props.getProperty("proxy_password"));
		}

		if ( props.containsKey("proxy_workstation")) {
			log.trace("AWS client config - proxy_workstation: {}", props.getProperty("proxy_workstation"));
			config.setProxyWorkstation(props.getProperty("proxy_workstation"));
		}

		if ( props.containsKey("signer_override")) {
			log.debug("AWS client config - signerOverride: {}", props.getProperty("signer_override"));
			config.setSignerOverride(props.getProperty("signer_override"));
		}

		if( props.containsKey("socket_send_buffer_size_hints") || props.containsKey("socket_recv_buffer_size_hints") ) {
			log.trace("AWS client config - socket_send_buffer_size_hints: {}, socket_recv_buffer_size_hints: {}", props.getProperty("socket_send_buffer_size_hints","0"), props.getProperty("socket_recv_buffer_size_hints", "0"));
			int send = Integer.parseInt(props.getProperty("socket_send_buffer_size_hints","0"));
			int recv = Integer.parseInt(props.getProperty("socket_recv_buffer_size_hints", "0"));
			config.setSocketBufferSizeHints(send,recv);
		}

		if( props.containsKey("socket_timeout")) {
			log.trace("AWS client config - socket_timeout: {}", props.getProperty("socket_timeout"));
			config.setSocketTimeout(Integer.parseInt(props.getProperty("socket_timeout")));
		}

		if( props.containsKey("user_agent")) {
			log.trace("AWS client config - user_agent: {}", props.getProperty("user_agent"));
			config.setUserAgent(props.getProperty("user_agent"));
		}

		return config;
	}

	// ~~
	/**
	 * Create the fileSystem
	 * @param uri URI
	 * @param accessKey Object maybe null for anonymous authentication
	 * @param secretKey Object maybe null for anonymous authentication
	 * @return S3FileSystem never null
	 */

	protected S3FileSystem createFileSystem(URI uri, Object accessKey, Object secretKey) {
		return createFileSystem0(uri, accessKey, secretKey, null);
	}

	/**
	 * Create the fileSystem
	 * @param uri URI
	 * @param accessKey Object maybe null for anonymous authentication
	 * @param secretKey Object maybe null for anonymous authentication
	 * @param sessionToken Object maybe null for anonymous authentication
	 * @return S3FileSystem never null
	 */
	protected S3FileSystem createFileSystem(URI uri, Object accessKey, Object secretKey, Object sessionToken) {
		return createFileSystem0(uri, accessKey, secretKey, sessionToken);
	}

	protected S3FileSystem createFileSystem0(URI uri, Object accessKey, Object secretKey, Object sessionToken) {
		AmazonS3Client client;
		ClientConfiguration config = createClientConfig(props);

		if (accessKey == null && secretKey == null) {
			client = new AmazonS3Client(new com.amazonaws.services.s3.AmazonS3Client(config));
		}
		else {
			AWSCredentials credentials = (sessionToken == null
						? new BasicAWSCredentials(accessKey.toString(), secretKey.toString())
						: new BasicSessionCredentials(accessKey.toString(), secretKey.toString(), sessionToken.toString()) );

			if( System.getenv("NXF_AWS_REGION")!=null ) {
				log.debug("Creating AWS S3 client with region: " + System.getenv("NXF_AWS_REGION") );
				client = new AmazonS3Client( config, credentials, System.getenv("NXF_AWS_REGION") );
			}
			else {
				client = new AmazonS3Client(new com.amazonaws.services.s3.AmazonS3Client(credentials,config));
			}

		}

		// note: path style access is going to be deprecated
		// https://aws.amazon.com/blogs/aws/amazon-s3-path-deprecation-plan-the-rest-of-the-story/
		boolean usePathStyle = "true".equals(props.getProperty("s_3_path_style_access")) || "true".equals(props.getProperty("s3_path_style_access"));
		if (usePathStyle) {
			S3ClientOptions options = S3ClientOptions.builder()
					.setPathStyleAccess(usePathStyle)
					.build();
			client.getClient().setS3ClientOptions(options);
		}

		// set the client acl
		client.setCannedAcl(getProp(props, "s_3_acl", "s3_acl", "s3Acl"));

		if (uri.getHost() != null) {
			client.setEndpoint(uri.getHost());
		}
		else if( props.getProperty("endpoint") != null ){
			client.setEndpoint(props.getProperty("endpoint"));
		}

		return new S3FileSystem(this, client, uri.getHost());
	}

	protected String getProp(Properties props, String... keys) {
		for( String k : keys ) {
			if( props.containsKey(k) ) {
				return props.getProperty(k);
			}
		}
		return null;
	}
	
	/**
	 * find /amazon.properties in the classpath
	 * @return Properties amazon.properties
	 */
	protected Properties loadAmazonProperties() {
		Properties props = new Properties();
		// http://www.javaworld.com/javaworld/javaqa/2003-06/01-qa-0606-load.html
		// http://www.javaworld.com/javaqa/2003-08/01-qa-0808-property.html
		try(InputStream in = Thread.currentThread().getContextClassLoader().getResourceAsStream("amazon.properties")){
			if (in != null){
				props.load(in);
			}
			
		} catch (IOException e) {}
		
		return props;
	}
	
	// ~~~

	private <T> void verifySupportedOptions(Set<? extends T> allowedOptions,
			Set<? extends T> actualOptions) {
		Sets.SetView<? extends T> unsupported = difference(actualOptions,
				allowedOptions);
		Preconditions.checkArgument(unsupported.isEmpty(),
				"the following options are not supported: %s", unsupported);
	}
	/**
	 * check that the paths exists or not
	 * @param path S3Path
	 * @return true if exists
	 */
	private boolean exists(S3Path path) {
		try {
            s3ObjectSummaryLookup.lookup(path);
			return true;
		}
        catch(NoSuchFileException e) {
			return false;
		}
	}

	/**
	 * Get the Control List, if the path not exists
     * (because the path is a directory and this key isnt created at amazon s3)
     * then return the ACL of the first child.
     *
	 * @param path {@link S3Path}
	 * @return AccessControlList
	 * @throws NoSuchFileException if not found the path and any child
	 */
	private AccessControlList getAccessControl(S3Path path) throws NoSuchFileException{
		S3ObjectSummary obj = s3ObjectSummaryLookup.lookup(path);
		// check first for file:
        return path.getFileSystem().getClient().getObjectAcl(obj.getBucketName(), obj.getKey());
	}

    /**
     * create a temporal directory to create streams
     * @return Path temporal folder
     * @throws IOException
     */
    protected Path createTempDir() throws IOException {
        return Files.createTempDirectory("temp-s3-");
    }

    public static void shutdown(boolean hard) {
		S3OutputStream.shutdownExecutor(hard);
		S3ParallelDownload.shutdown(hard);
	}
}
