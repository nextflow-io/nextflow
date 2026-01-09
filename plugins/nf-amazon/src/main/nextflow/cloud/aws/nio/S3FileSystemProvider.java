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
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.NoSuchFileException;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributeView;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileAttribute;
import java.nio.file.attribute.FileAttributeView;
import java.nio.file.attribute.FileTime;
import java.nio.file.spi.FileSystemProvider;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import software.amazon.awssdk.services.s3.model.*;
import software.amazon.awssdk.services.s3.model.S3Object;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import nextflow.cloud.aws.AwsClientFactory;
import nextflow.cloud.aws.config.AwsConfig;
import nextflow.cloud.aws.nio.util.IOUtils;
import nextflow.cloud.aws.nio.util.S3MultipartOptions;
import nextflow.cloud.aws.nio.util.S3ObjectId;
import nextflow.cloud.aws.nio.util.S3ObjectSummaryLookup;
import nextflow.extension.FilesEx;
import nextflow.file.CopyOptions;
import nextflow.file.FileHelper;
import nextflow.file.FileSystemTransferAware;
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
public class S3FileSystemProvider extends FileSystemProvider implements FileSystemTransferAware {

	private static final Logger log = LoggerFactory.getLogger(S3FileSystemProvider.class);

	final Map<String, S3FileSystem> fileSystems = new HashMap<>();

    private final S3ObjectSummaryLookup s3ObjectSummaryLookup = new S3ObjectSummaryLookup();

	@Override
	public String getScheme() {
		return "s3";
	}

	@Override
	public FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
		Preconditions.checkNotNull(uri, "uri is null");
		Preconditions.checkArgument(uri.getScheme().equals("s3"), "uri scheme must be 's3': '%s'", uri);

		final String bucketName = S3Path.bucketName(uri);
		synchronized (fileSystems) {
			if( fileSystems.containsKey(bucketName))
				throw new FileSystemAlreadyExistsException("S3 filesystem already exists. Use getFileSystem() instead");

			final AwsConfig awsConfig = new AwsConfig(env);
			//
			final S3FileSystem result = createFileSystem(uri, awsConfig);
			fileSystems.put(bucketName, result);
			return result;
		}
	}

	@Override
	public FileSystem getFileSystem(URI uri) {
		final String bucketName = S3Path.bucketName(uri);
		final FileSystem fileSystem = this.fileSystems.get(bucketName);

		if (fileSystem == null) {
			throw new FileSystemNotFoundException("S3 filesystem not yet created. Use newFileSystem() instead");
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
		Preconditions.checkArgument(uri.getScheme().equals(getScheme()),"URI scheme must be %s", getScheme());
		return getFileSystem(uri).getPath(uri.getPath());
	}

    @Override
    public DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {

        Preconditions.checkArgument(dir instanceof S3Path,"path must be an instance of %s", S3Path.class.getName());
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
				"cannot create InputStream for root directory: %s", FilesEx.toUriString(s3Path));

		InputStream result = s3Path
					.getFileSystem()
					.getClient()
					.getObject(s3Path.getBucket(), s3Path.getKey());

        if (result == null)
			throw new IOException(String.format("The specified path is a directory: %s", FilesEx.toUriString(s3Path)));

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
						throw new FileAlreadyExistsException(FilesEx.toUriString(s3Path));
					}
				} else {
					if (!createNew && !create) {
						throw new NoSuchFileException(FilesEx.toUriString(s3Path));
					}
				}
			}
		}

		return createUploaderOutputStream(s3Path);
	}

	@Override
	public boolean canUpload(Path source, Path target) {
		return FileSystems.getDefault().equals(source.getFileSystem()) && target instanceof S3Path;
	}

	@Override
	public boolean canDownload(Path source, Path target) {
		return source instanceof S3Path && FileSystems.getDefault().equals(target.getFileSystem());
	}

	@Override
	public void download(Path remoteFile, Path localDestination, CopyOption... options) throws IOException {
		final S3Path source = (S3Path)remoteFile;

		final CopyOptions opts = CopyOptions.parse(options);
		// delete target if it exists and REPLACE_EXISTING is specified
		if (opts.replaceExisting()) {
			FileHelper.deletePath(localDestination);
		}
		else if (Files.exists(localDestination))
			throw new FileAlreadyExistsException(localDestination.toString());

		// Read S3 file attributes (metadata) for the source path, returns Optional.empty() if file doesn't exist
		final Optional<S3FileAttributes> attrs = readAttr1(source);
		// Extract directory status from attributes, defaulting to false if no attributes found
		final boolean isDir = attrs.map(S3FileAttributes::isDirectory).orElse(false);
		// Get file size only for non-directories (directories have size 0), defaulting to 0L if no attributes
		final long size = attrs.filter(a -> !a.isDirectory()).map(S3FileAttributes::size).orElse(0L);
		final String type = isDir ? "directory": "file";
		final S3Client s3Client = source.getFileSystem().getClient();
		log.debug("S3 download {} from={} to={} size={}", type, FilesEx.toUriString(source), localDestination, size);
		if( isDir ) {
			s3Client.downloadDirectory(source, localDestination.toFile());
		}
		else {
			s3Client.downloadFile(source, localDestination.toFile(), size);
		}
	}

	@Override
	public void upload(Path localFile, Path remoteDestination, CopyOption... options) throws IOException {
		final S3Path target = (S3Path) remoteDestination;

		CopyOptions opts = CopyOptions.parse(options);
		LinkOption[] linkOptions = (opts.followLinks()) ? new LinkOption[0] : new LinkOption[] { LinkOption.NOFOLLOW_LINKS };

		// attributes of source file
		if (Files.readAttributes(localFile, BasicFileAttributes.class, linkOptions).isSymbolicLink())
			throw new IOException("Uploading of symbolic links not supported - offending path: " + localFile);

		final Optional<S3FileAttributes> attrs = readAttr1(target);
		final boolean exits = attrs.isPresent();

		// delete target if it exists and REPLACE_EXISTING is specified
		if (opts.replaceExisting()) {
			FileHelper.deletePath(target);
		}
		else if ( exits )
			throw new FileAlreadyExistsException(target.toString());

		final boolean isDir = Files.isDirectory(localFile);
		final String type = isDir ? "directory": "file";
		log.debug("S3 upload {} from={} to={}", type, localFile, FilesEx.toUriString(target));
		final S3Client s3Client = target.getFileSystem().getClient();
		if( isDir ) {
			s3Client.uploadDirectory(localFile.toFile(), target);
		}
		else {
			s3Client.uploadFile(localFile.toFile(), target);
		}
	}

	private S3OutputStream createUploaderOutputStream( S3Path fileToUpload ) {
		S3Client s3 = fileToUpload.getFileSystem().getClient();
		Properties props = fileToUpload.getFileSystem().properties();

		final String storageClass = fileToUpload.getStorageClass()!=null ? fileToUpload.getStorageClass() : props.getProperty("upload_storage_class");
		final S3MultipartOptions opts = props != null ? new S3MultipartOptions(props) : new S3MultipartOptions();
		final S3ObjectId objectId = fileToUpload.toS3ObjectId();
		S3OutputStream stream = new S3OutputStream(s3.getClient(), objectId, opts)
				.setCannedAcl(s3.getCannedAcl())
				.setStorageClass(storageClass)
				.setStorageEncryption(props.getProperty("storage_encryption"))
				.setKmsKeyId(props.getProperty("storage_kms_key_id"))
				.setContentType(fileToUpload.getContentType())
				.setTags(fileToUpload.getTagsList());
		return stream;
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
					.getObject(s3Path.getBucket(), s3Path.getKey());

			if (is == null)
				throw new IOException(String.format("The specified path is a directory: %s", path));

			Files.write(tempFile, IOUtils.toByteArray(is));
		}
        catch (NoSuchFileException e){
            // Intentionally ignored. File could not exist.
		}

        // and we can use the File SeekableByteChannel implementation
		final SeekableByteChannel seekable = Files .newByteChannel(tempFile, options);
		final List<Tag> tags = ((S3Path) path).getTagsList();
		final String contentType = ((S3Path) path).getContentType();

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
                    try (InputStream stream = Files.newInputStream(tempFile)) {
                        /*
                         FIXME: if the stream is {@link InputStream#markSupported()} i can reuse the same stream
                         and evict the close and open methods of probeContentType. By this way:
                         metadata.setContentType(new Tika().detect(stream, tempFile.getFileName().toString()));
                        */
                        s3Path.getFileSystem()
                                .getClient()
                                .putObject(s3Path.getBucket(), s3Path.getKey(), stream, tags, contentType, Files.size(tempFile));
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

		List<Tag> tags = s3Path.getTagsList();

		String keyName = s3Path.getKey()
				+ (s3Path.getKey().endsWith("/") ? "" : "/");

		s3Path.getFileSystem()
				.getClient()
				.putObject(s3Path.getBucket(), keyName, new ByteArrayInputStream(new byte[0]), tags, null, 0);
	}

	@Override
	public void delete(Path path) throws IOException {
		Preconditions.checkArgument(path instanceof S3Path,
				"path must be an instance of %s", S3Path.class.getName());

		S3Path s3Path = (S3Path) path;

        if (Files.notExists(path)){
            throw new NoSuchFileException("the path: " + FilesEx.toUriString(s3Path) + " does not exist");
        }

        if (Files.isDirectory(path)){
            try (DirectoryStream<Path> stream = Files.newDirectoryStream(path)){
                if (stream.iterator().hasNext()){
                    throw new DirectoryNotEmptyException("the path: " + FilesEx.toUriString(s3Path) + " is a directory and is not empty");
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
						"target already exists: %s", FilesEx.toUriString(s3Target)));
			}
		}

		S3Client client = s3Source.getFileSystem() .getClient();
		final List<Tag> tags = ((S3Path) target).getTagsList();
		final String contentType = ((S3Path) target).getContentType();
		final String storageClass = ((S3Path) target).getStorageClass();

        //TransferManager alternative
        CopyObjectRequest.Builder reqBuilder = CopyObjectRequest.builder()
                .sourceBucket(s3Source.getBucket())
                .sourceKey(s3Source.getKey())
                .destinationBucket(s3Target.getBucket())
                .destinationKey(s3Target.getKey());
			log.trace("Copy file via copy object - source: source={}, target={}, tags={}, storageClass={}", s3Source, s3Target, tags, storageClass);
			client.copyFile(reqBuilder, tags, contentType, storageClass);
	}


	@Override
	public void move(Path source, Path target, CopyOption... options) throws IOException {
		for( CopyOption it : options ) {
			if( it==StandardCopyOption.ATOMIC_MOVE )
				throw new IllegalArgumentException("Atomic move not supported by S3 file system provider");
		}
		copy(source,target,options);
		delete(source);
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

		S3Client client = s3Path.getFileSystem().getClient();

		if( modes==null || modes.length==0 ) {
			// when no modes are given, the method is invoked
			// by `Files.exists` method, therefore just use summary lookup
			s3ObjectSummaryLookup.lookup((S3Path)path);
			return;
		}

		// get ACL and check if the file exists as a side-effect
		AccessControlPolicy acl = getAccessControl(s3Path);
        String caller = client.getCallerAccount();
		for (AccessMode accessMode : modes) {
			switch (accessMode) {
			case EXECUTE:
				throw new AccessDeniedException(s3Path.toString(), null,
						"file is not executable");
			case READ:
				if (caller == null) {
                    //if we cannot get the user's canonical ID, try read the object;
                    s3ObjectSummaryLookup.lookup((S3Path) path);
                }
                else if (!hasPermissions(acl, caller,
						EnumSet.of(Permission.FULL_CONTROL, Permission.READ))) {
					throw new AccessDeniedException(s3Path.toString(), null,
							"file is not readable");
				}
				break;
			case WRITE:
				if (caller == null) {
                    log.warn("User's Canonical Id cannot be retrieved. We can not check the access.");
                }
                else if (!hasPermissions(acl, caller,
						EnumSet.of(Permission.FULL_CONTROL, Permission.WRITE))) {
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
	private boolean hasPermissions(AccessControlPolicy acl, String owner,
			EnumSet<Permission> permissions) {
		boolean result = false;
		for (Grant grant : acl.grants()) {
			if (grant.grantee().id().equals(owner)
					&& permissions.contains(grant.permission())) {
				result = true;
				break;
			}
		}
		return result;
	}

	@Override
	public <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
		Preconditions.checkArgument(path instanceof S3Path,
				"path must be an instance of %s", S3Path.class.getName());
		S3Path s3Path = (S3Path) path;
		if (type.isAssignableFrom(BasicFileAttributeView.class)) {
			try {
				return (V) new S3FileAttributesView(readAttr0(s3Path));
			}
			catch (IOException e) {
				throw new RuntimeException("Unable read attributes for file: " + FilesEx.toUriString(s3Path), e);
			}
		}
		log.trace("Unsupported S3 file system provider file attribute view: " + type.getName());
        return null;
	}


	@Override
	public <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
		Preconditions.checkArgument(path instanceof S3Path,
				"path must be an instance of %s", S3Path.class.getName());
		S3Path s3Path = (S3Path) path;
		if (type.isAssignableFrom(BasicFileAttributes.class)) {
			return (A) ("".equals(s3Path.getKey())
					// the root bucket is implicitly a directory
					? new S3FileAttributes("/", null, 0, true, false)
					// read the target path attributes
					: readAttr0(s3Path));
		}
		// not support attribute class
		throw new UnsupportedOperationException(format("only %s supported", BasicFileAttributes.class));
	}

	private Optional<S3FileAttributes> readAttr1(S3Path s3Path) throws IOException {
		try {
			return Optional.of(readAttr0(s3Path));
		}
		catch (NoSuchFileException e) {
			return Optional.<S3FileAttributes>empty();
		}
	}

	private S3FileAttributes readAttr0(S3Path s3Path) throws IOException {
		S3Object objectSummary = s3ObjectSummaryLookup.lookup(s3Path);

		// parse the data to BasicFileAttributes.
		FileTime lastModifiedTime = null;
		if( objectSummary.lastModified() != null ) {
			lastModifiedTime = FileTime.from(objectSummary.lastModified().toEpochMilli(), TimeUnit.MILLISECONDS);
		}

		long size =  objectSummary.size();
		boolean directory = false;
		boolean regularFile = false;
		String key = objectSummary.key();
		// check if is a directory and the key of this directory exists in amazon s3
		if (objectSummary.key().equals(s3Path.getKey() + "/") && objectSummary.key().endsWith("/")) {
			directory = true;
		}
		// is a directory but does not exist in amazon s3
		else if ((!objectSummary.key().equals(s3Path.getKey()) || "".equals(s3Path.getKey())) && objectSummary.key().startsWith(s3Path.getKey())){
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

		return new S3FileAttributes(key, lastModifiedTime, size, directory, regularFile);
	}

	@Override
	public Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setAttribute(Path path, String attribute, Object value,
			LinkOption... options) throws IOException {
		throw new UnsupportedOperationException();
	}

	protected S3FileSystem createFileSystem(URI uri, AwsConfig awsConfig) {
		// try to load amazon props
		Properties props = loadAmazonProperties();
		// add properties for legacy compatibility
		props.putAll(awsConfig.getS3LegacyProperties());

		final String bucketName = S3Path.bucketName(uri);
        // do not use `global` flag for custom endpoint because
        // when enabling that flag, it overrides S3 endpoints with AWS global endpoint
        // see https://github.com/nextflow-io/nextflow/pull/5779
		final boolean global = bucketName!=null && !awsConfig.getS3Config().isCustomEndpoint();
		final AwsClientFactory factory = new AwsClientFactory(awsConfig, awsConfig.resolveS3Region());
		final S3Client client = new S3Client(factory, props, global);

		// set the client acl
		client.setCannedAcl(getProp(props, "s_3_acl", "s3_acl", "s3acl", "s3Acl"));
		client.setStorageEncryption(props.getProperty("storage_encryption"));
		client.setKmsKeyId(props.getProperty("storage_kms_key_id"));
        client.setRequesterPaysEnabled(props.getProperty("requester_pays"));

		if( props.getProperty("glacier_auto_retrieval") != null )
			log.warn("Glacier auto-retrieval is no longer supported, config option `aws.client.glacierAutoRetrieval` will be ignored");

		return new S3FileSystem(this, client, uri, props);
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
        catch(IOException e) {
			return false;
		}
	}

	/**
	 * Get the Control List, if the path does not exist
     * (because the path is a directory and this key isn't created at amazon s3)
     * then return the ACL of the first child.
     *
	 * @param path {@link S3Path}
	 * @return AccessControlList
	 * @throws IOException if error getting access control
	 */
	private AccessControlPolicy getAccessControl(S3Path path) throws IOException{
        String key = path.getKey();
        if (key == null || key.isEmpty())
            return path.getFileSystem().getClient().getBucketAcl(path.getBucket());
        return path.getFileSystem().getClient().getObjectAcl(path.getBucket(), key);
	}

    /**
     * create a temporal directory to create streams
     * @return Path temporal folder
     * @throws IOException
     */
    protected Path createTempDir() throws IOException {
        return Files.createTempDirectory("temp-s3-");
    }
}
