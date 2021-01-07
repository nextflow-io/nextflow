// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.CoreUtils;
import com.azure.core.util.logging.ClientLogger;
import com.azure.core.util.polling.SyncPoller;
import com.azure.storage.blob.models.BlobRequestConditions;
import com.azure.storage.blob.models.BlobStorageException;
import com.azure.storage.blob.models.BlobCopyInfo;
import com.azure.storage.blob.models.BlobErrorCode;
import com.azure.storage.blob.models.ParallelTransferOptions;
import reactor.core.publisher.Flux;
import reactor.core.publisher.Mono;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.UncheckedIOException;
import java.net.HttpURLConnection;
import java.net.URI;
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
import java.nio.file.LinkOption;
import java.nio.file.NoSuchFileException;
import java.nio.file.NotDirectoryException;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.BasicFileAttributeView;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileAttribute;
import java.nio.file.attribute.FileAttributeView;
import java.nio.file.spi.FileSystemProvider;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.function.Consumer;
import java.util.function.Supplier;

/**
 * The {@code AzureFileSystemProvider} is Azure Storage's implementation of the nio interface on top of Azure Blob
 * Storage.
 * <p>
 * Particular care should be taken when working with a remote storage service. This implementation makes no guarantees
 * on behavior or state should other processes operate on the same data concurrently; file systems from this provider
 * will assume they have exclusive access to their data and will behave without regard for potential of interfering
 * applications. Moreover, remote file stores introduce higher latencies. Therefore, additional consideration should be
 * given to managing concurrency: race conditions are more likely to manifest and network failures occur more frequently
 * than disk failures. These and other such distributed application scenarios must be considered when working with this
 * file system. While the {@code AzureFileSystem} will ensure it takes appropriate steps towards robustness and
 * reliability, the application developer must design around these failure scenarios and have fallback and retry options
 * available.
 * <p>
 * The Azure Blob Storage service backing these APIs is not a true FileSystem, nor is it the goal of this implementation
 * to force Azure Blob Storage to act like a full-fledged file system. Some APIs and scenarios will remain unsupported
 * indefinitely until they may be sensibly implemented. Other APIs may experience lower performance than is expected
 * because of the number of network requests needed to ensure correctness. The javadocs for each type and method should
 * also be read carefully to understand what guarantees are made and how they may differ from the contract defined by
 * {@link FileSystemProvider}.
 * <p>
 * The scheme for this provider is {@code "azb"}, and the format of the URI to identify an {@code AzureFileSystem} is
 * {@code "azb://?account=<accountName>"}. The name of the Storage account is used to uniquely identify the file
 * system.
 * <p>
 * An {@link AzureFileSystem} is backed by an account. An {@link AzureFileStore} is backed by a container. Any number of
 * containers may be specified as file stores upon creation of the file system. When a file system is created,
 * it will try to retrieve the properties of each container to ensure connection to the account. If any of the
 * containers does not exist, it will be created. Failure to access or create containers as necessary will result in
 * an exception and failure to create the file system. Any data existing in the containers will be preserved and
 * accessible via the file system, though customers should be aware that it must be in a format understandable by
 * the types in this package or behavior will be undefined.
 * <p>
 * {@link #newFileSystem(URI, Map)} will check for the following keys in the configuration map and expect the named
 * types. Any entries not listed here will be ignored. Note that {@link AzureFileSystem} has public constants defined
 * for each of the keys for convenience.
 * <ul>
 *     <li>{@code AzureStorageAccountKey:}{@link String}</li>
 *     <li>{@code AzureStorageSasToken:}{@link String}</li>
 *     <li>{@code AzureStorageHttpLogDetailLevel:}{@link com.azure.core.http.policy.HttpLogDetailLevel}</li>
 *     <li>{@code AzureStorageMaxTries:}{@link Integer}</li>
 *     <li>{@code AzureStorageTryTimeout:}{@link Integer}</li>
 *     <li>{@code AzureStorageRetryDelayInMs:}{@link Long}</li>
 *     <li>{@code AzureStorageMaxRetryDelayInMs:}{@link Long}</li>
 *     <li>{@code AzureStorageRetryPolicyType:}{@link com.azure.storage.common.policy.RetryPolicyType}</li>
 *     <li>{@code AzureStorageSecondaryHost:}{@link String}</li>
 *     <li>{@code AzureStorageSecondaryHost:}{@link Integer}</li>
 *     <li>{@code AzureStorageBlockSize:}{@link Long}</li>
 *     <li>{@code AzureStoragePutBlobThreshold:}{@link Long}</li>
 *     <li>{@code AzureStorageMaxConcurrencyPerRequest:}{@link Integer}</li>
 *     <li>{@code AzureStorageDownloadResumeRetries:}{@link Integer}</li>
 *     <li>{@code AzureStorageUseHttps:}{@link Boolean}</li>
 *     <li>{@code AzureStorageFileStores:}{@link String}</li>
 * </ul>
 * <p>
 * Either an account key or a sas token must be specified. If both are provided, the account key will be preferred. If
 * a sas token is specified, the customer must take care that it has appropriate permissions to perform the actions
 * demanded of the file system in a given workflow, including the initial connection check specified above. Furthermore,
 * it must have an expiry time that lasts at least until the file system is closed as there is no token refresh offered
 * at this time. The same token will be applied to all containers.
 * <p>
 * An iterable of file stores must also be provided; each entry should simply be the name of a container. The first
 * container listed will be considered the default file store and the root directory of which will be the file system's
 * default directory. All other values listed are used to configure the underlying
 * {@link com.azure.storage.blob.BlobServiceClient}. Please refer to that type for more information on these values.
 *
 * @see FileSystemProvider
 */
public final class AzureFileSystemProvider extends FileSystemProvider {
    private final ClientLogger logger = new ClientLogger(AzureFileSystemProvider.class);

    /**
     * A helper for setting the HTTP properties when creating a directory.
     */
    public static final String CONTENT_TYPE = "Content-Type";

    /**
     * A helper for setting the HTTP properties when creating a directory.
     */
    public static final String CONTENT_DISPOSITION = "Content-Disposition";

    /**
     * A helper for setting the HTTP properties when creating a directory.
     */
    public static final String CONTENT_LANGUAGE = "Content-Language";

    /**
     * A helper for setting the HTTP properties when creating a directory.
     */
    public static final String CONTENT_ENCODING = "Content-Encoding";

    /**
     * A helper for setting the HTTP properties when creating a directory.
     */
    public static final String CONTENT_MD5 = "Content-MD5";

    /**
     * A helper for setting the HTTP properties when creating a directory.
     */
    public static final String CACHE_CONTROL = "Cache-Control";

    private static final String ACCOUNT_QUERY_KEY = "account";
    private static final int COPY_TIMEOUT_SECONDS = 30;

    private final ConcurrentMap<String, FileSystem> openFileSystems;


    // Specs require a public zero argument constructor.
    /**
     * Creates an AzureFileSystemProvider.
     */
    public AzureFileSystemProvider() {
        this.openFileSystems = new ConcurrentHashMap<>();
    }

    /**
     * Returns the URI scheme that identifies this provider: {@code "azb".}
     *
     * @return {@code "azb"}
     */
    @Override
    public String getScheme() {
        return "azb";
    }

    /**
     * Constructs a new FileSystem object identified by a URI.
     * <p>
     * The format of a {@code URI} identifying a file system is {@code "azb://?account=<accountName>"}.
     * <p>
     * Once closed, a file system with the same identifier may be reopened.
     *
     * @param uri URI reference
     * @param config A map of provider specific properties to configure the file system
     * @return a new file system.
     * @throws IllegalArgumentException If the pre-conditions for the uri parameter aren't met, or the env parameter
     * does not contain properties required by the provider, or a property value is invalid.
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     * @throws FileSystemAlreadyExistsException If the file system has already been created.
     */
    @Override
    public FileSystem newFileSystem(URI uri, Map<String, ?> config) throws IOException {
        String accountName = extractAccountName(uri);

        if (this.openFileSystems.containsKey(accountName)) {
            throw LoggingUtility.logError(this.logger, new FileSystemAlreadyExistsException("Name: " + accountName));
        }

        AzureFileSystem afs = new AzureFileSystem(this, accountName, config);
        this.openFileSystems.put(accountName, afs);

        return afs;
    }

    /**
     * Returns an existing FileSystem created by this provider.
     * <p>
     * The format of a {@code URI} identifying an file system is {@code "azb://?account=&lt;accountName&gt;"}.
     * <p>
     * Trying to retrieve a closed file system will throw a {@link FileSystemNotFoundException}. Once closed, a
     * file system with the same identifier may be reopened.
     *
     * @param uri URI reference
     * @return the file system
     * @throws IllegalArgumentException If the pre-conditions for the uri parameter aren't met
     * @throws FileSystemNotFoundException If the file system already exists
     * @throws SecurityException never
     */
    @Override
    public FileSystem getFileSystem(URI uri) {
        String accountName = extractAccountName(uri);
        if (!this.openFileSystems.containsKey(accountName)) {
            throw LoggingUtility.logError(this.logger, new FileSystemNotFoundException("Name: " + accountName));
        }
        return this.openFileSystems.get(accountName);
    }

    /**
     * Return a Path object by converting the given URI. The resulting Path is associated with a FileSystem that already
     * exists.
     *
     * @param uri The URI to convert
     * @return The path identified by the URI.
     * @throws IllegalArgumentException If the URI scheme does not identify this provider or other preconditions on the
     * uri parameter do not hold
     * @throws FileSystemNotFoundException if the file system identified by the query does not exist
     * @throws SecurityException never
     *
     * @see #getFileSystem(URI) for information on the URI format
     */
    @Override
    public Path getPath(URI uri) {
        return getFileSystem(uri).getPath(uri.getPath());
    }

    /**
     * Unsupported. Use {@link #newInputStream(Path, OpenOption...)} or {@link #newOutputStream(Path, OpenOption...)}
     * instead.
     *
     * @param path the Path
     * @param set open options
     * @param fileAttributes attributes
     * @return a new seekable byte channel
     * @throws UnsupportedOperationException Operation is not supported.
     * @throws IllegalArgumentException if the set contains an invalid combination of options
     * @throws FileAlreadyExistsException if a file of that name already exists and the CREATE_NEW option is specified
     * (optional specific exception)
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> set,
            FileAttribute<?>... fileAttributes) throws IOException {
        throw LoggingUtility.logError(logger, new UnsupportedOperationException());
    }

    /**
     * Opens an {@link InputStream} to the given path.
     * <p>
     * The stream will not attempt to read or buffer the entire file. However, when fetching data, it will always
     * request the same size chunk of several MB to prevent network thrashing on small reads. Mark and reset are
     * supported.
     * <p>
     * Only {@link StandardOpenOption#READ} is supported. Any other option will throw.
     *
     * @param path the path to the file to open
     * @param options options specifying how the file is opened
     * @return a new input stream
     * @throws IllegalArgumentException if an invalid combination of options is specified
     * @throws UnsupportedOperationException if an unsupported option is specified
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        // Validate options. Only read is supported.
        if (options.length > 1 || (options.length > 0 && !options[0].equals(StandardOpenOption.READ))) {
            throw LoggingUtility.logError(logger,
                new UnsupportedOperationException("Only the read option is supported."));
        }

        AzureResource resource = new AzureResource(path);
        AzurePath.ensureFileSystemOpen(resource.getPath());

        // Ensure the path points to a file.
        if (!resource.checkDirStatus().equals(DirectoryStatus.NOT_A_DIRECTORY)) {
            throw LoggingUtility.logError(logger, new IOException("Path either does not exist or points to a directory."
                + "Path must point to a file. Path: " + path.toString()));
        }

        // Note that methods on BlobInputStream are already synchronized.
        return new NioBlobInputStream(resource.getBlobClient().openInputStream(), resource.getPath());
    }

    /**
     * Opens an {@link OutputStream} to the given path. The resulting file will be stored as a block blob.
     * <p>
     * The only supported options are {@link StandardOpenOption#CREATE}, {@link StandardOpenOption#CREATE_NEW},
     * {@link StandardOpenOption#WRITE}, {@link StandardOpenOption#TRUNCATE_EXISTING}. Any other options will throw an
     * {@link UnsupportedOperationException}. {@code WRITE} and {@code TRUNCATE_EXISTING} must be specified or an
     * {@link IllegalArgumentException} will be thrown. Hence, files cannot be updated, only overwritten completely.
     * <p>
     * This stream will not attempt to buffer the entire file, however some buffering will be done for potential
     * optimizations and to avoid network thrashing. Specifically, up to
     * {@link AzureFileSystem#AZURE_STORAGE_PUT_BLOB_THRESHOLD} bytes will be buffered initially. If that threshold is
     * exceeded, the data will be broken into chunks and sent in blocks, and writes will be buffered into sizes of
     * {@link AzureFileSystem#AZURE_STORAGE_UPLOAD_BLOCK_SIZE}. The maximum number of buffers of this size to be
     * allocated is defined by {@link AzureFileSystem#AZURE_STORAGE_MAX_CONCURRENCY_PER_REQUEST}, which also configures
     * the level of parallelism with which we may write and thus may affect write speeds as well.
     * <p>
     * The data is only committed when the steam is closed. Hence data cannot be read from the destination until the
     * stream is closed. When the close method returns, it is guaranteed that, barring any errors, the data is finalized
     * and available for reading.
     * <p>
     * Writing happens asynchronously. Bytes passed for writing are stored until either the threshold or block size are
     * met at which time they are sent to the service. When the write method returns, there is no guarantee about which
     * phase of this process the data is in other than it has been accepted and will be written. Again, closing will
     * guarantee that the data is written and available.
     * <p>
     * Flush is a no-op as regards data transfers, but it can be used to check the state of the stream for errors.
     * This can be a useful tool because writing happens asynchronously, and therefore an error from a previous write
     * may not otherwise be thrown unless the stream is flushed, closed, or written to again.
     *
     * @param path the path to the file to open or create
     * @param options options specifying how the file is opened
     * @return a new output stream
     * @throws IllegalArgumentException if an invalid combination of options is specified
     * @throws UnsupportedOperationException if an unsupported option is specified
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public OutputStream newOutputStream(Path path, OpenOption... options) throws IOException {
        // If options are empty, add Create, Write, TruncateExisting as defaults per nio docs.
        if (options == null || options.length == 0) {
            options = new OpenOption[] {
                StandardOpenOption.CREATE,
                StandardOpenOption.WRITE,
                StandardOpenOption.TRUNCATE_EXISTING };
        }
        List<OpenOption> optionsList = Arrays.asList(options);

        // Check for unsupported options.
        List<OpenOption> supportedOptions = Arrays.asList(
            StandardOpenOption.CREATE_NEW,
            StandardOpenOption.CREATE,
            StandardOpenOption.WRITE,
            // Though we don't actually truncate, the same result is achieved by overwriting the destination.
            StandardOpenOption.TRUNCATE_EXISTING);
        for (OpenOption option : optionsList) {
            if (!supportedOptions.contains(option)) {
                throw LoggingUtility.logError(logger, new UnsupportedOperationException("Unsupported option: "
                    + option.toString()));
            }
        }

        // Write and truncate must be specified
        if (!optionsList.contains(StandardOpenOption.WRITE)
            || !optionsList.contains(StandardOpenOption.TRUNCATE_EXISTING)) {
            throw LoggingUtility.logError(logger,
                new IllegalArgumentException("Write and TruncateExisting must be specified to open an OutputStream"));
        }

        AzureResource resource = new AzureResource(path);
        AzurePath.ensureFileSystemOpen(resource.getPath());
        DirectoryStatus status = resource.checkDirStatus();

        // Cannot write to a directory.
        if (DirectoryStatus.isDirectory(status)) {
            throw LoggingUtility.logError(logger, new IOException("Cannot open an OutputStream to a directory. Path: "
                + path.toString()));
        }

        // Writing to an empty location requires a create option.
        if (status.equals(DirectoryStatus.DOES_NOT_EXIST)
            && !(optionsList.contains(StandardOpenOption.CREATE)
            || optionsList.contains(StandardOpenOption.CREATE_NEW))) {
            throw LoggingUtility.logError(logger, new IOException("Writing to an empty location requires a create "
                + "option. Path: " + path.toString()));
        }

        // Cannot write to an existing file if create new was specified.
        if (status.equals(DirectoryStatus.NOT_A_DIRECTORY) && optionsList.contains(StandardOpenOption.CREATE_NEW)) {
            throw LoggingUtility.logError(logger, new IOException("A file already exists at this location and "
                + "CREATE_NEW was specified. Path: " + path.toString()));
        }

        // Create options based on file system config
        AzureFileSystem fs = (AzureFileSystem) (path.getFileSystem());
        Integer blockSize = fs.getBlockSize() == null ? null : fs.getBlockSize().intValue();
        Integer putBlobThreshold = fs.getPutBlobThreshold() == null ? null : fs.getPutBlobThreshold().intValue();
        ParallelTransferOptions pto = new ParallelTransferOptions(blockSize, fs.getMaxConcurrencyPerRequest(), null,
            putBlobThreshold);

        // Add an extra etag check for create new
        BlobRequestConditions rq = null;
        if (optionsList.contains(StandardOpenOption.CREATE_NEW)) {
            rq = new BlobRequestConditions().setIfNoneMatch("*");
        }

        return new NioBlobOutputStream(resource.getBlobClient().getBlockBlobClient().getBlobOutputStream(pto, null,
            null, null, rq), resource.getPath());
    }

    /**
     * Returns an {@link AzureDirectoryStream} for iterating over the contents of a directory. The elements returned by
     * the directory stream's iterator are of type Path, each one representing an entry in the directory. The Path
     * objects are obtained as if by resolving the name of the directory entry against dir. The entries returned by the
     * iterator are filtered by the given filter.
     * <p>
     * When not using the try-with-resources construct, then directory stream's close method should be invoked after
     * iteration is completed so as to free any resources held for the open directory.
     * <p>
     * Where the filter terminates due to an uncaught error or runtime exception then it is propagated to the hasNext or
     * next method. Where an IOException is thrown, it results in the hasNext or next method throwing a
     * DirectoryIteratorException with the IOException as the cause.
     *
     * @param path the path to the directory
     * @param filter the directory stream filter
     * @return a new and open {@code DirectoryStream} object
     * @throws IllegalArgumentException If the path type is not an instance of {@link AzurePath}.
     * @throws NotDirectoryException if the file could not otherwise be opened because it is not a directory
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public DirectoryStream<Path> newDirectoryStream(Path path, DirectoryStream.Filter<? super Path> filter)
        throws IOException {
        if (!(path instanceof AzurePath)) {
            throw LoggingUtility.logError(logger, new IllegalArgumentException("This provider cannot operate on "
                + "subtypes of Path other than AzurePath"));
        }
        AzurePath.ensureFileSystemOpen(path);

        /*
        Ensure the path is a directory. Note that roots are always directories. The case of an invalid root will be
        caught in instatiating the stream below.

        Possible optimization later is to save the result of the list call to use as the first list call inside the
        stream rather than a list call for checking the status and a list call for listing.
         */
        if (!((AzurePath) path).isRoot() && !(new AzureResource(path).checkDirectoryExists())) {
            throw LoggingUtility.logError(logger, new NotDirectoryException(path.toString()));
        }

        return new AzureDirectoryStream((AzurePath) path, filter);
    }

    /**
     * Creates a new directory at the specified path.
     * <p>
     * The existence of a directory in the {@code AzureFileSystem} is defined on two levels. <i>Weak existence</i> is
     * defined by the presence of a non-zero number of blobs prefixed with the directory's path. This concept is also
     * known as a  <i>virtual directory</i> and enables the file system to work with containers that were pre-loaded
     * with data by another source but need to be accessed by this file system. <i>Strong existence</i> is defined as
     * the presence of an actual storage resource at the given path, which in the case of directories, is a zero-length
     * blob whose name is the directory path with a particular metadata field indicating the blob's status as a
     * directory. This is also known as a <i>concrete directory</i>. Directories created by this file system will
     * strongly exist. Operations targeting directories themselves as the object (e.g. setting properties) will target
     * marker blobs underlying concrete directories. Other operations (e.g. listing) will operate on the blob-name
     * prefix.
     * <p>
     * This method fulfills the nio contract of: "The check for the existence of the file and the creation of the
     * directory if it does not exist are a single operation that is atomic with respect to all other filesystem
     * activities that might affect the directory." More specifically, this method will atomically check for <i>strong
     * existence</i> of another file or directory at the given path and fail if one is present. On the other hand, we
     * only check for <i>weak existence</i> of the parent to determine if the given path is valid. Additionally, the
     * action of checking whether the parent exists, is <i>not</i> atomic with the creation of the directory. Note that
     * while it is possible that the parent may be deleted between when the parent is determined to exist and the
     * creation of the child, the creation of the child will always ensure the existence of a virtual parent, so the
     * child will never be left floating and unreachable. The different checks on parent and child is due to limitations
     * in the Storage service API.
     * <p>
     * There may be some unintuitive behavior when working with directories in this file system, particularly virtual
     * directories (usually those not created by this file system). A virtual directory will disappear as soon as all
     * its children have been deleted. Furthermore, if a directory with the given path weakly exists at the time of
     * calling this method, this method will still return success and create a concrete directory at the target
     * location. In other words, it is possible to "double create" a directory if it first weakly exists and then is
     * strongly created. This is both because it is impossible to atomically check if a virtual directory exists while
     * creating a concrete directory and because such behavior will have minimal side effects--no files will be
     * overwritten and the directory will still be available for writing as intended, though it may not be empty. This
     * is not a complete list of such unintuitive behavior.
     * <p>
     * This method will attempt to extract standard HTTP content headers from the list of file attributes to set them
     * as blob headers. All other attributes will be set as blob metadata. The value of every attribute will be
     * converted to a {@code String} with the exception of the Content-MD5 attribute which expects a {@code byte[]}.
     * When extracting the content headers, the following strings will be used for comparison (constants for these
     * values can be found on this type):
     * <ul>
     *     <li>{@code Content-Type}</li>
     *     <li>{@code Content-Disposition}</li>
     *     <li>{@code Content-Language}</li>
     *     <li>{@code Content-Encoding}</li>
     *     <li>{@code Content-MD5}</li>
     *     <li>{@code Cache-Control}</li>
     * </ul>
     * Note that these properties also have a particular semantic in that if one is specified, all are updated. In other
     * words, if any of the above is set, all those that are not set will be cleared. See the
     * <a href="https://docs.microsoft.com/en-us/rest/api/storageservices/set-blob-properties">Azure Docs</a> for more
     * information.
     *
     * @param path the directory to create
     * @param fileAttributes an optional list of file attributes to set atomically when creating the directory
     * @throws IllegalArgumentException If the path type is not an instance of {@link AzurePath}.
     * @throws UnsupportedOperationException if the array contains an attribute that cannot be set atomically when
     * creating the directory
     * @throws FileAlreadyExistsException if a directory could not otherwise be created because a file of that name
     * already exists
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public void createDirectory(Path path, FileAttribute<?>... fileAttributes) throws IOException {
        fileAttributes = fileAttributes == null ? new FileAttribute<?>[0] : fileAttributes;

        // Get the destination for the directory. Will throw if path is a root.
        AzureResource azureResource = new AzureResource(path);
        AzurePath.ensureFileSystemOpen(azureResource.getPath());

        // Check if parent exists. If it does, atomically check if a file already exists and create a new dir if not.
        if (azureResource.checkParentDirectoryExists()) {
            try {
                azureResource.setFileAttributes(Arrays.asList(fileAttributes))
                    .putDirectoryBlob(new BlobRequestConditions().setIfNoneMatch("*"));
            } catch (BlobStorageException e) {
                if (e.getStatusCode() == HttpURLConnection.HTTP_CONFLICT
                    && e.getErrorCode().equals(BlobErrorCode.BLOB_ALREADY_EXISTS)) {
                    throw LoggingUtility.logError(logger,
                        new FileAlreadyExistsException(azureResource.getPath().toString()));
                } else {
                    throw LoggingUtility.logError(logger,
                        new IOException("An error occurred when creating the directory", e));
                }
            }
        } else {
            throw LoggingUtility.logError(logger, new IOException("Parent directory does not exist for path: "
                + azureResource.getPath().toString()));
        }
    }

    /**
     * Deletes the specified resource.
     * <p>
     * This method is not atomic with respect to other file system operations. It is possible to delete a file in use by
     * another process, and doing so will not immediately invalidate any channels open to that file--they will simply
     * start to fail. Root directories cannot be deleted even when empty.
     *
     * @param path the path to the file to delete
     * @throws IllegalArgumentException If the path type is not an instance of {@link AzurePath}.
     * @throws NoSuchFileException if the file does not exist
     * @throws DirectoryNotEmptyException if the file is a directory and could not otherwise be deleted because the
     * directory is not empty
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public void delete(Path path) throws IOException {
        // Basic validation. Must be an AzurePath. Cannot be a root.
        AzureResource azureResource = new AzureResource(path);
        AzurePath.ensureFileSystemOpen(azureResource.getPath());

        // Check directory status--possibly throw DirectoryNotEmpty or NoSuchFile.
        DirectoryStatus dirStatus = azureResource.checkDirStatus();
        if (dirStatus.equals(DirectoryStatus.DOES_NOT_EXIST)) {
            throw LoggingUtility.logError(logger, new NoSuchFileException(path.toString()));
        }
        if (dirStatus.equals(DirectoryStatus.NOT_EMPTY)) {
            throw LoggingUtility.logError(logger, new DirectoryNotEmptyException(path.toString()));
        }

        // After all validation has completed, delete the resource.
        try {
            azureResource.getBlobClient().delete();
        } catch (BlobStorageException e) {
            if (e.getErrorCode().equals(BlobErrorCode.BLOB_NOT_FOUND)) {
                throw LoggingUtility.logError(logger, new NoSuchFileException(path.toString()));
            }
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    /**
     * Copies the resource at the source location to the destination.
     * <p>
     * This method is not atomic with respect to other file system operations. More specifically, the checks necessary
     * to validate the inputs and state of the file system are not atomic with the actual copying of data. If the copy
     * is triggered, the copy itself is atomic and only a complete copy will ever be left at the destination.
     * <p>
     * In addition to those in the docs for {@link FileSystemProvider#copy(Path, Path, CopyOption...)}, this method has
     * the following requirements for successful completion. {@link StandardCopyOption#COPY_ATTRIBUTES} must be passed
     * as it is impossible not to copy blob properties; if this option is not passed, an
     * {@link UnsupportedOperationException} will be thrown. Neither the source nor the destination can be a root
     * directory; if either is a root directory, an {@link IllegalArgumentException} will be thrown. The parent
     * directory of the destination must at least weakly exist; if it does not, an {@link IOException} will be thrown.
     * The only supported option other than {@link StandardCopyOption#COPY_ATTRIBUTES} is
     * {@link StandardCopyOption#REPLACE_EXISTING}; the presence of any other option will result in an
     * {@link UnsupportedOperationException}.
     * <p>
     * This method supports both virtual and concrete directories as both the source and destination. Unlike when
     * creating a directory, the existence of a virtual directory at the destination will cause this operation to fail.
     * This is in order to prevent the possibility of overwriting a non-empty virtual directory with a file. Still, as
     * mentioned above, this check is not atomic with the creation of the resultant directory.
     *
     * @param source the path to the file to copy
     * @param destination the path to the target file
     * @param copyOptions specifying how the copy should be done
     * @throws UnsupportedOperationException if the array contains a copy option that is not supported
     * @throws FileAlreadyExistsException if the target file exists but cannot be replaced because the REPLACE_EXISTING
     * option is not specified
     * @throws DirectoryNotEmptyException the REPLACE_EXISTING option is specified but the file cannot be replaced
     * because it is a non-empty directory
     * @throws IOException If an I/O error occurs.
     * @throws IllegalArgumentException If the path type is not an instance of {@link AzurePath}.
     * @throws SecurityException never
     * @see #createDirectory(Path, FileAttribute[]) for more information about directory existence.
     */
    @Override
    public void copy(Path source, Path destination, CopyOption... copyOptions) throws IOException {
        // If paths point to the same file, operation is a no-op.
        if (source.equals(destination)) {
            return;
        }

        // Read and validate options.
        // Remove accepted options as we find them. Anything left we don't support.
        boolean replaceExisting = false;
        List<CopyOption> optionsList = new ArrayList<>(Arrays.asList(copyOptions));
        if (!optionsList.contains(StandardCopyOption.COPY_ATTRIBUTES)) {
            throw LoggingUtility.logError(logger, new UnsupportedOperationException(
                "StandardCopyOption.COPY_ATTRIBUTES must be specified as the service will always copy "
                    + "file attributes."));
        }
        optionsList.remove(StandardCopyOption.COPY_ATTRIBUTES);
        if (optionsList.contains(StandardCopyOption.REPLACE_EXISTING)) {
            replaceExisting = true;
            optionsList.remove(StandardCopyOption.REPLACE_EXISTING);
        }
        if (!optionsList.isEmpty()) {
            throw LoggingUtility.logError(logger, new UnsupportedOperationException("Unsupported copy option found. "
                + "Only StandardCopyOption.COPY_ATTRIBUTES and StandareCopyOption.REPLACE_EXISTING are supported."));
        }

        // Validate paths. Build resources.
        // Copying a root directory or attempting to create/overwrite a root directory is illegal.
        AzureResource sourceRes = new AzureResource(source);
        AzurePath.ensureFileSystemOpen(sourceRes.getPath());
        AzureResource destinationRes = new AzureResource(destination);
        AzurePath.ensureFileSystemOpen(destinationRes.getPath());

        // Check destination is not a directory with children.
        DirectoryStatus destinationStatus = destinationRes.checkDirStatus();
        if (destinationStatus.equals(DirectoryStatus.NOT_EMPTY)) {
            throw LoggingUtility.logError(logger, new DirectoryNotEmptyException(destination.toString()));
        }

        /*
        Set request conditions if we should not overwrite. We can error out here if we know something already exists,
        but we will also create request conditions as a safeguard against overwriting something that was created
        between our check and put.
         */
        BlobRequestConditions requestConditions = null;
        if (!replaceExisting) {
            if (!destinationStatus.equals(DirectoryStatus.DOES_NOT_EXIST)) {
                throw LoggingUtility.logError(logger,
                    new FileAlreadyExistsException(destinationRes.getPath().toString()));
            }
            requestConditions = new BlobRequestConditions().setIfNoneMatch("*");
        }

        /*
        More path validation

        Check that the parent for the destination exists. We only need to perform this check if there is nothing
        currently at the destination, for if the destination exists, its parent at least weakly exists and we
        can skip a service call.
         */
        if (destinationStatus.equals(DirectoryStatus.DOES_NOT_EXIST) && !destinationRes.checkParentDirectoryExists()) {
            throw LoggingUtility.logError(logger, new IOException("Parent directory of destination location does not "
                + "exist. The destination path is therefore invalid. Destination: "
                + destinationRes.getPath().toString()));
        }

        /*
        Try to copy the resource at the source path.

        There is an optimization here where we try to do the copy first and only check for a virtual directory if
        there's a 404. In the cases of files and concrete directories, this only requires one request. For virtual
        directories, however, this requires three requests: failed copy, check status, create directory. Depending on
        customer scenarios and how many virtual directories they copy, it could be better to check the directory status
        first and then do a copy or createDir, which would always be two requests for all resource types.
         */
        try {
            SyncPoller<BlobCopyInfo, Void> pollResponse =
                destinationRes.getBlobClient().beginCopy(sourceRes.getBlobClient().getBlobUrl(), null, null, null,
                    null, requestConditions, null);
            pollResponse.waitForCompletion(Duration.ofSeconds(COPY_TIMEOUT_SECONDS));
        } catch (BlobStorageException e) {
            // If the source was not found, it could be because it's a virtual directory. Check the status.
            // If a non-dir resource existed, it would have been copied above. This check is therefore sufficient.
            if (e.getErrorCode().equals(BlobErrorCode.BLOB_NOT_FOUND)
                && !sourceRes.checkDirStatus().equals(DirectoryStatus.DOES_NOT_EXIST)) {
                /*
                We already checked that the parent exists and validated the paths above, so we can put the blob
                directly.
                 */
                destinationRes.putDirectoryBlob(requestConditions);
            } else {
                throw LoggingUtility.logError(logger, new IOException(e));
            }
        } catch (RuntimeException e) { // To better log possible timeout from poller.
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    // Used for checking the status of the root directory. To be implemented later when needed.
    /*int checkRootDirStatus(BlobContainerClient rootClient) {

    }*/

    /**
     * Unsupported.
     *
     * @param path path
     * @param path1 path
     * @param copyOptions options
     * @throws UnsupportedOperationException Operation is not supported.
     */
    @Override
    public void move(Path path, Path path1, CopyOption... copyOptions) throws IOException {
        throw LoggingUtility.logError(logger, new UnsupportedOperationException());
    }

    /**
     * Unsupported.
     *
     * @param path path
     * @param path1 path
     * @throws UnsupportedOperationException Operation is not supported.
     */
    @Override
    public boolean isSameFile(Path path, Path path1) throws IOException {
        throw LoggingUtility.logError(logger, new UnsupportedOperationException());
    }

    /**
     * Always returns false as hidden files are not supported.
     *
     * @param path the path
     * @return false
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public boolean isHidden(Path path) throws IOException {
        return false;
    }

    /**
     * Unsupported.
     *
     * @param path path
     * @return the file store where the file is stored.
     * @throws UnsupportedOperationException Operation is not supported.
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public FileStore getFileStore(Path path) throws IOException {
        throw LoggingUtility.logError(logger, new UnsupportedOperationException());
    }

    /**
     * Checks the existence, and optionally the accessibility, of a file.
     * <p>
     * This method may only be used to check the existence of a file. It is not possible to determine the permissions
     * granted to a given client, so if any mode argument is specified, an {@link UnsupportedOperationException} will be
     * thrown.
     *
     * @param path the path to the file to check
     * @param accessModes The access modes to check; may have zero elements
     * @throws NoSuchFileException if a file does not exist
     * @throws java.nio.file.AccessDeniedException the requested access would be denied or the access cannot be
     * determined because the Java virtual machine has insufficient privileges or other reasons
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public void checkAccess(Path path, AccessMode... accessModes) throws IOException {
        if (accessModes != null && accessModes.length != 0) {
            throw LoggingUtility.logError(logger, new AccessDeniedException("The access cannot be determined."));
        }
        AzurePath.ensureFileSystemOpen(path);

        /*
        Some static utility methods in the jdk require checking access on a root. ReadAttributes is not supported on
        roots as they are containers. Furthermore, we always assume that roots exist as they are verified at creation
        and cannot be deleted by the file system. Thus, we prefer a short circuit for roots.
         */
        if (path instanceof AzurePath && ((AzurePath) path).isRoot()) {
            return;
        }

        // Read attributes already wraps BlobStorageException in an IOException.
        try {
            readAttributes(path, BasicFileAttributes.class);
        } catch (IOException e) {
            Throwable cause = e.getCause();
            if (cause instanceof BlobStorageException
                && BlobErrorCode.BLOB_NOT_FOUND.equals(((BlobStorageException) cause).getErrorCode())) {
                throw LoggingUtility.logError(logger, new NoSuchFileException(path.toString()));
            } else {
                throw LoggingUtility.logError(logger, e);
            }
        }
    }

    /**
     * Returns a file attribute view of a given type.
     * <p>
     * See {@link AzureBasicFileAttributeView} and {@link AzureBlobFileAttributeView} for more information.
     * <p>
     * Reading or setting attributes on a virtual directory is not supported and will throw an {@link IOException}. See
     * {@link #createDirectory(Path, FileAttribute[])} for more information on virtual directories.
     *
     * @param path the path to the file
     * @param type the Class object corresponding to the file attribute view
     * @param linkOptions ignored
     * @return a file attribute view of the specified type, or null if the attribute view type is not available
     */
    @Override
    @SuppressWarnings("unchecked")
    public <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... linkOptions) {
        /*
        No resource validation is necessary here. That can happen at the time of making a network requests internal to
        the view object.
         */
        if (type == BasicFileAttributeView.class || type == AzureBasicFileAttributeView.class) {
            return (V) new AzureBasicFileAttributeView(path);
        } else if (type == AzureBlobFileAttributeView.class) {
            return (V) new AzureBlobFileAttributeView(path);
        } else {
            return null;
        }
    }

    /**
     * Reads a file's attributes as a bulk operation.
     * <p>
     * See {@link AzureBasicFileAttributes} and {@link AzureBlobFileAttributes} for more information.
     * <p>
     * Reading attributes on a virtual directory is not supported and will throw an {@link IOException}. See
     * {@link #createDirectory(Path, FileAttribute[])} for more information on virtual directories.
     *
     * @param path the path to the file
     * @param type the Class of the file attributes required to read
     * @param linkOptions ignored
     * @return the file attributes
     * @throws UnsupportedOperationException if an attributes of the given type are not supported
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    @SuppressWarnings("unchecked")
    public <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... linkOptions)
        throws IOException {
        AzurePath.ensureFileSystemOpen(path);

        Class<? extends BasicFileAttributeView> view;
        if (type == BasicFileAttributes.class || type == AzureBasicFileAttributes.class) {
            view = AzureBasicFileAttributeView.class;
        } else if (type == AzureBlobFileAttributes.class) {
            view = AzureBlobFileAttributeView.class;
        } else {
            throw LoggingUtility.logError(logger, new UnsupportedOperationException());
        }

        /*
        Resource validation will happen in readAttributes of the view. We don't want to double check, and checking
        internal to the view ensures it is always checked no matter which code path is taken.
         */
        return (A) getFileAttributeView(path, view, linkOptions).readAttributes();
    }

    /**
     * Reads a set of file attributes as a bulk operation.
     * <p>
     * See {@link AzureBasicFileAttributes} and {@link AzureBlobFileAttributes} for more information.
     * <p>
     * Reading attributes on a virtual directory is not supported and will throw an {@link IOException}. See
     * {@link #createDirectory(Path, FileAttribute[])} for more information on virtual directories.
     *
     * @param path the path to the file
     * @param attributes the attributes to read
     * @param linkOptions ignored
     * @return a map of the attributes returned; may be empty. The map's keys are the attribute names, its values are
     * the attribute values
     * @throws UnsupportedOperationException if an attributes of the given type are not supported
     * @throws IllegalArgumentException if no attributes are specified or an unrecognized attributes is specified
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public Map<String, Object> readAttributes(Path path, String attributes, LinkOption... linkOptions)
        throws IOException {
        if (attributes == null) {
            throw LoggingUtility.logError(logger, new IllegalArgumentException("Attribute string cannot be null."));
        }

        AzurePath.ensureFileSystemOpen(path);

        Map<String, Object> results = new HashMap<>();

        /*
        AzureBlobFileAttributes can do everything the basic attributes can do and more. There's no need to instantiate
        one of each if both are specified somewhere in the list as that will waste a network call. This can be
        generified later if we need to add more attribute types, but for now we can stick to just caching the supplier
        for a single attributes object.
         */
        Map<String, Supplier<Object>> attributeSuppliers = null; // Initialized later as needed.
        String viewType;
        String attributeList;
        String[] parts = attributes.split(":");

        if (parts.length > 2) {
            throw LoggingUtility.logError(logger,
                new IllegalArgumentException("Invalid format for attribute string: " + attributes));
        }

        if (parts.length == 1) {
            viewType = "basic"; // Per jdk docs.
            attributeList = attributes;
        } else {
            viewType = parts[0];
            attributeList = parts[1];
        }

        /*
        For specificity, our basic implementation of BasicFileAttributes uses the name azureBasic. However, the docs
        state that "basic" must be supported, so we funnel to azureBasic.
         */
        if (viewType.equals("basic")) {
            viewType = AzureBasicFileAttributeView.NAME;
        }
        if (!viewType.equals(AzureBasicFileAttributeView.NAME) && !viewType.equals(AzureBlobFileAttributeView.NAME)) {
            throw LoggingUtility.logError(logger,
                new UnsupportedOperationException("Invalid attribute view: " + viewType));
        }

        for (String attributeName : attributeList.split(",")) {
            /*
            We rely on the azureBlobFAV to actually do the work here as mentioned above, but if basic is specified, we
            should at least validate that the attribute is available on a basic view.
             */
            // TODO: Put these strings in constants
            if (viewType.equals(AzureBasicFileAttributeView.NAME)) {
                if (!AzureBasicFileAttributes.ATTRIBUTE_STRINGS.contains(attributeName) && !attributeName.equals("*")) {
                    throw LoggingUtility.logError(logger,
                        new IllegalArgumentException("Invalid attribute. View: " + viewType
                            + ". Attribute: " + attributeName));
                }
            }

            // As mentioned, azure blob can fulfill requests to both kinds of views.
            // Populate the supplier if we haven't already.
            if (attributeSuppliers == null) {
                attributeSuppliers = AzureBlobFileAttributes.getAttributeSuppliers(
                    this.readAttributes(path, AzureBlobFileAttributes.class, linkOptions));
            }

            // If "*" is specified, add all of the attributes from the specified set.
            if (attributeName.equals("*")) {
                Set<String> attributesToAdd;
                if (viewType.equals(AzureBasicFileAttributeView.NAME)) {
                    for (String attr : AzureBasicFileAttributes.ATTRIBUTE_STRINGS) {
                        results.put(attr, attributeSuppliers.get(attr).get());
                    }
                } else {
                    // attributeSuppliers is guaranteed to have been set by this point.
                    for (Map.Entry<String, Supplier<Object>> entry: attributeSuppliers.entrySet()) {
                        results.put(entry.getKey(), entry.getValue().get());
                    }
                }

            } else if (!attributeSuppliers.containsKey(attributeName)) {
                // Validate that the attribute is legal and add the value returned by the supplier to the results.
                throw LoggingUtility.logError(logger,
                    new IllegalArgumentException("Invalid attribute. View: " + viewType
                        + ". Attribute: " + attributeName));
            } else {
                results.put(attributeName, attributeSuppliers.get(attributeName).get());

            }
        }

        // Throw if nothing specified per jdk docs.
        if (results.isEmpty()) {
            throw LoggingUtility.logError(logger,
                new IllegalArgumentException("No attributes were specified. Attributes: " + attributes));
        }

        return results;
    }

    /**
     * Sets the value of a file attribute.
     * <p>
     * See {@link AzureBlobFileAttributeView} for more information.
     * <p>
     * Setting attributes on a virtual directory is not supported and will throw an {@link IOException}. See
     * {@link #createDirectory(Path, FileAttribute[])} for more information on virtual directories.
     *
     * @param path the path to the file
     * @param attributes the attribute to set
     * @param value the attribute value
     * @param linkOptions ignored
     * @throws UnsupportedOperationException if an attribute view is not available
     * @throws IllegalArgumentException if the attribute name is not specified, or is not recognized, or the attribute
     * value is of the correct type but has an inappropriate value
     * @throws ClassCastException If the attribute value is not of the expected type or is a collection containing
     * elements that are not of the expected type
     * @throws IOException If an I/O error occurs.
     * @throws SecurityException never
     */
    @Override
    public void setAttribute(Path path, String attributes, Object value, LinkOption... linkOptions) throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        String viewType;
        String attributeName;
        String[] parts = attributes.split(":");
        if (parts.length > 2) {
            throw LoggingUtility.logError(logger,
                new IllegalArgumentException("Invalid format for attribute string: " + attributes));
        }
        if (parts.length == 1) {
            viewType = "basic"; // Per jdk docs.
            attributeName = attributes;
        } else {
            viewType = parts[0];
            attributeName = parts[1];
        }

        /*
        For specificity, our basic implementation of BasicFileAttributes uses the name azureBasic. However, the docs
        state that "basic" must be supported, so we funnel to azureBasic.
         */
        if (viewType.equals("basic")) {
            viewType = AzureBasicFileAttributeView.NAME;
        }

        // We don't actually support any setters on the basic view.
        if (viewType.equals(AzureBasicFileAttributeView.NAME)) {
            throw LoggingUtility.logError(logger,
                new IllegalArgumentException("Invalid attribute. View: " + viewType
                    + ". Attribute: " + attributeName));
        } else if (viewType.equals(AzureBlobFileAttributeView.NAME)) {
            Map<String, Consumer<Object>> attributeConsumers = AzureBlobFileAttributeView.setAttributeConsumers(
                this.getFileAttributeView(path, AzureBlobFileAttributeView.class, linkOptions));
            if (!attributeConsumers.containsKey(attributeName)) {
                // Validate that the attribute is legal and add the value returned by the supplier to the results.
                throw LoggingUtility.logError(logger,
                    new IllegalArgumentException("Invalid attribute. View: " + viewType
                        + ". Attribute: " + attributeName));
            }
            try {
                attributeConsumers.get(attributeName).accept(value);
            } catch (UncheckedIOException e) {
                if (e.getMessage().equals(AzureBlobFileAttributeView.ATTR_CONSUMER_ERROR)) {
                    throw LoggingUtility.logError(logger, e.getCause());
                }
            }
        } else {
            throw LoggingUtility.logError(logger,
                new UnsupportedOperationException("Invalid attribute view: " + viewType));
        }
    }

    void closeFileSystem(String fileSystemName) {
        this.openFileSystems.remove(fileSystemName);
    }

    private String extractAccountName(URI uri) {
        if (!uri.getScheme().equals(this.getScheme())) {
            throw LoggingUtility.logError(this.logger, new IllegalArgumentException(
                "URI scheme does not match this provider"));
        }
        if (CoreUtils.isNullOrEmpty(uri.getQuery())) {
            throw LoggingUtility.logError(this.logger, new IllegalArgumentException("URI does not contain a query "
                + "component. FileSystems require a URI of the format \"azb://?account=<account_name>\"."));
        }

        String accountName = Flux.fromArray(uri.getQuery().split("&"))
                .filter(s -> s.startsWith(ACCOUNT_QUERY_KEY + "="))
                .switchIfEmpty(Mono.error(LoggingUtility.logError(this.logger, new IllegalArgumentException(
                        "URI does not contain an \"" + ACCOUNT_QUERY_KEY + "=\" parameter. FileSystems require a URI "
                            + "of the format \"azb://?account=<account_name>\""))))
                .map(s -> s.substring(ACCOUNT_QUERY_KEY.length() + 1))
                .blockLast();

        if (CoreUtils.isNullOrEmpty(accountName)) {
            throw LoggingUtility.logError(logger, new IllegalArgumentException("No account name provided in URI"
                + " query."));
        }

        return accountName;
    }
}
