/*
 * Copyright 2021, Microsoft Corp
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
package nextflow.cloud.azure.nio

import static java.nio.file.StandardCopyOption.*
import static java.nio.file.StandardOpenOption.*

import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryStream
import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.FileSystemAlreadyExistsException
import java.nio.file.FileSystemNotFoundException
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributeView
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import com.azure.storage.blob.BlobServiceClient
import com.azure.storage.blob.models.BlobStorageException
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.azure.batch.AzHelper
/**
 * Implements NIO File system provider for Azure Blob Storage
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzFileSystemProvider extends FileSystemProvider {

    public static final String AZURE_STORAGE_ACCOUNT_NAME = 'AZURE_STORAGE_ACCOUNT_NAME'
    public static final String AZURE_STORAGE_ACCOUNT_KEY = 'AZURE_STORAGE_ACCOUNT_KEY'
    public static final String AZURE_STORAGE_SAS_TOKEN = 'AZURE_STORAGE_SAS_TOKEN'

    public static final String AZURE_CLIENT_ID = 'AZURE_CLIENT_ID'
    public static final String AZURE_CLIENT_SECRET = 'AZURE_CLIENT_SECRET'
    public static final String AZURE_TENANT_ID = 'AZURE_TENANT_ID'

    public static final String AZURE_MANAGED_IDENTITY_USER = 'AZURE_MANAGED_IDENTITY_USER'
    public static final String AZURE_MANAGED_IDENTITY_SYSTEM = 'AZURE_MANAGED_IDENTITY_SYSTEM'

    public static final String SCHEME = 'az'

    private Map<String,String> env = new HashMap<>(System.getenv())
    private Map<String,AzFileSystem> fileSystems = [:]
    private String sasToken = null
    private String accountKey = null

    /**
     * @inheritDoc
     */
    @Override
    String getScheme() {
        return SCHEME
    }

    String getSasToken() {
        return this.sasToken
    }

    String getAccountKey() {
        return this.accountKey
    }

    static private AzPath asAzPath(Path path) {
        if( path !instanceof AzPath )
            throw new IllegalArgumentException("Not a valid Azure blob storage path object: `$path` [${path?.class?.name?:'-'}]" )
        return (AzPath)path
    }

    protected String getContainerName(URI uri) {
        assert uri
        if( !uri.scheme )
            throw new IllegalArgumentException("Missing URI scheme")

        if( uri.scheme.toLowerCase() != SCHEME )
            throw new IllegalArgumentException("Mismatch provider URI scheme: `$scheme`")

        if( !uri.authority ) {
            if( uri.path == '/' )
                return '/'
            else
                throw new IllegalArgumentException("Missing Azure blob storage container name")
        }

        return uri.authority.toLowerCase()
    }

    protected BlobServiceClient createBlobServiceWithKey(String accountName, String accountKey) {
        AzHelper.getOrCreateBlobServiceWithKey(accountName, accountKey)
    }

    protected BlobServiceClient createBlobServiceWithToken(String accountName, String sasToken) {
        AzHelper.getOrCreateBlobServiceWithToken(accountName, sasToken)
    }

    protected BlobServiceClient createBlobServiceWithServicePrincipal(String accountName, String clientId, String clientSecret, String tenantId) {
        AzHelper.getOrCreateBlobServiceWithServicePrincipal(accountName, clientId, clientSecret, tenantId)
    }

    protected BlobServiceClient createBlobServiceWithManagedIdentity(String accountName, String clientId) {
        AzHelper.getOrCreateBlobServiceWithManagedIdentity(accountName, clientId)
    }

    /**
     * Constructs a new {@code FileSystem} object identified by a URI. This
     * method is invoked by the {@link java.nio.file.FileSystems#newFileSystem(URI, Map)}
     * method to open a new file system identified by a URI.
     *
     * <p> The {@code uri} parameter is an absolute, hierarchical URI, with a
     * scheme equal (without regard to case) to the scheme supported by this
     * provider. The exact form of the URI is highly provider dependent. The
     * {@code env} parameter is a map of provider specific properties to configure
     * the file system.
     *
     * <p> This method throws {@link java.nio.file.FileSystemAlreadyExistsException} if the
     * file system already exists because it was previously created by an
     * invocation of this method. Once a file system is {@link
     * java.nio.file.FileSystem#close closed} it is provider-dependent if the
     * provider allows a new file system to be created with the same URI as a
     * file system it previously created.
     *
     * @param   uri
     *          URI reference
     * @param   config
     *          A map of provider specific properties to configure the file system;
     *          may be empty
     *
     * @return  A new file system
     *
     * @throws  IllegalArgumentException
     *          If the pre-conditions for the {@code uri} parameter aren't met,
     *          or the {@code env} parameter does not contain properties required
     *          by the provider, or a property value is invalid
     * @throws  IOException
     *          An I/O error occurs creating the file system
     * @throws  SecurityException
     *          If a security manager is installed and it denies an unspecified
     *          permission required by the file system provider implementation
     * @throws  java.nio.file.FileSystemAlreadyExistsException
     *          If the file system has already been created
     */
    @Override
    AzFileSystem newFileSystem(URI uri, Map<String, ?> config) throws IOException {
        final bucket = getContainerName(uri)
        newFileSystem0(bucket, config)
    }

    /**
     * Creates a new {@link AzFileSystem} for the given `bucket`.
     *
     * @param bucket The bucket name for which the file system will be created
     * @param config
     *          A {@link Map} object holding the file system configuration settings. Valid keys:
     *          - credentials: path of the file
     *          - projectId
     *          - location
     *          - storageClass
     * @return
     * @throws IOException
     */
    synchronized AzFileSystem newFileSystem0(String bucket, Map<String, ?> config) throws IOException {

        if( fileSystems.containsKey(bucket) )
            throw new FileSystemAlreadyExistsException("File system already exists for Azure blob container: `$bucket`")

        final accountName = config.get(AZURE_STORAGE_ACCOUNT_NAME) as String
        final accountKey = config.get(AZURE_STORAGE_ACCOUNT_KEY) as String
        final sasToken = config.get(AZURE_STORAGE_SAS_TOKEN) as String

        final servicePrincipalId = config.get(AZURE_CLIENT_ID) as String
        final servicePrincipalSecret = config.get(AZURE_CLIENT_SECRET) as String
        final tenantId = config.get(AZURE_TENANT_ID) as String

        final managedIdentityUser = config.get(AZURE_MANAGED_IDENTITY_USER) as String
        final managedIdentitySystem = config.get(AZURE_MANAGED_IDENTITY_SYSTEM) as Boolean

        if( !accountName )
            throw new IllegalArgumentException("Missing AZURE_STORAGE_ACCOUNT_NAME")

        BlobServiceClient client

        if( managedIdentityUser || managedIdentitySystem ) {
            client = createBlobServiceWithManagedIdentity(accountName, managedIdentityUser)
        }
        else if( servicePrincipalSecret && servicePrincipalId && tenantId ) {
            client = createBlobServiceWithServicePrincipal(accountName, servicePrincipalId, servicePrincipalSecret, tenantId)
        }
        else if( sasToken ) {
            client = createBlobServiceWithToken(accountName, sasToken)
            this.sasToken = sasToken
        }
        else if( accountKey ) {
            client = createBlobServiceWithKey(accountName, accountKey)
            this.accountKey = accountKey
        }
        else {
            throw new IllegalArgumentException("Missing Azure storage credentials: please specify a managed identity, service principal, or storage account key")
        }

        final result = createFileSystem(client, bucket, config)
        fileSystems[bucket] = result
        return result
    }

    /**
     * Creates a new {@link AzFileSystem} object.
     *
     * @param client
     * @param bucket
     * @param config
     * @return
     */
    protected AzFileSystem createFileSystem(BlobServiceClient client, String bucket, Map<String,?> config) {
        def result = new AzFileSystem(this, client, bucket)
        return result
    }

    /**
     * Returns an existing {@code FileSystem} created by this provider.
     *
     * <p> This method returns a reference to a {@code FileSystem} that was
     * created by invoking the {@link #newFileSystem(URI,Map) newFileSystem(URI,Map)}
     * method. File systems created the {@link #newFileSystem(Path,Map)
     * newFileSystem(Path,Map)} method are not returned by this method.
     * The file system is identified by its {@code URI}. Its exact form
     * is highly provider dependent. In the case of the default provider the URI's
     * path component is {@code "/"} and the authority, query and fragment components
     * are undefined (Undefined components are represented by {@code null}).
     *
     * <p> Once a file system created by this provider is {@link
     * java.nio.file.FileSystem#close closed} it is provider-dependent if this
     * method returns a reference to the closed file system or throws {@link
     * java.nio.file.FileSystemNotFoundException}. If the provider allows a new file system to
     * be created with the same URI as a file system it previously created then
     * this method throws the exception if invoked after the file system is
     * closed (and before a new instance is created by the {@link #newFileSystem
     * newFileSystem} method).
     *
     * @param   uri
     *          URI reference
     *
     * @return  The file system
     *
     * @throws  IllegalArgumentException
     *          If the pre-conditions for the {@code uri} parameter aren't met
     * @throws  java.nio.file.FileSystemNotFoundException
     *          If the file system does not exist
     * @throws  SecurityException
     *          If a security manager is installed and it denies an unspecified
     *          permission.
     */
    @Override
    FileSystem getFileSystem(URI uri) {
        final bucket = getContainerName(uri)
        getFileSystem0(bucket,false)
    }

    protected AzFileSystem getFileSystem0(String bucket, boolean canCreate) {

        def fs = fileSystems.get(bucket)
        if( !fs ) {
            if( canCreate )
                fs = newFileSystem0(bucket, env)
            else
                throw new FileSystemNotFoundException("Missing Azure storage blob file system for bucket: `$bucket`")
        }

        return fs
    }

    /**
     * Return a {@code Path} object by converting the given {@link URI}. The
     * resulting {@code Path} is associated with a {@link FileSystem} that
     * already exists or is constructed automatically.
     *
     * <p> The exact form of the URI is file system provider dependent. In the
     * case of the default provider, the URI scheme is {@code "file"} and the
     * given URI has a non-empty path component, and undefined query, and
     * fragment components. The resulting {@code Path} is associated with the
     * default {@link java.nio.file.FileSystems#getDefault default} {@code FileSystem}.
     *
     * @param   uri
     *          The URI to convert
     *
     * @return  The resulting {@code Path}
     *
     * @throws  IllegalArgumentException
     *          If the URI scheme does not identify this provider or other
     *          preconditions on the uri parameter do not hold
     * @throws  java.nio.file.FileSystemNotFoundException
     *          The file system, identified by the URI, does not exist and
     *          cannot be created automatically
     * @throws  SecurityException
     *          If a security manager is installed and it denies an unspecified
     *          permission.
     */
    @Override
    AzPath getPath(URI uri) {
        final bucket = getContainerName(uri)
        bucket=='/' ? getPath('/') : getPath("$bucket/${uri.path}")
    }

    /**
     * Get a {@link AzPath} from an object path string
     *
     * @param path A path in the form {@code containerName/blobName}
     * @return A {@link AzPath} object
     */
    AzPath getPath(String path) {
        assert path

        // -- special root bucket
        if( path == '/' ) {
            final fs = getFileSystem0('/',true)
            return new AzPath(fs, "/")
        }

        // -- remove first slash, if any
        while( path.startsWith("/") )
            path = path.substring(1)

        // -- find the first component ie. the container name
        int p = path.indexOf('/')
        final bucket = p==-1 ? path : path.substring(0,p)

        // -- get the file system
        final fs = getFileSystem0(bucket,true)

        // create a new path
        new AzPath(fs, "/$path")
    }

    static private FileSystemProvider provider( Path path ) {
        path.getFileSystem().provider()
    }

    @Deprecated
    static private BlobServiceClient storage( Path path ) {
        ((AzPath)path).getFileSystem().getBlobServiceClient()
    }

    private void checkRoot(Path path) {
        if( path.toString() == '/' )
            throw new UnsupportedOperationException("Operation 'checkRoot' not supported on root path")
    }

    @Override
    SeekableByteChannel newByteChannel(Path obj, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        checkRoot(obj)

        final modeWrite = options.contains(WRITE) || options.contains(APPEND)
        final modeRead = options.contains(READ) || !modeWrite

        if( modeRead && modeWrite ) {
            throw new IllegalArgumentException("Azure Blob Storage file cannot be opened in R/W mode at the same time")
        }
        if( options.contains(APPEND) ) {
            throw new IllegalArgumentException("Azure Blob Storage file system does not support `APPEND` mode")
        }
        if( options.contains(SYNC) ) {
            throw new IllegalArgumentException("Azure Blob Storage file system does not support `SYNC` mode")
        }
        if( options.contains(DSYNC) ) {
            throw new IllegalArgumentException("Azure Blob Storage file system does not support `DSYNC` mode")
        }

        final path = asAzPath(obj)
        final fs = path.getFileSystem()
        if( modeRead ) {
            return fs.newReadableByteChannel(path)
        }

        // -- mode write
        if( options.contains(CREATE_NEW) ) {
            if( fs.exists(path) )
                throw new FileAlreadyExistsException(path.toUriString())
        }
        else if( !options.contains(CREATE)  ) {
            if( !fs.exists(path) )
                throw new NoSuchFileException(path.toUriString())
        }
        if( options.contains(APPEND) ) {
            throw new IllegalArgumentException("File APPEND mode is not supported by Azure Blob Storage")
        }
        return fs.newWritableByteChannel(path)
    }


    @Override
    DirectoryStream<Path> newDirectoryStream(Path obj, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final path = asAzPath(obj)
        path.fileSystem.newDirectoryStream(path, filter)
    }

    @Override
    void createDirectory(Path dir, FileAttribute<?>... attrs) throws IOException {
        checkRoot(dir)
        final path = asAzPath(dir)
        try {
            path.fileSystem.createDirectory(path)
        }
        catch (BlobStorageException e) {
            // 409 (CONFLICT) is returned when the path already
            // exists, ignore it
            if( e.statusCode!=409 )
                throw new IOException("Unable to create Azure blob directory: ${dir.toUriString()} - cause: ${e.message}", e)
        }
    }

    @Override
    void delete(Path obj) throws IOException {
        checkRoot(obj)
        final path = asAzPath(obj)
        path.fileSystem.delete(path)
    }


    @Override
    void copy(Path from, Path to, CopyOption... options) throws IOException {
        assert provider(from) == provider(to)
        if( from == to )
            return // nothing to do -- just return

        checkRoot(from); checkRoot(to)
        final source = asAzPath(from)
        final target = asAzPath(to)
        final fs = source.getFileSystem()

        if( options.contains(REPLACE_EXISTING) && fs.exists(target) ) {
            delete(target)
        }

        fs.copy(source, target)
    }

    @Override
    void move(Path source, Path target, CopyOption... options) throws IOException {
        copy(source,target,options)
        delete(source)
    }

    @Override
    boolean isSameFile(Path path, Path path2) throws IOException {
        return path == path2
    }

    @Override
    boolean isHidden(Path path) throws IOException {
        return path.getFileName()?.toString()?.startsWith('.')
    }

    @Override
    FileStore getFileStore(Path path) throws IOException {
        throw new UnsupportedOperationException("Operation 'getFileStore' is not supported by AzFileSystem")
    }

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        checkRoot(path)
        final az = asAzPath(path)
        readAttributes(az, AzFileAttributes.class)
        if( AccessMode.EXECUTE in modes)
            throw new AccessDeniedException(az.toUriString(), null, 'Execute permission not allowed')
    }

    @Override
    def <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        checkRoot(path)
        if( type == BasicFileAttributeView || type == AzFileAttributesView ) {
            def azPath = asAzPath(path)
            return (V)azPath.fileSystem.getFileAttributeView(azPath)
        }
        throw new UnsupportedOperationException("Not a valid Azure Blob Storage file attribute view: $type")
    }

    @Override
    def <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        if( type == BasicFileAttributes || type == AzFileAttributes ) {
            def azPath = asAzPath(path)
            def result = (A)azPath.fileSystem.readAttributes(azPath)
            if( result )
                return result
            throw new NoSuchFileException(azPath.toUriString())
        }
        throw new UnsupportedOperationException("Not a valid Azure Blob Storage file attribute type: $type")
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Operation 'readAttributes' is not supported by AzFileSystem")
    }

    @Override
    void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Operation 'setAttribute' is not supported by AzFileSystem")
    }

}
