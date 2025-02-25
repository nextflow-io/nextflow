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

import java.nio.channels.Channels
import java.nio.channels.SeekableByteChannel
import java.nio.file.DirectoryNotEmptyException
import java.nio.file.DirectoryStream
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.WatchService
import java.nio.file.attribute.UserPrincipalLookupService
import java.time.Duration
import java.time.temporal.ChronoUnit
import java.util.concurrent.TimeoutException
import java.util.function.Predicate

import com.azure.core.util.polling.SyncPoller
import com.azure.storage.blob.BlobServiceClient
import com.azure.storage.blob.models.BlobContainerItem
import com.azure.storage.blob.models.BlobCopyInfo
import com.azure.storage.blob.models.BlobItem
import com.azure.storage.blob.models.BlobStorageException
import com.azure.storage.blob.models.ListBlobsOptions
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.azure.config.AzConfig
/**
 * Implements a file system for Azure Blob Storage service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzFileSystem extends FileSystem {

    static class GuessPath {
        boolean exists
        boolean directory
        boolean empty
    }

    private static String EMPTY_DIR_MARKER = '.azure_blob_dir'

    private static String SLASH = '/'

    private AzFileSystemProvider provider

    private BlobServiceClient storageClient

    private String containerName

    private int maxCopyDurationSecs = 3600


    @PackageScope AzFileSystem() {}

    @PackageScope
    AzFileSystem(AzFileSystemProvider provider, BlobServiceClient storageClient, String bucket) {
        this.provider = provider
        this.containerName = bucket
        this.storageClient = storageClient
    }

    String getContainerName() { containerName }

    BlobServiceClient getBlobServiceClient() { storageClient }

    @Override
    AzFileSystemProvider provider() {
        return provider
    }

    @Override
    void close() throws IOException {
        // nothing to do
    }

    @Override
    boolean isOpen() {
        return true
    }

    @Override
    boolean isReadOnly() {
        return containerName == SLASH
    }

    @Override
    String getSeparator() {
        return SLASH
    }

    Iterable<? extends Path> getRootDirectories() {
        return containerName == SLASH ? listContainers() : [new AzPath(this, "/$containerName/") ]
    }

    private Iterable<? extends Path> listContainers() {
        return apply(()-> listContainers0())
    }

    private Iterable<? extends Path> listContainers0() {
        final containers = new ArrayList()
        storageClient
                .listBlobContainers()
                .forEach { BlobContainerItem it -> provider.getPath(it.getName()) }
        return containers
    }

    @Override
    Iterable<FileStore> getFileStores() {
        throw new UnsupportedOperationException("Operation 'getFileStores' is not supported by AzFileSystem")
    }

    @Override
    Set<String> supportedFileAttributeViews() {
        return Collections.unmodifiableSet( ['basic'] as Set )
    }

    /**
     * Get a {@link AzPath} given a string path. When a path starts with `/`
     * is interpreted as an absolute path in which the first component
     * is the bucket Azure storage bucket name.
     *
     * Otherwise it is interpreted as relative object name in the current
     * file system.
     *
     */
    @Override
    AzPath getPath(String path, String... more) {
        assert path

        if( more ) {
            path = concat(path,more)
        }

        if( path.startsWith('/') ) {
            return provider.getPath(path.substring(1))
        }
        else {
            return new AzPath(this, path)
        }
    }

    private String concat(String path, String... more) {
        def concat = []
        while( path.length()>1 && path.endsWith('/') )
            path = path.substring(0,path.length()-1)
        concat << path
        concat.addAll( more.collect {String it -> trimSlash(it)} )
        return concat.join('/')
    }

    private String trimSlash(String str) {
        while( str.startsWith('/') )
            str = str.substring(1)
        while( str.endsWith('/') )
            str = str.substring(0,str.length()-1)
        return str
    }

    @Override
    PathMatcher getPathMatcher(String syntaxAndPattern) {
        throw new UnsupportedOperationException("Operation 'getPathMatcher' is not supported by AzFileSystem")
    }

    @Override
    UserPrincipalLookupService getUserPrincipalLookupService() {
        throw new UnsupportedOperationException("Operation 'getUserPrincipalLookupService' is not supported by AzFileSystem")
    }

    @Override
    WatchService newWatchService() throws IOException {
        throw new UnsupportedOperationException("Operation 'newWatchService' is not supported by AzFileSystem")
    }

    @PackageScope
    SeekableByteChannel newReadableByteChannel(AzPath path) {
        if( path.isContainer() )
            throw new IllegalArgumentException("Operation not allowed on blob container path: ${path.toUriString()}")
        final client = path.blobClient()

        try {
            final size = client.getProperties().getBlobSize()
            final channel = Channels.newChannel( client.openInputStream() )
            return new AzReadableByteChannel(channel, size)
        }
        catch (BlobStorageException e) {
            if( e.statusCode == 404 )
                throw new NoSuchFileException(path.toUriString())
            else
                throw new IOException("Unexpected exception processing Azure container blob: ${path.toUriString()}", e)
        }
    }

    @PackageScope
    SeekableByteChannel newWritableByteChannel(AzPath path) {
        if( path.isContainer() )
            throw new IllegalArgumentException("Operation not allowed on blob container path: ${path.toUriString()}")

        try {
            final client = path.blobClient()
            final outStream = client.getBlockBlobClient().getBlobOutputStream(true)
            final writer = Channels.newChannel(outStream)
            return new AzWriteableByteChannel(writer)
        }
        catch (BlobStorageException e) {
            if( e.statusCode == 404 )
                throw new NoSuchFileException(path.toUriString())
            else
                throw new IOException("Unexpected exception processing Azure container blob: ${path.toUriString()}", e)
        }
    }

    @PackageScope
    DirectoryStream<Path> newDirectoryStream(AzPath dir, DirectoryStream.Filter<? super Path> filter) {
        if( !dir.containerName )
            throw new NoSuchFileException("Missing Azure storage container name: ${dir.toUriString()}")

        if( dir.containerName == SLASH ) {
            return listContainers(dir, filter)
        }
        else {
            listFiles(dir, filter)
        }
    }


    /**
     * Create a {@link DirectoryStream} that allows the iteration over the files in a Azure Blob Storage
     * file system system
     *
     * @param fs The underlying {@link AzFileSystem} instance
     * @param filter A {@link java.nio.file.DirectoryStream.Filter} object to select which files to include in the file traversal
     * @return A {@link DirectoryStream} object to traverse the associated objects
     */
    private DirectoryStream<Path> listFiles(AzPath dir, DirectoryStream.Filter<? super Path> filter) {
        return apply(()-> listFiles0(dir,filter))
    }

    private DirectoryStream<Path> listFiles0(AzPath dir, DirectoryStream.Filter<? super Path> filter) {

        // -- create the list operation options
        def prefix = dir.blobName()
        if( prefix && !prefix.endsWith('/') ) {
            prefix += '/'
        }

        // -- list the bucket content
        final blobs = dir.containerClient()
                .listBlobsByHierarchy(prefix)
                .iterator()

        // wrap the result with a directory stream
        return new DirectoryStream<Path>() {

            @Override
            Iterator<Path> iterator() {
                return new AzPathIterator.ForBlobs(dir, blobs, filter)
            }

            @Override void close() throws IOException { }
        }
    }

    /**
     * Create a {@link DirectoryStream} that allows the iteration over the buckets in a Azure Blob Storage Storage
     * file system system
     *
     * @param fs The underlying {@link AzFileSystem} instance
     * @param filter A {@link java.nio.file.DirectoryStream.Filter} object to select which buckets to include in the file traversal
     * @return A {@link DirectoryStream} object to traverse the associated objects
     */
    private DirectoryStream<Path> listContainers(AzPath path, DirectoryStream.Filter<? super Path> filter ) {
        return apply(()-> listContainers0(path, filter))
    }

    private DirectoryStream<Path> listContainers0(AzPath path, DirectoryStream.Filter<? super Path> filter) {

        Iterator<BlobContainerItem> containers = storageClient.listBlobContainers().iterator()

        return new DirectoryStream<Path>() {

            @Override
            Iterator<Path> iterator() {
                return new AzPathIterator.ForContainers(path, containers, filter)
            }

            @Override void close() throws IOException { }
        }
    }

    @PackageScope
    void createDirectory(AzPath path) {
        if( isReadOnly() )
            throw new UnsupportedOperationException("Operation 'createDirectory' not supported in root path")

        if( !path.containerName )
            throw new IllegalArgumentException("Missing Azure storage blob container name")

        if( path.isContainer() ) {
            storageClient.createBlobContainer(path.checkContainerName())
        }
        else {
            path.directory = true
            // Create an hidden file blob to as a placeholder for a new empty directory
            // NOTE: this is *not* required by this implementation however third party tools
            // such as azcopy get confused creation a blob name ending with a slash
            // therefore this library creates an hidden file to simulate the creation
            // of a directory
            path.resolve(EMPTY_DIR_MARKER)
                    .blobClient()
                    .getAppendBlobClient()
                    .create()
        }
    }

    @PackageScope
    void delete(AzPath path)  {
        if( !path.containerName )
            throw new IllegalArgumentException("Missing Azure blob container name")

        if( path.isContainer() ) {
            deleteBucket(path)
        }
        else {
            deleteFile(path)
        }
    }

    private void checkContainerExistsOrEmpty(AzPath path) {
        try {
            final container = path.containerClient()
            final opts = new ListBlobsOptions().setMaxResultsPerPage(10)
            final blobs = apply(()-> container.listBlobs(opts, null))
            if( apply(()-> blobs.iterator().hasNext()) ) {
                throw new DirectoryNotEmptyException(path.toUriString())
            }
        }
        catch (BlobStorageException e) {
            if( e.statusCode == 404 )
                throw new NoSuchFileException(path.toUriString())
            else
                throw new IOException("Unexpected exception accessing blob storage path: ${path.toUriString()}")
        }
    }

    private void checkPathExistOrEmpty(AzPath path) {
        final result = guessPath(path)

        if( result.directory && !result.empty )
            throw new DirectoryNotEmptyException(path.toUriString())

        if( !result.exists )
            throw new NoSuchFileException(path.toUriString())
    }

    private GuessPath guessPath(AzPath path) {
        boolean exists = false
        boolean isDirectory = false

        final opts = new ListBlobsOptions()
                .setPrefix(path.blobName())
                .setMaxResultsPerPage(10)
        try {
            final values = apply(()-> path.containerClient().listBlobs(opts,null).iterator())

            final char SLASH = '/'
            final String name = path.blobName()

            int count=0
            while( apply(()-> values.hasNext()) ) {
                BlobItem blob = apply(()-> values.next())
                if( blob.name == name )
                    exists = true
                else if( blob.name.startsWith(name) && blob.name.charAt(name.length())==SLASH ) {
                    exists = true
                    isDirectory = true
                }
                // ignore empty dir marker file
                if( blob.name.endsWith('/'+EMPTY_DIR_MARKER))
                    continue
                count++
            }

            if( isDirectory && count>=1 )
                return new GuessPath(exists: true, directory: true, empty: false)
            else
                return new GuessPath(exists: exists, directory: isDirectory, empty: true)

        }
        catch (BlobStorageException e) {
            if( e.statusCode==404 )
                return new GuessPath(exists: false)
            throw  e
        }
    }

    private void deleteFile(AzPath path) {
        checkPathExistOrEmpty(path)
        apply(()-> path.blobClient().delete())
    }

    private void deleteBucket(AzPath path) {
        checkContainerExistsOrEmpty(path)
        apply(()-> path.containerClient().delete())
    }

    @PackageScope
    void copy(AzPath source, AzPath target) {
        final sasToken = provider.getSasToken()
        String sourceUrl = source.blobClient().getBlobUrl()

        if (sasToken != null) {
            if (sourceUrl.contains('?')){
                sourceUrl = String.format("%s&%s", sourceUrl, sasToken);
            } else {
                sourceUrl = String.format("%s?%s", sourceUrl, sasToken);
            }
        }

        SyncPoller<BlobCopyInfo, Void> pollResponse =
                target.blobClient().beginCopy( sourceUrl, null )
        pollResponse.waitForCompletion(Duration.ofSeconds(maxCopyDurationSecs))
    }

    @PackageScope
    AzFileAttributes readAttributes(AzPath path) {
        final cache = path.attributesCache()
        if( cache )
            return cache

        if( path.toString() == SLASH ) {
            return AzFileAttributes.root()
        }

        if( path.containerName && !path.blobName()) {
            return readContainerAttrs0(path)
        }
        else
            return readBlobAttrs0(path)
    }

    private AzFileAttributes readBlobAttrs0(AzPath path) {
        try {
            return new AzFileAttributes(path.blobClient())
        }
        catch (BlobStorageException e) {
            if( e.statusCode != 404 )
                throw e

            // if the previous check failed and
            // the following pass, it can only be a directory
            final result = guessPath(path)
            if( result.exists ) {
                def name = path.blobName()
                if( !name.endsWith('/') ) name += '/'
                return new AzFileAttributes(path.containerClient(), name)
            }
            else {
                return null
            }
        }
    }

    private AzFileAttributes readContainerAttrs0(AzPath path) {
        try {
            return new AzFileAttributes(path.containerClient())
        }
        catch (BlobStorageException e) {
            if( e.statusCode==404 )
                throw new NoSuchFileException(path.toUriString())
            throw e
        }
    }

    @PackageScope
    AzFileAttributesView getFileAttributeView(AzPath path) {
        try {
            return new AzFileAttributesView(path.blobClient())
        }
        catch (BlobStorageException e) {
            if( e.statusCode==404 )
                throw new NoSuchFileException(path.toUriString())
            else
                throw new IOException("Unable to get attributes for file: ${path.toUriString()}", e)
        }
    }

    @PackageScope
    boolean exists( AzPath path ) {
        try {
            return readAttributes(path) != null
        }
        catch( IOException e ){
            return false
        }
    }


    /**
     * Creates a retry policy using the configuration specified by {@link nextflow.cloud.azure.config.AzRetryConfig}
     *
     * @param cond A predicate that determines when a retry should be triggered
     * @return The {@link dev.failsafe.RetryPolicy} instance
     */
    protected <T> RetryPolicy<T> retryPolicy(Predicate<? extends Throwable> cond) {
        // this is needed because apparently bytebuddy used by testing framework is not able
        // to handle properly this method signature using both generics and `@Memoized` annotation.
        // therefore the `@Memoized` has been moved to the inner method invocation
        return (RetryPolicy<T>) retryPolicy0(cond)
    }

    @Memoized
    protected RetryPolicy retryPolicy0(Predicate<? extends Throwable> cond) {
        final cfg = AzConfig.getConfig().retryConfig()
        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                log.debug("Azure I/O exception - attempt: ${event.attemptCount}; cause: ${event.lastFailure?.message}")
            }
        }
        return RetryPolicy.builder()
                .handleIf(cond)
                .withBackoff(cfg.delay.toMillis(), cfg.maxDelay.toMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(cfg.maxAttempts)
                .withJitter(cfg.jitter)
                .onRetry(listener)
                .build()
    }

    /**
     * Carry out the invocation of the specified action using a retry policy
     * when {@code TooManyRequests} Azure Batch error is returned
     *
     * @param action A {@link dev.failsafe.function.CheckedSupplier} instance modeling the action to be performed in a safe manner
     * @return The result of the supplied action
     */
    protected <T> T apply(CheckedSupplier<T> action) {
        final policy = retryPolicy((Throwable t) -> t instanceof IOException || t.cause instanceof IOException || t instanceof TimeoutException || t.cause instanceof TimeoutException)
        return Failsafe.with(policy).get(action)
    }
}
