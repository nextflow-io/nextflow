/*
 * Copyright 2013-2025, Seqera Labs
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

package io.seqera.tower.plugin.fs;

import io.seqera.http.HxClient;
import io.seqera.tower.plugin.TowerConfig;
import io.seqera.tower.plugin.TowerHxClientFactory;
import io.seqera.tower.plugin.datalink.DataLink;
import io.seqera.tower.plugin.datalink.DataLinkUtils;
import nextflow.SysEnv;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.net.http.HttpClient;
import java.nio.ByteBuffer;
import java.nio.channels.NonWritableChannelException;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.AccessDeniedException;
import java.nio.file.AccessMode;
import java.nio.file.CopyOption;
import java.nio.file.NoSuchFileException;
import java.nio.file.DirectoryStream;
import java.nio.file.FileStore;
import java.nio.file.FileSystem;
import java.nio.file.FileSystemAlreadyExistsException;
import java.nio.file.FileSystemNotFoundException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.ProviderMismatchException;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributeView;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileAttribute;
import java.nio.file.attribute.FileAttributeView;
import java.nio.file.attribute.FileTime;
import java.nio.file.spi.FileSystemProvider;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Base64;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * File System Provider for Seqera Platform Data-Link Paths
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class SeqeraFileSystemProvider extends FileSystemProvider {

    private static final Logger log = LoggerFactory.getLogger(SeqeraFileSystemProvider.class);
    public static final String SCHEME = "seqera";

    private final Map<String, SeqeraFileSystem> fileSystems = new HashMap<>();

    @Override
    public String getScheme() {
        return SCHEME;
    }

    protected SeqeraPath toSeqeraPath(Path path) {
        if (!(path instanceof SeqeraPath)) {
            throw new ProviderMismatchException();
        }
        return (SeqeraPath) path;
    }

    private void checkScheme(URI uri) {
        final String scheme = uri.getScheme().toLowerCase();
        if (!scheme.equals(getScheme())) {
            throw new IllegalArgumentException("Not a valid " + getScheme().toUpperCase() + " scheme: " + scheme);
        }
    }

    @Override
    public synchronized FileSystem newFileSystem(URI uri, Map<String, ?> config) throws IOException {
        checkScheme(uri);

        // URI format: seqera://datalink-name/optional/path
        final String dataLinkName = uri.getHost();
        if (dataLinkName == null) {
            throw new IllegalArgumentException("Data-link name must be specified in URI host");
        }

        final String key = buildFileSystemKey(uri, config);
        if (fileSystems.containsKey(key)) {
            throw new FileSystemAlreadyExistsException("File system already exists for: " + key);
        }

        // Extract configuration
        final TowerConfig twConfig = new TowerConfig(config, SysEnv.get());

        if (twConfig.getAccessToken() == null) {
            throw new IllegalArgumentException("Access token must be provided via 'token' config or TOWER_ACCESS_TOKEN environment variable");
        }
        // Create HTTP client
        final HxClient httpClient = TowerHxClientFactory.httpClient(twConfig.getAccessToken(), SysEnv.get("TOWER_REFRESH_TOKEN"),twConfig.getEndpoint(), twConfig.getRetryPolicy() );

        final DataLink dataLink = DataLinkUtils.findDataLink(httpClient, twConfig.getEndpoint(), twConfig.getWorkspaceId(), dataLinkName);

        final SeqeraFileSystem fs = new SeqeraFileSystem(this, dataLink, twConfig.getEndpoint(), twConfig.getWorkspaceId(), httpClient);
        fileSystems.put(key, fs);
        return fs;
    }

    private String buildFileSystemKey(URI uri, Map<String, ?> config) {
        final String dataLinkName = uri.getHost();
        final String endpoint = config.containsKey("endpoint")
            ? (String) config.get("endpoint")
            : (System.getenv("TOWER_API_ENDPOINT") != null ? System.getenv("TOWER_API_ENDPOINT") : "https://api.cloud.seqera.io");
        final String workspaceId = config.containsKey("workspaceId")
            ? (String) config.get("workspaceId")
            : (System.getenv("TOWER_WORKSPACE_ID") != null ? System.getenv("TOWER_WORKSPACE_ID") : "");
        return endpoint + ":" + workspaceId + ":" + dataLinkName;
    }

    @Override
    public synchronized FileSystem getFileSystem(URI uri) throws FileSystemNotFoundException {
        checkScheme(uri);
        final String key = buildFileSystemKey(uri, new HashMap<>());
        final SeqeraFileSystem fs = fileSystems.get(key);
        if (fs == null) {
            throw new FileSystemNotFoundException("File system not found for: " + uri);
        }
        return fs;
    }

    public synchronized FileSystem getFileSystemOrCreate(URI uri, Map<String, ?> config) throws IOException {
        checkScheme(uri);
        final String key = buildFileSystemKey(uri, config);
        if (!fileSystems.containsKey(key)) {
            return newFileSystem(uri, config);
        }
        return fileSystems.get(key);
    }

    public synchronized FileSystem getFileSystemOrCreate(URI uri) throws IOException {
        return getFileSystemOrCreate(uri, new HashMap<>());
    }

    @Override
    public SeqeraPath getPath(URI uri) {
        try {
            final SeqeraFileSystem fs = (SeqeraFileSystem) getFileSystemOrCreate(uri);
            return (SeqeraPath) fs.getPath(uri);
        } catch (IOException e) {
            throw new RuntimeException("Failed to get path for URI: " + uri, e);
        }
    }

    @Override
    public OutputStream newOutputStream(Path path, OpenOption... options) throws IOException {
        final SeqeraPath seqPath = toSeqeraPath(path);
        final SeqeraFileSystem fs = (SeqeraFileSystem) seqPath.getFileSystem();

        // Check if write is requested
        boolean write = false;
        for (OpenOption opt : options) {
            if (opt == StandardOpenOption.WRITE || opt == StandardOpenOption.CREATE ||
                opt == StandardOpenOption.CREATE_NEW || opt == StandardOpenOption.TRUNCATE_EXISTING) {
                write = true;
                break;
            }
        }

        if (!write) {
            throw new IllegalArgumentException("Write option must be specified");
        }

        // Create a temporary file and return an OutputStream that will upload on close
        final Path tempFile = Files.createTempFile("seqera-upload-", ".tmp");
        return new SeqeraOutputStream(tempFile, seqPath, fs);
    }

    private static class SeqeraOutputStream extends OutputStream {
        private final OutputStream delegate;
        private final Path tempFile;
        private final SeqeraPath targetPath;
        private final SeqeraFileSystem fs;

        SeqeraOutputStream(Path tempFile, SeqeraPath targetPath, SeqeraFileSystem fs) throws IOException {
            this.tempFile = tempFile;
            this.targetPath = targetPath;
            this.fs = fs;
            this.delegate = Files.newOutputStream(tempFile);
        }

        @Override
        public void write(int b) throws IOException {
            delegate.write(b);
        }

        @Override
        public void write(byte[] b) throws IOException {
            delegate.write(b);
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException {
            delegate.write(b, off, len);
        }

        @Override
        public void flush() throws IOException {
            delegate.flush();
        }

        @Override
        public void close() throws IOException {
            delegate.close();
            try {
                // Upload the file to Seqera Platform
                DataLinkUtils.uploadFile(
                    fs.getHttpClient(),
                    fs.getEndpoint(),
                    fs.getDataLink(),
                    targetPath.getPathForApi(),
                    tempFile,
                    fs.getWorkspaceId()
                );
            } finally {
                // Clean up temp file
                Files.deleteIfExists(tempFile);
            }
        }
    }

    @Override
    public InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        final SeqeraPath seqPath = toSeqeraPath(path);
        final SeqeraFileSystem fs = (SeqeraFileSystem) seqPath.getFileSystem();

        // Download to temp file and return input stream
        final Path tempFile = Files.createTempFile("seqera-download-", ".tmp");
        try {
            DataLinkUtils.downloadFile(
                fs.getHttpClient(),
                fs.getEndpoint(),
                fs.getDataLink(),
                seqPath.getPathForApi(),
                tempFile,
                fs.getWorkspaceId()
            );
            return new DeleteOnCloseInputStream(tempFile);
        } catch (Exception e) {
            Files.deleteIfExists(tempFile);
            throw e;
        }
    }

    private static class DeleteOnCloseInputStream extends InputStream {
        private final InputStream delegate;
        private final Path tempFile;

        DeleteOnCloseInputStream(Path tempFile) throws IOException {
            this.tempFile = tempFile;
            this.delegate = Files.newInputStream(tempFile);
        }

        @Override
        public int read() throws IOException {
            return delegate.read();
        }

        @Override
        public int read(byte[] b) throws IOException {
            return delegate.read(b);
        }

        @Override
        public int read(byte[] b, int off, int len) throws IOException {
            return delegate.read(b, off, len);
        }

        @Override
        public void close() throws IOException {
            try {
                delegate.close();
            } finally {
                Files.deleteIfExists(tempFile);
            }
        }
    }

    @Override
    public SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        final SeqeraPath seqPath = toSeqeraPath(path);

        // Check if this is a read or write operation
        boolean write = options.contains(StandardOpenOption.WRITE) ||
                       options.contains(StandardOpenOption.CREATE) ||
                       options.contains(StandardOpenOption.CREATE_NEW);

        if (write) {
            throw new UnsupportedOperationException("SeekableByteChannel write not supported - use newOutputStream instead");
        }

        // For read, create a temp file and return a channel
        final SeqeraFileSystem fs = (SeqeraFileSystem) seqPath.getFileSystem();
        final Path tempFile = Files.createTempFile("seqera-channel-", ".tmp");

        try {
            DataLinkUtils.downloadFile(
                fs.getHttpClient(),
                fs.getEndpoint(),
                fs.getDataLink(),
                seqPath.getPathForApi(),
                tempFile,
                fs.getWorkspaceId()
            );
            return new DeleteOnCloseSeekableByteChannel(tempFile);
        } catch (Exception e) {
            Files.deleteIfExists(tempFile);
            throw e;
        }
    }

    private static class DeleteOnCloseSeekableByteChannel implements SeekableByteChannel {
        private final SeekableByteChannel delegate;
        private final Path tempFile;

        DeleteOnCloseSeekableByteChannel(Path tempFile) throws IOException {
            this.tempFile = tempFile;
            this.delegate = Files.newByteChannel(tempFile);
        }

        @Override
        public int read(ByteBuffer dst) throws IOException {
            return delegate.read(dst);
        }

        @Override
        public int write(ByteBuffer src) throws IOException {
            throw new NonWritableChannelException();
        }

        @Override
        public long position() throws IOException {
            return delegate.position();
        }

        @Override
        public SeekableByteChannel position(long newPosition) throws IOException {
            delegate.position(newPosition);
            return this;
        }

        @Override
        public long size() throws IOException {
            return delegate.size();
        }

        @Override
        public SeekableByteChannel truncate(long size) throws IOException {
            throw new NonWritableChannelException();
        }

        @Override
        public boolean isOpen() {
            return delegate.isOpen();
        }

        @Override
        public void close() throws IOException {
            try {
                delegate.close();
            } finally {
                Files.deleteIfExists(tempFile);
            }
        }
    }

    @Override
    public DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final SeqeraPath seqPath = toSeqeraPath(dir);
        final SeqeraFileSystem fs = (SeqeraFileSystem) seqPath.getFileSystem();

        // Use browse API to list files
        final List<Path> files = listFiles(fs, seqPath.getPathForApi());

        return new DirectoryStream<Path>() {
            @Override
            public Iterator<Path> iterator() {
                return files.iterator();
            }

            @Override
            public void close() throws IOException {
                // Nothing to close
            }
        };
    }

    private List<Path> listFiles(SeqeraFileSystem fs, String path) throws IOException {
        // Use utility method to get file names
        final List<String> fileNames = DataLinkUtils.listFiles(
            fs.getHttpClient(),
            fs.getEndpoint(),
            fs.getDataLink(),
            path,
            fs.getWorkspaceId()
        );

        // Convert file names to Path objects
        final List<Path> result = new ArrayList<>();
        for (String name : fileNames) {
            final String fullPath = path != null && !path.isEmpty() ? path + "/" + name : name;
            result.add(new SeqeraPath(fs, fullPath));
        }
        return result;
    }

    @Override
    public void createDirectory(Path dir, FileAttribute<?>... attrs) throws IOException {
        // Data-links don't require explicit directory creation
        log.debug("Directory creation not required for Seqera paths: {}", dir);
    }

    @Override
    public void delete(Path path) throws IOException {
        final SeqeraPath seqPath = toSeqeraPath(path);
        final SeqeraFileSystem fs = (SeqeraFileSystem) seqPath.getFileSystem();

        // Use utility method to delete file
        DataLinkUtils.deleteFile(
            fs.getHttpClient(),
            fs.getEndpoint(),
            fs.getDataLink(),
            seqPath.getPathForApi(),
            fs.getWorkspaceId()
        );
    }

    @Override
    public void copy(Path source, Path target, CopyOption... options) throws IOException {
        // Determine if source and target are Seqera paths or local paths
        final boolean sourceIsSeqera = source instanceof SeqeraPath;
        final boolean targetIsSeqera = target instanceof SeqeraPath;

        if (!sourceIsSeqera && !targetIsSeqera) {
            // Both are local paths - not our responsibility
            throw new UnsupportedOperationException("Copy between local paths should use Files.copy()");
        }

        if (sourceIsSeqera && targetIsSeqera) {
            // Both are Seqera paths - server-side copy not supported yet
            throw new UnsupportedOperationException("Copy between Seqera paths not supported - download and re-upload instead");
        }

        // Handle copy from local to Seqera or Seqera to local
        if (sourceIsSeqera) {
            // Download from Seqera to local
            copyFromSeqeraToLocal((SeqeraPath) source, target);
        } else {
            // Upload from local to Seqera
            copyFromLocalToSeqera(source, (SeqeraPath) target);
        }
    }

    private void copyFromSeqeraToLocal(SeqeraPath source, Path target) throws IOException {
        final SeqeraFileSystem fs = (SeqeraFileSystem) source.getFileSystem();

        // Check if source is a directory or file
        final BasicFileAttributes attrs = readAttributes(source, BasicFileAttributes.class);

        if (attrs.isDirectory()) {
            // Download folder recursively
            DataLinkUtils.downloadFolder(
                fs.getHttpClient(),
                fs.getEndpoint(),
                fs.getDataLink(),
                source.getPathForApi(),
                target,
                fs.getWorkspaceId()
            );
        } else {
            // Download single file
            DataLinkUtils.downloadFile(
                fs.getHttpClient(),
                fs.getEndpoint(),
                fs.getDataLink(),
                source.getPathForApi(),
                target,
                fs.getWorkspaceId()
            );
        }
    }

    private void copyFromLocalToSeqera(Path source, SeqeraPath target) throws IOException {
        final SeqeraFileSystem fs = (SeqeraFileSystem) target.getFileSystem();

        if (Files.isDirectory(source)) {
            // Upload folder recursively
            DataLinkUtils.uploadFolder(
                fs.getHttpClient(),
                fs.getEndpoint(),
                fs.getDataLink(),
                target.getPathForApi(),
                source,
                fs.getWorkspaceId()
            );
        } else {
            // Upload single file
            DataLinkUtils.uploadFile(
                fs.getHttpClient(),
                fs.getEndpoint(),
                fs.getDataLink(),
                target.getPathForApi(),
                source,
                fs.getWorkspaceId()
            );
        }
    }

    @Override
    public void move(Path source, Path target, CopyOption... options) throws IOException {
        throw new UnsupportedOperationException("Move not supported");
    }

    @Override
    public boolean isSameFile(Path path, Path path2) throws IOException {
        return path.equals(path2);
    }

    @Override
    public boolean isHidden(Path path) throws IOException {
        return false;
    }

    @Override
    public FileStore getFileStore(Path path) throws IOException {
        throw new UnsupportedOperationException("File store not supported");
    }

    @Override
    public void checkAccess(Path path, AccessMode... modes) throws IOException {
        final SeqeraPath seqPath = toSeqeraPath(path);
        final SeqeraFileSystem fs = (SeqeraFileSystem) seqPath.getFileSystem();

        for (AccessMode m : modes) {
            if (m == AccessMode.EXECUTE) {
                throw new AccessDeniedException("Execute mode not supported");
            }
        }

        // Check if the file exists using getFileDetails
        final DataLinkUtils.DataLinkItem item = DataLinkUtils.getFileDetails(
            fs.getHttpClient(),
            fs.getEndpoint(),
            fs.getDataLink(),
            seqPath.getPathForApi(),
            fs.getWorkspaceId()
        );

        if (item == null) {
            throw new NoSuchFileException(path.toString());
        }
    }

    @Override
    @SuppressWarnings("unchecked")
    public <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        if (type == BasicFileAttributeView.class) {
            return (V) new SeqeraFileAttributeView(toSeqeraPath(path));
        }
        return null;
    }

    @Override
    @SuppressWarnings("unchecked")
    public <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        if (type == BasicFileAttributes.class || BasicFileAttributes.class.isAssignableFrom(type)) {
            return (A) new SeqeraFileAttributes(toSeqeraPath(path));
        }
        throw new UnsupportedOperationException("Attributes type not supported: " + type);
    }

    @Override
    public Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Read file attributes not supported");
    }

    @Override
    public void setAttribute(Path path, String attribute, Object value, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("Set file attributes not supported");
    }

    /**
     * Check if a file attribute view is supported
     */
    public boolean supportsFileAttributeView(String name) {
        return "basic".equals(name);
    }

    private static class SeqeraFileAttributeView implements BasicFileAttributeView {
        private final SeqeraPath path;

        SeqeraFileAttributeView(SeqeraPath path) {
            this.path = path;
        }

        @Override
        public String name() {
            return "basic";
        }

        @Override
        public BasicFileAttributes readAttributes() throws IOException {
            return new SeqeraFileAttributes(path);
        }

        @Override
        public void setTimes(FileTime lastModifiedTime, FileTime lastAccessTime, FileTime createTime) throws IOException {
            throw new UnsupportedOperationException("Set times not supported");
        }
    }

    private static class SeqeraFileAttributes implements BasicFileAttributes {
        private final SeqeraPath path;
        private final FileTime creationTime;
        private final FileTime lastModifiedTime;
        private final boolean isDirectory;
        private final long size;

        SeqeraFileAttributes(SeqeraPath path) throws IOException {
            this.path = path;
            this.creationTime = FileTime.from(Instant.now());
            this.lastModifiedTime = FileTime.from(Instant.now());

            final String pathStr = path.getPathForApi();

            // If path is empty, it's the root which is a directory
            if (pathStr == null || pathStr.isEmpty()) {
                this.isDirectory = true;
                this.size = 0L;
            } else {
                // Fetch file details from the API
                final SeqeraFileSystem fs = (SeqeraFileSystem) path.getFileSystem();
                final DataLinkUtils.DataLinkItem item = DataLinkUtils.getFileDetails(
                    fs.getHttpClient(),
                    fs.getEndpoint(),
                    fs.getDataLink(),
                    pathStr,
                    fs.getWorkspaceId()
                );

                if (item == null) {
                    throw new NoSuchFileException(path.toString());
                }

                this.isDirectory = item.getType() == DataLinkUtils.DataLinkItemType.FOLDER;
                this.size = item.getSize() != null ? item.getSize().longValue() : 0L;
            }
        }

        @Override
        public FileTime lastModifiedTime() {
            return lastModifiedTime;
        }

        @Override
        public FileTime lastAccessTime() {
            return lastModifiedTime;
        }

        @Override
        public FileTime creationTime() {
            return creationTime;
        }

        @Override
        public boolean isRegularFile() {
            return !isDirectory;
        }

        @Override
        public boolean isDirectory() {
            return isDirectory;
        }

        @Override
        public boolean isSymbolicLink() {
            return false;
        }

        @Override
        public boolean isOther() {
            return false;
        }

        @Override
        public long size() {
            return size;
        }

        @Override
        public Object fileKey() {
            return null;
        }
    }
}
