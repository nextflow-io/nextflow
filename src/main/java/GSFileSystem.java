import com.google.api.gax.paging.Page;
import com.google.cloud.storage.Bucket;
import com.google.cloud.storage.Storage;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;

import java.io.IOException;
import java.nio.file.*;
import java.nio.file.attribute.UserPrincipalLookupService;
import java.nio.file.spi.FileSystemProvider;
import java.util.Set;

public class GSFileSystem extends FileSystem {

    private final GSFileSystemProvider provider;
    private final String endpoint;
    private final Storage storage;


    public GSFileSystem(GSFileSystemProvider provider, Storage storage, String endpoint) {
        this.provider = provider;
        this.storage = storage;
        this.endpoint = endpoint;
    }

    @Override
    public FileSystemProvider provider() {
        return provider;
    }

    @Override
    public void close() throws IOException {
        this.provider.fileSystem.compareAndSet(this, null);
    }

    @Override
    public boolean isOpen() {
        return this.provider.fileSystem.get() != null;
    }

    @Override
    public boolean isReadOnly() {
        return false;
    }

    @Override
    public String getSeparator() {
        return GSPath.PATH_SEPARATOR;
    }

    @Override
    public Iterable<Path> getRootDirectories() {
        ImmutableList.Builder<Path> builder = ImmutableList.builder();
        Page<Bucket> buckets = storage.list();

        for (Bucket bucket : buckets.iterateAll()) {
            builder.add(new GSPath(this, bucket.getName()));
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
            return new GSPath(this, first);
        }

        return new GSPath(this, first, more);
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

    public Storage getStorage() {
        return storage;
    }

    public String getEndpoint() {
        return endpoint;
    }
}