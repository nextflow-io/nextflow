import com.google.cloud.storage.Bucket;
import com.google.cloud.storage.BucketInfo;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import com.upplication.s3fs.util.IOUtils;
import com.upplication.s3fs.util.S3MultipartOptions;
import com.upplication.s3fs.util.S3ObjectSummaryLookup;
import com.upplication.s3fs.util.S3UploadRequest;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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

import static com.google.common.collect.Sets.difference;
import static java.lang.String.format;

public class GSFileSystemProvider extends FileSystemProvider {

    final AtomicReference<GSFileSystem> fileSystem = new AtomicReference<>();

    @Override
    public String getScheme() {
        return "gs";
    }

    @Override
    public FileSystem newFileSystem(URI uri, Map<String, ?> env)
            throws IOException {
        Preconditions.checkNotNull(uri, "uri is null");
        Preconditions.checkArgument(uri.getScheme().equals("gs"),
                "uri scheme must be 'gs': '%s'", uri);

    }

    @Override
    public FileSystem getFileSystem(URI uri) {
        FileSystem fileSystem = this.fileSystem.get();

        if (fileSystem == null) {
            throw new FileSystemNotFoundException(
                    String.format("GC fileSystem not yet created. Use newFileSystem() instead"));
        }

        return fileSystem;
    }

    @Override
    public Path getPath(URI uri) {
        Preconditions.checkArgument(uri.getScheme().equals(getScheme()),
                "URI scheme must be %s", getScheme());

        if (uri.getHost() != null && !uri.getHost().isEmpty() &&
                !uri.getHost().equals(fileSystem.get().getEndpoint())) {
            throw new IllegalArgumentException(format(
                    "only empty URI host that matching the current fileSystem: %s",
                    fileSystem.get().getEndpoint()));
        }
        return getFileSystem(uri).getPath(uri.getPath());
    }

    @Override
    public DirectoryStream<Path> newDirectoryStream(Path dir,
                                                    DirectoryStream.Filter<? super Path> filter) throws IOException {

        Preconditions.checkArgument(dir instanceof GSPath,
                "path must be an instance of %s", GSPath.class.getName());
        final GSPath gsPath = (GSPath) dir;

        return new DirectoryStream<Path>() {
            @Override
            public void close() throws IOException {
                // nothing to do here
            }

            @Override
            public Iterator<Path> iterator() {
                return new GSIterator
            }
        }
    }
}
