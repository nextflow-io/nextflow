package software.amazon.nio.spi.s3;

import static software.amazon.nio.spi.s3.NextflowS3Path.*;
import nextflow.SysEnv;
import nextflow.cloud.aws.config.AwsConfig;
import nextflow.cloud.aws.config.AwsS3Config;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider;
import software.amazon.awssdk.core.exception.SdkClientException;
import software.amazon.awssdk.regions.Region;
import software.amazon.awssdk.regions.providers.InstanceProfileRegionProvider;
import software.amazon.awssdk.services.s3.S3AsyncClient;
import software.amazon.awssdk.services.s3.S3CrtAsyncClientBuilder;
import software.amazon.nio.spi.s3.config.S3NioSpiConfiguration;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileAttribute;
import java.nio.file.attribute.FileAttributeView;
import java.util.*;

public class NextflowS3FileSystemProvider extends S3FileSystemProvider {

    private static final Logger log = LoggerFactory.getLogger(NextflowS3FileSystemProvider.class);
    final Map<String, S3FileSystem> fileSystems = new HashMap<>();

    @Override
	public FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        if (!uri.getScheme().equals(getScheme())) {
            throw new IllegalArgumentException("URI scheme must be " + getScheme());
        }

        String bucketName = fileSystemInfo(uri).bucket();
        log.debug("Creating filesystem for S3 bucket {}", bucketName);
        synchronized (fileSystems) {
			if( fileSystems.containsKey(bucketName))
				throw new FileSystemAlreadyExistsException("S3 filesystem already exists. Use getFileSystem() instead");

			final AwsConfig awsConfig = new AwsConfig(env);
			final S3FileSystem result = createFileSystem(uri, awsConfig);
			fileSystems.put(bucketName, result);
			return result;
		}

	}

    private S3FileSystem createFileSystem(URI uri, AwsConfig awsConfig) {
        var info = fileSystemInfo(uri);
        var config = new S3NioSpiConfiguration().withEndpoint(info.endpoint()).withBucketName(info.bucket());
        if (info.accessKey() != null) {
            config.withCredentials(info.accessKey(), info.accessSecret());
        }
        S3ClientProvider clientProvider = createS3clientProvider( awsConfig, config);
        S3FileSystem fs = new S3FileSystem(this, config);
        fs.clientProvider(clientProvider);
        return fs;
    }

    private S3ClientProvider createS3clientProvider(AwsConfig awsConfig, S3NioSpiConfiguration spiConfig) {
        // try to load amazon props
        Properties props = loadAmazonProperties();
        // add properties for legacy compatibility
        props.putAll(awsConfig.getS3LegacyProperties());

        S3AsyncClientConfiguration clientConfig = S3AsyncClientConfiguration.create(props);
        var region = getRegion(awsConfig);
        var s3Config = awsConfig.getS3Config();
        if (s3Config != null)
            addS3ConfigurationToSpiConfig(s3Config, region, spiConfig);
        spiConfig.withOpenOptions(clientConfig.getOpenOptions());
        var clientProvider = new S3ClientProvider(spiConfig);
        clientProvider.asyncClientBuilder(getAsyncClientBuilder(clientConfig));
        return clientProvider;
    }

    private S3CrtAsyncClientBuilder getAsyncClientBuilder(S3AsyncClientConfiguration clientConfig) {
        S3CrtAsyncClientBuilder builder = S3AsyncClient.crtBuilder().crossRegionAccessEnabled(true);
        if( clientConfig.hasHttpConfiguration() )
            builder.httpConfiguration(clientConfig.getHttpConfiguration());
        if( clientConfig.hasRetryConfiguration() )
            builder.retryConfiguration(clientConfig.getRetryConfiguration());
        if( clientConfig.hasMaxConcurrency() )
            builder.maxConcurrency(clientConfig.getMaxConcurrency());
        return builder;
    }

    private String getRegion(AwsConfig awsConfig){
        if( awsConfig.getRegion() != null && !awsConfig.getRegion().isBlank() )
            return awsConfig.getRegion();
        if( SysEnv.get("AWS_REGION") != null )
            return SysEnv.get("AWS_REGION");
        if( SysEnv.get("AWS_DEFAULT_REGION") != null )
            return SysEnv.get("AWS_DEFAULT_REGION");
        try {
            return new InstanceProfileRegionProvider().getRegion().id();
        } catch ( SdkClientException e) {
            log.warn("Cannot fetch AWS region either from configuration environment and instance. Setting default US_EAST_1");
            return Region.US_EAST_1.id();
        }
    }

    private void addS3ConfigurationToSpiConfig(AwsS3Config s3config, String region, S3NioSpiConfiguration spiConfig){
        spiConfig.withForcePathStyle(s3config.getPathStyleAccess());
        if (s3config.getAnonymous() != null)
            spiConfig.withCredentials(AnonymousCredentialsProvider.create().resolveCredentials());

        if (s3config.getEndpoint() != null && !s3config.getEndpoint().isBlank()) {
            spiConfig.withEndpoint(s3config.getEndpoint());
        } else {
            spiConfig.withRegion(region);
        }

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

    @Override
	public FileSystem getFileSystem(URI uri) {
        final String bucketName = fileSystemInfo(uri).bucket();
        log.debug("Getting filesystem for S3 bucket {}", bucketName);
        synchronized (fileSystems) {
            final FileSystem fileSystem = this.fileSystems.get(bucketName);

            if (fileSystem == null) {
                throw new FileSystemNotFoundException("S3 filesystem not yet created. Use newFileSystem() instead");
            }
            return fileSystem;
        }
	}

	/**
     * Deviation from spec: throws FileSystemNotFoundException if FileSystem
     * hasn't yet been initialized. Call newFileSystem() first.
     * Need credentials. Maybe set credentials after? how?
     */
	@Override
	public Path getPath(URI uri) {
		if (!uri.getScheme().equals(getScheme())) {
            throw new IllegalArgumentException("URI scheme must be " + getScheme());
        }
		return new NextflowS3Path((S3Path) getFileSystem(uri).getPath(uri.getPath()));
	}

    @Override
    public InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        return super.newInputStream(unwrapS3Path(path), options);
    }

    @Override
    public SeekableByteChannel newByteChannel(
        Path path,
        Set<? extends OpenOption> options,
        FileAttribute<?>... attrs
    ) throws IOException {
        return super.newByteChannel(unwrapS3Path(path), options, attrs);
    }

    @Override
    public OutputStream newOutputStream(Path path, OpenOption... options) throws IOException {
        if( path instanceof NextflowS3Path ){
            final NextflowS3Path nxfS3Path = (NextflowS3Path) path;
            if( options != null && options.length > 0)
                return super.newOutputStream(nxfS3Path.toS3Path(), updateOptions(options, nxfS3Path.getOpenOptions()));
            else
                return super.newOutputStream(nxfS3Path.toS3Path(), nxfS3Path.getOpenOptions());
        }
        return super.newOutputStream(path, options);
    }

    private OpenOption[] updateOptions(OpenOption[] options, NextflowS3PathOpenOptions newOption) {
        final OpenOption[] newOptions = Arrays.copyOf(options, options.length + 1 );
        newOptions[options.length] = newOption;
        return newOptions;
    }

    @Override
    public void createSymbolicLink(Path link, Path target, FileAttribute<?>... attrs) throws IOException {
        super.createSymbolicLink(unwrapS3Path(link), unwrapS3Path(target), attrs);
    }

    @Override
    public void createLink(Path link, Path existing) throws IOException {
        super.createLink(unwrapS3Path(link), unwrapS3Path(existing));
    }

    @Override
    public boolean deleteIfExists(Path path) throws IOException {
        return super.deleteIfExists(unwrapS3Path(path));
    }

    @Override
    public Path readSymbolicLink(Path link) throws IOException {
        return super.readSymbolicLink(unwrapS3Path(link));
    }

    @Override
    public <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        return super.readAttributes(unwrapS3Path(path), type, options);
    }

    @Override
    public void copy(Path source, Path target, CopyOption... options) throws IOException {
        super.copy(unwrapS3Path(source), unwrapS3Path(target), options);
    }

    @Override
    public void move(Path source, Path target, CopyOption... options) throws IOException {
        super.move(unwrapS3Path(source), unwrapS3Path(target), options);
    }

    @Override
    public <V extends FileAttributeView> V getFileAttributeView(Path path, Class<V> type, LinkOption... options) {
        return super.getFileAttributeView(unwrapS3Path(path), type, options);
    }

    @Override
    public void createDirectory(Path dir, FileAttribute<?>... attrs) throws IOException {
        super.createDirectory(unwrapS3Path(dir), attrs);
    }

    @Override
    public void delete(Path path) throws IOException {
        super.delete(unwrapS3Path(path));
    }

    @Override
    public DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        return super.newDirectoryStream(unwrapS3Path(dir), filter);
    }
}
