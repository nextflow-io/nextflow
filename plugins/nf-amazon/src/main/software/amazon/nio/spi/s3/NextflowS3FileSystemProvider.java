package software.amazon.nio.spi.s3;

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
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.FileSystemAlreadyExistsException;
import java.nio.file.FileSystemNotFoundException;
import java.nio.file.Path;
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
			//
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
        return new S3FileSystem(this, config);
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
		if (!uri.getScheme().equals(getScheme())) {
            throw new IllegalArgumentException("URI scheme must be " + getScheme());
        }
		return getFileSystem(uri).getPath(uri.getPath());
	}


}
