/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.client

import com.google.api.gax.rpc.StatusCode
import com.google.auth.oauth2.GoogleCredentials
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.google.GoogleOpts
import nextflow.util.MemoryUnit
import nextflow.util.Duration
import java.util.stream.Collectors
/**
 * Model Google Batch config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchConfig {

    static final private int DEFAULT_MAX_SPOT_ATTEMPTS = 0
    
    static final private List<Integer> DEFAULT_RETRY_LIST = List.of(50001)

    // These are the default values provided by BatchServiceStubSettings.java
    // https://github.com/googleapis/google-cloud-java/blob/01a581887a8d5e8300208d64472bccddccf333ca/java-batch/google-cloud-batch/src/main/java/com/google/cloud/batch/v1/stub/BatchServiceStubSettings.java#L507C1-L507C1
    static private final Duration DEFAULT_GRPC_INITIAL_RETRY_DELAY = Duration.of("1s")
    static private final Duration DEFAULT_GRPC_MAX_RETRY_DELAY = Duration.of("10s")
    static private final Duration DEFAULT_GRPC_INITIAL_RPC_TIMEOUT = Duration.of("1m")
    static private final Duration DEFAULT_GRPC_MAX_RPC_TIMEOUT = Duration.of("1m")
    static private final Duration DEFAULT_GRPC_TOTAL_TIMEOUT = Duration.of("1m")
    static private final double DEFAULT_GRPC_TIMEOUT_MULTIPLIER = 1.0d
    static private final double DEFAULT_GRPC_RETRY_DELAY_MULTIPLIER = 1.3d
    // This is the default set of retryable status codes defined
    // https://github.com/googleapis/google-cloud-java/blob/01a581887a8d5e8300208d64472bccddccf333ca/java-batch/google-cloud-batch/src/main/java/com/google/cloud/batch/v1/stub/BatchServiceStubSettings.java#L487
    static private final Set<StatusCode.Code> DEFAULT_GRPC_RETRYABLE_STATUS_CODES = Set.of(StatusCode.Code.UNAVAILABLE)

    private GoogleOpts googleOpts
    private GoogleCredentials credentials
    private List<String> allowedLocations
    private String bootDiskImage
    private MemoryUnit bootDiskSize
    private String cpuPlatform
    private int maxSpotAttempts
    private boolean installGpuDrivers
    private boolean preemptible
    private boolean spot
    private boolean usePrivateAddress
    private String network
    private String subnetwork
    private String serviceAccountEmail
    private BatchRetryConfig retryConfig
    private List<Integer> autoRetryExitCodes
    private Set<StatusCode.Code> grpcRetryableCodes = DEFAULT_GRPC_RETRYABLE_STATUS_CODES
    private Duration grpcInitialRetryDelay = DEFAULT_GRPC_INITIAL_RETRY_DELAY
    private Duration grpcMaxRetryDelay = DEFAULT_GRPC_MAX_RETRY_DELAY
    private Duration grpcInitialRpcTimeout = DEFAULT_GRPC_INITIAL_RPC_TIMEOUT
    private Duration grpcMaxRpcTimeout = DEFAULT_GRPC_MAX_RPC_TIMEOUT
    private Duration grpcTotalTimeout = DEFAULT_GRPC_TOTAL_TIMEOUT
    private double grpcTimeoutMultiplier = DEFAULT_GRPC_TIMEOUT_MULTIPLIER
    private double grpcRetryDelayMultiplier = DEFAULT_GRPC_RETRY_DELAY_MULTIPLIER

    GoogleOpts getGoogleOpts() { return googleOpts }
    GoogleCredentials getCredentials() { return credentials }
    List<String> getAllowedLocations() { allowedLocations }
    String getBootDiskImage() { bootDiskImage }
    MemoryUnit getBootDiskSize() { bootDiskSize }
    String getCpuPlatform() { cpuPlatform }
    int getMaxSpotAttempts() { maxSpotAttempts }
    boolean getInstallGpuDrivers() { installGpuDrivers }
    boolean getPreemptible() { preemptible }
    boolean getSpot() { spot }
    boolean getUsePrivateAddress() { usePrivateAddress }
    String getNetwork() { network }
    String getSubnetwork() { subnetwork }
    String getServiceAccountEmail() { serviceAccountEmail }
    BatchRetryConfig getRetryConfig() { retryConfig }
    List<Integer> getAutoRetryExitCodes() { autoRetryExitCodes }
    Set<StatusCode.Code> getGrpcRetryableCodes() { grpcRetryableCodes }
    Duration getGrpcInitialRetryDelay() { grpcInitialRetryDelay }
    Duration getGrpcMaxRetryDelay() { grpcMaxRetryDelay }
    Duration getGrpcInitialRpcTimeout() { grpcInitialRpcTimeout }
    Duration getGrpcMaxRpcTimeout() { grpcMaxRpcTimeout }
    Duration getGrpcTotalTimeout() { grpcTotalTimeout }
    double getGrpcTimeoutMultiplier() { grpcTimeoutMultiplier }
    double getGrpcRetryDelayMultiplier() { grpcRetryDelayMultiplier }

    static BatchConfig create(Session session) {
        final result = new BatchConfig()
        result.googleOpts = GoogleOpts.create(session)
        result.credentials = result.googleOpts.credentials
        result.allowedLocations = session.config.navigate('google.batch.allowedLocations', List.of()) as List<String>
        result.bootDiskImage = session.config.navigate('google.batch.bootDiskImage')
        result.bootDiskSize = session.config.navigate('google.batch.bootDiskSize') as MemoryUnit
        result.cpuPlatform = session.config.navigate('google.batch.cpuPlatform')
        result.maxSpotAttempts = session.config.navigate('google.batch.maxSpotAttempts', DEFAULT_MAX_SPOT_ATTEMPTS) as int
        result.installGpuDrivers = session.config.navigate('google.batch.installGpuDrivers',false)
        result.preemptible = session.config.navigate('google.batch.preemptible',false)
        result.spot = session.config.navigate('google.batch.spot',false)
        result.usePrivateAddress = session.config.navigate('google.batch.usePrivateAddress',false)
        result.network = session.config.navigate('google.batch.network')
        result.subnetwork = session.config.navigate('google.batch.subnetwork')
        result.serviceAccountEmail = session.config.navigate('google.batch.serviceAccountEmail')
        result.retryConfig = new BatchRetryConfig( session.config.navigate('google.batch.retryPolicy') as Map ?: Map.of() )
        result.autoRetryExitCodes = session.config.navigate('google.batch.autoRetryExitCodes', DEFAULT_RETRY_LIST) as List<Integer>
        result.grpcRetryableCodes = convertToEnum(
            session.config.navigate("google.grpc.retryableCodes", List.of("UNAVAILABLE")) as List<String>,
            StatusCode.Code.class) as Set<StatusCode.Code>
        result.grpcInitialRetryDelay = session.config.navigate("google.grpc.initialRetryDelay", DEFAULT_GRPC_INITIAL_RETRY_DELAY) as Duration
        result.grpcMaxRetryDelay = session.config.navigate("google.grpc.maxRetryDelay", DEFAULT_GRPC_MAX_RETRY_DELAY) as Duration
        result.grpcInitialRpcTimeout = session.config.navigate("google.grpc.initialRpcTimeout", DEFAULT_GRPC_INITIAL_RPC_TIMEOUT) as Duration
        result.grpcMaxRpcTimeout = session.config.navigate("google.grpc.maxRpcTimeout", DEFAULT_GRPC_MAX_RPC_TIMEOUT) as Duration
        result.grpcTotalTimeout = session.config.navigate("google.grpc.totalTimeout", DEFAULT_GRPC_TOTAL_TIMEOUT) as Duration
        result.grpcTimeoutMultiplier = session.config.navigate("google.grpc.timeoutMultiplier", DEFAULT_GRPC_TIMEOUT_MULTIPLIER) as double
        result.grpcRetryDelayMultiplier = session.config.navigate("google.grpc.retryDelayMultiplier", DEFAULT_GRPC_RETRY_DELAY_MULTIPLIER) as double
        return result
    }

    static private Set<? extends Enum> convertToEnum(final Collection<String> enumNames, final Class<? extends Enum> type) {
        return enumNames.stream()
                .map( (name) -> { return Enum.valueOf(type, name)})
                .collect(Collectors.toUnmodifiableSet())
    }

    @Override
    String toString(){
        return "BatchConfig[googleOpts=$googleOpts"
    }

}
