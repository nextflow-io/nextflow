/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.container;

import java.nio.file.Path;
import java.util.UUID;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;

import nextflow.dataflow.DataflowThreads;
import nextflow.util.Duration;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Reads of the singularity and charliecloud caches on dataflow threads, including the loop thread
 * of {@code NXF_DATAFLOW_POOL_SIZE=1}
 */
class ContainerCacheDataflowTest {

    interface Download { Path apply(String url); }

    static Function<String, Path> cache(String type, Download download) {
        if( type.equals("singularity") ) {
            var cache = new SingularityCache("/cache", null, false, Duration.of("1 min")) {
                @Override
                Path downloadContainerImage(String imageUrl) {
                    return download.apply(imageUrl);
                }
            };
            return cache::getCachePathFor;
        }
        var cache = new CharliecloudCache() {
            @Override
            Path downloadCharliecloudImage(String imageUrl) {
                return download.apply(imageUrl);
            }
        };
        return cache::getCachePathFor;
    }

    @ParameterizedTest
    @ValueSource(strings = {"singularity", "charliecloud"})
    void shouldPullTheImageOnceForConcurrentReaders(String type) throws Exception {
        var count = new AtomicInteger();
        var cache = cache(type, url -> {
            count.incrementAndGet();
            sleep();
            return Path.of("/cache", url);
        });
        var url = "docker://foo/" + UUID.randomUUID();

        var paths = DataflowThreads.call(8, () -> cache.apply(url));

        assertThat(paths).hasSize(8).containsOnly(Path.of("/cache", url));
        assertThat(count).hasValue(1);
    }

    @ParameterizedTest
    @ValueSource(strings = {"singularity", "charliecloud"})
    void shouldPropagateThePullError(String type) {
        var count = new AtomicInteger();
        var cache = cache(type, url -> {
            count.incrementAndGet();
            throw new IllegalArgumentException(type + " boom");
        });
        var url = "docker://foo/" + UUID.randomUUID();

        for( int i = 0; i < 2; i++ ) {
            assertThatThrownBy(() -> DataflowThreads.call(1, () -> cache.apply(url)))
                .isInstanceOf(IllegalStateException.class)
                .hasCauseInstanceOf(IllegalArgumentException.class)
                .hasRootCauseMessage(type + " boom");
        }
        assertThat(count).hasValue(1);
    }

    static void sleep() {
        try {
            Thread.sleep(100);
        }
        catch( InterruptedException e ) {
            throw new RuntimeException(e);
        }
    }
}
