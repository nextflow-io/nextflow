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

package nextflow.spack;

import java.nio.file.Path;
import java.util.Map;
import java.util.UUID;
import java.util.concurrent.atomic.AtomicInteger;

import nextflow.dataflow.DataflowThreads;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Reads of the spack cache on dataflow threads, including the loop thread of {@code NXF_DATAFLOW_POOL_SIZE=1}
 */
class SpackCacheDataflowTest {

    @TempDir
    Path dir;

    @Test
    void shouldCreateTheEnvironmentOnceForConcurrentReaders() throws Exception {
        var count = new AtomicInteger();
        var cache = new SpackCache(new SpackConfig(Map.of("cacheDir", dir.toString()), Map.of())) {
            @Override
            Path createLocalSpackEnv(String spackEnv, String arch) {
                count.incrementAndGet();
                sleep();
                return dir.resolve(spackEnv);
            }
        };
        var env = "foo@" + UUID.randomUUID();

        var paths = DataflowThreads.call(8, () -> cache.getCachePathFor(env, null));

        assertThat(paths).hasSize(8).containsOnly(dir.resolve(env));
        assertThat(count).hasValue(1);
    }

    @Test
    void shouldPropagateTheCreationError() {
        var count = new AtomicInteger();
        var cache = new SpackCache(new SpackConfig(Map.of("cacheDir", dir.toString()), Map.of())) {
            @Override
            Path createLocalSpackEnv(String spackEnv, String arch) {
                count.incrementAndGet();
                throw new IllegalArgumentException("spack boom");
            }
        };
        var env = "foo@" + UUID.randomUUID();

        for( int i = 0; i < 2; i++ ) {
            assertThatThrownBy(() -> DataflowThreads.call(1, () -> cache.getCachePathFor(env, null)))
                .isInstanceOf(IllegalStateException.class)
                .hasCauseInstanceOf(IllegalArgumentException.class)
                .hasRootCauseMessage("spack boom");
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
