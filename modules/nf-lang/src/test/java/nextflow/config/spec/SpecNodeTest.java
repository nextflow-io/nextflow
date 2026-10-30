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
package nextflow.config.spec;

import java.util.List;
import java.util.Map;

import nextflow.script.types.Duration;
import nextflow.script.types.MemoryUnit;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Specifies the config options derived for the built-in config scopes.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class SpecNodeTest {

    private static SpecNode.Option option(String... names) {
        return SpecNode.ROOT.getOption(List.of(names));
    }

    @ParameterizedTest
    @ValueSource(strings = { "arch", "cache", "clusterOptions", "cpus", "disk", "publishDir" })
    void describeEveryProcessDirectiveEvenIfOverloaded(String name) {
        var node = SpecNode.ROOT.getChild(List.of("process", name));
        assertThat(node.description()).isNotEmpty();
    }

    @Test
    void inferProcessConfigOptionTypesFromDirectives() {
        assertThat(option("process", "clusterOptions").types()).containsExactlyInAnyOrder(String.class, List.class, String[].class);
        assertThat(option("process", "cpus").types()).containsExactly(Integer.class);
        assertThat(option("process", "errorStrategy").types()).containsExactly(String.class);
        assertThat(option("process", "executor").types()).containsExactly(String.class);
        assertThat(option("process", "ext").types()).containsExactly(Map.class);
        assertThat(option("process", "memory").types()).containsExactly(MemoryUnit.class);
        assertThat(option("process", "publishDir").types()).containsExactlyInAnyOrder(Map.class, String.class, List.class);
        assertThat(option("process", "resourceLimits").types()).containsExactly(Map.class);
        assertThat(option("process", "secret").types()).containsExactly(List.class);
        assertThat(option("process", "time").types()).containsExactly(Duration.class);
    }

    @Test
    void defineNextflowRetryPolicyOptions() {
        assertThat(option("nextflow", "retryPolicy", "delay").types()).containsExactly(Duration.class);
        assertThat(option("nextflow", "retryPolicy", "maxDelay").types()).containsExactly(Duration.class);
        assertThat(option("nextflow", "retryPolicy", "maxAttempts").types()).containsExactly(Integer.class);
        assertThat(option("nextflow", "retryPolicy", "jitter").types()).containsExactly(Double.class);
        assertThat(option("nextflow", "retryPolicy", "multiplier").types()).containsExactly(Double.class);
        assertThat(option("nextflow", "retryPolicy", "delay").description()).isNotEmpty();
    }

}
