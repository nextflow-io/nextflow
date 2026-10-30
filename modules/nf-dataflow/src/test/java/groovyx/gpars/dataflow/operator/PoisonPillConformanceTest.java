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
package groovyx.gpars.dataflow.operator;

import java.util.List;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Pins the {@link PoisonPill} identity that {@code Channel.STOP} and hundreds of
 * {@code == Channel.STOP} assertions rely on.
 */
class PoisonPillConformanceTest {

    @Test
    void shouldBeAControlMessageSingleton() {
        var pill = PoisonPill.getInstance();
        assertThat(pill).isSameAs(PoisonPill.getInstance());
        assertThat(pill).isInstanceOf(ControlMessage.class);
        assertThat(pill.isImmediate()).isFalse();
        assertThat(pill).isEqualTo(PoisonPill.getInstance());
    }

    @Test
    void shouldNotBeEqualToOtherValues() {
        var pill = PoisonPill.getInstance();
        assertThat(pill).isNotEqualTo(null);
        assertThat(pill).isNotEqualTo("STOP");
        assertThat(pill).isNotSameAs(PoisonPill.getImmediateInstance());
        assertThat(List.of(1, 2, pill).indexOf(pill)).isEqualTo(2);
    }

}
