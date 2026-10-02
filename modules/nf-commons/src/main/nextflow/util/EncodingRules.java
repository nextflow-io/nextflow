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

package nextflow.util;

/**
 * The per-key encoding half of a task hash version: how {@link HashBuilder} turns
 * any object into bytes. Frozen per {@code TaskHashSpec} — never mutated in place,
 * because historical rule sets must keep producing historical hashes.
 */
public class EncodingRules {

    /** Behaviour introduced by #6679 "Record types" (2026-03-09) — current master. */
    public static final EncodingRules RECORD_TYPES = new EncodingRules(true, true);

    /** Behaviour before #6679. */
    public static final EncodingRules LEGACY = new EncodingRules(false, false);

    private final boolean orderIndependentMaps;

    private final boolean cacheFunnelFirst;

    public EncodingRules(boolean orderIndependentMaps, boolean cacheFunnelFirst) {
        this.orderIndependentMaps = orderIndependentMaps;
        this.cacheFunnelFirst = cacheFunnelFirst;
    }

    public HashBuilder apply(HashBuilder builder) {
        return builder
            .withOrderIndependentMaps(orderIndependentMaps)
            .withCacheFunnelFirst(cacheFunnelFirst);
    }

    /** Stable text form, contributing to the spec fingerprint. Never reformat this. */
    public String canonicalForm() {
        return "orderIndependentMaps=" + orderIndependentMaps
            + ";cacheFunnelFirst=" + cacheFunnelFirst;
    }
}
