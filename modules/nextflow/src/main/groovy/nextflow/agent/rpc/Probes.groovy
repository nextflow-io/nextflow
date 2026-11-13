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

package nextflow.agent.rpc

/**
 * Everything about the DRIVER HOST the ladder needs to observe, behind an interface so the rows
 * can be driven deterministically from a test. Nothing here takes an argument that varies per
 * agent definition: each answer is a property of the host, which is what makes the memoization
 * below sound.
 */
interface Probes {
    /**
     * The local address the kernel selects for the default route, or {@code null} when the host
     * has no route at all. A routing-table lookup -- no packet is sent.
     */
    String outboundAddress()

    /** Every non-loopback, non-link-local interface address, for the multi-homed warning. */
    List<String> interfaceAddresses()
}
