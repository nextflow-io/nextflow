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

import groovy.transform.CompileStatic

/**
 * The address a containerized agent task dials to reach the driver's RPC broker, and where it
 * came from.
 *
 * <p>Carrying the SOURCE alongside the host is not decoration: the expensive failure mode of the
 * whole feature is a plausible-but-unroutable address, and the only cheap way to tell an operator
 * that the driver GUESSED is to name the rung that answered. {@link #warnings} carries the cases
 * the ladder resolves but is not certain about (a multi-homed driver, a containerized driver whose
 * task may land on another docker network), which the broker prints at the registration line.
 *
 * <p>An unresolved instance is an ERROR ROW, not a null: {@link #error} is the message the
 * pre-ignition guard rejects the run with, already naming what was tried and what to set.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentRpcHost {

    final String host
    final String source
    final List<String> warnings
    /** The error-row code -- {@code E1}..{@code E7} -- or {@code null} when resolved. */
    final String code
    /** The rejection message, or {@code null} when resolved. */
    final String error

    private AgentRpcHost(String host, String source, List<String> warnings, String code, String error) {
        this.host = host
        this.source = source
        this.warnings = warnings != null ? Collections.unmodifiableList(warnings) : Collections.<String>emptyList()
        this.code = code
        this.error = error
    }

    static AgentRpcHost of(String host, String source, List<String> warnings = null) {
        return new AgentRpcHost(host, source, warnings, null, null)
    }

    static AgentRpcHost error(String code, String message) {
        return new AgentRpcHost(null, null, null, code, message)
    }

    boolean isResolved() { host != null }

    String getHost() { host }

    String getSource() { source }

    List<String> getWarnings() { warnings }

    String getCode() { code }

    String getError() { error }

    /** The address with the rung that produced it, e.g. {@code 10.0.3.17 (inferred from default route)}. */
    String describe() {
        return resolved ? "${host} (${source})".toString() : "unresolved (${code}: ${error})".toString()
    }

    @Override
    String toString() { describe() }
}
