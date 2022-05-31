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

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic
/**
 * Compute Engine VM instance provisioning model
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
enum ProvisioningModel {
    /**
     * <pre>
     * Unspecified.
     * </pre>
     *
     * <code>PROVISIONING_MODEL_UNSPECIFIED = 0;</code>
     */
    PROVISIONING_MODEL_UNSPECIFIED(0),
    /**
     * <pre>
     * Standard VM.
     * </pre>
     *
     * <code>STANDARD = 1;</code>
     */
    STANDARD(1),
    /**
     * <pre>
     * SPOT VM.
     * </pre>
     *
     * <code>SPOT = 2;</code>
     */
    SPOT(2),
    /**
     * <pre>
     * Preemptible VM (PVM).
     * Above SPOT VM is the preferable model for preemptible VM instances: the
     * old preemptible VM model (indicated by this field) is the older model,
     * and has been migrated to use the SPOT model as the underlying technology.
     * This old model will still be supported.
     * </pre>
     *
     * <code>PREEMPTIBLE = 3;</code>
     */
    PREEMPTIBLE(3),
    UNRECOGNIZED(-1);
    

    /**
     * @param value The numeric wire value of the corresponding enum entry.
     * @return The enum associated with the given numeric wire value.
     * @deprecated Use {@link #forNumber(int)} instead.
     */
    static ProvisioningModel valueOf(int value) {
        return forNumber(value);
    }

    /**
     * @param value The numeric wire value of the corresponding enum entry.
     * @return The enum associated with the given numeric wire value.
     */
    static ProvisioningModel forNumber(int value) {
        switch (value) {
            case 0: return PROVISIONING_MODEL_UNSPECIFIED;
            case 1: return STANDARD;
            case 2: return SPOT;
            case 3: return PREEMPTIBLE;
            default: return null;
        }
    }

    private final int value;

    private ProvisioningModel(int value) {
        this.value = value;
    }

    @Override
    String toString() {
        return this.name()
    }
}
