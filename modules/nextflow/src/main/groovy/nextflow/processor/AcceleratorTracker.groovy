/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.processor

import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.util.TrackingSemaphore

/**
 * Specialized semaphore that keeps track of accelerators by
 * id. The id can be an integer or a UUID.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class AcceleratorTracker {

    private static final List<String> DEVICE_ENV_NAMES = [
        'CUDA_VISIBLE_DEVICES',
        'HIP_VISIBLE_DEVICES',
        'ROCR_VISIBLE_DEVICES'
    ]

    static AcceleratorTracker create() {
        return create(SysEnv.get())
    }

    static AcceleratorTracker create(Map<String,String> env) {
        return DEVICE_ENV_NAMES.stream()
            .filter(name -> env.containsKey(name))
            .map((name) -> {
                final ids = env.get(name).tokenize(',')
                return new AcceleratorTracker(name, ids)
            })
            .findFirst().orElse(new AcceleratorTracker())
    }

    private final String name
    private final TrackingSemaphore semaphore

    private AcceleratorTracker(String name, List<String> ids) {
        this.name = name
        this.semaphore = new TrackingSemaphore(ids)
    }

    private AcceleratorTracker() {
        this(null, [])
    }

    String name() {
        return name
    }

    int total() {
        return semaphore.totalPermits()
    }

    int available() {
        return semaphore.availablePermits()
    }

    List<String> acquire(int permits) {
        return semaphore.acquire(permits)
    }

    void release(List<String> ids) {
        semaphore.release(ids)
    }

}
