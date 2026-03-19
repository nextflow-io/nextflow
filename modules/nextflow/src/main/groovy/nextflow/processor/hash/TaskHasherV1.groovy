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
package nextflow.processor.hash

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.HashBuilder
/**
 * V1 task hash computation strategy.
 *
 * This is the original hashing behavior before the record types change.
 * Maps are hashed by values only (order-dependent) and CacheFunnel
 * is checked after Map and SerializableMarker.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskHasherV1 extends AbstractTaskHasher {

    TaskHasherV1(TaskRun task) {
        super(task)
    }

    @Override
    protected HashBuilder createHashBuilder() {
        return new HashBuilder()
            .withOrderIndependentMaps(false)
            .withCacheFunnelFirst(false)
    }
}
