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
package nextflow.processor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Implement the v2 task hash computation strategy.
 *
 * This version uses order-independent Map hashing (via entrySet)
 * and checks CacheFunnel before Map in the hash builder.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskHasherV2 extends TaskHasherV1 {

    TaskHasherV2(TaskRun task) {
        super(task)
    }

    @Override
    protected int hashVersion() { 2 }
}
