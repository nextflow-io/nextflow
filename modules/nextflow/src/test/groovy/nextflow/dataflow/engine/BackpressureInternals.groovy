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
package nextflow.dataflow.engine

import java.util.function.BooleanSupplier

import groovy.transform.CompileStatic

/**
 * Test helpers for backpressure.
 */
@CompileStatic
final class BackpressureInternals {

    private BackpressureInternals() {}

    /**
     * Wait until a condition holds, failing after a deadline.
     */
    static void await(BooleanSupplier condition) throws InterruptedException {
        final deadline = System.currentTimeMillis() + 5_000
        while( !condition.getAsBoolean() ) {
            if( System.currentTimeMillis() > deadline )
                throw new AssertionError("Condition not met in time")
            Thread.sleep(5)
        }
    }
}
