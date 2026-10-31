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

import java.util.Map;

import groovy.lang.Closure;
import groovyx.gpars.scheduler.Pool;

/**
 * Operator helpers for tests in other packages.
 */
public final class TestOperators {

    private TestOperators() {}

    /**
     * Start an operator on the given pool instead of the dataflow pool.
     */
    public static DataflowProcessor start(Pool pool, Map<String, ?> channels, Closure<?> code) {
        return new DataflowOperator(pool, channels, code).start();
    }
}
