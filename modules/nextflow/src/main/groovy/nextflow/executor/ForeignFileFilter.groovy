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

package nextflow.executor

import java.nio.file.Path

import org.pf4j.ExtensionPoint

/**
 * Plugin extension point that lets a plugin declare a foreign input file readable by every task
 * at the same absolute path, so it is not staged into the work directory.
 *
 * <p>It can only veto the executor's decision: it is asked only for files the executor
 * considers foreign, and it cannot make a local file foreign.
 *
 * <p><b>Implementations must be thread-safe</b>: a single instance is shared by every task of
 * every run in the JVM. The run context is available from the given executor.
 */
interface ForeignFileFilter extends ExtensionPoint {

    /**
     * @return {@code true} when every task can read {@code path} at the same absolute path
     */
    boolean isTaskAccessible(Path path, Executor executor)

}
