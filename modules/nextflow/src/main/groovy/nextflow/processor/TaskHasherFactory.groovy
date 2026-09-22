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

import org.pf4j.ExtensionPoint

/**
 * Plugin extension point through which a plugin supplies the {@link TaskHasher} for the tasks
 * of this run, i.e. decides what the task cache key covers.
 *
 * <p>{@link TaskProcessor} resolves the registered factories once per processor, in
 * {@link nextflow.plugin.Priority} order, and asks each one for a hasher; the first non-null
 * answer wins. With no factory registered, or all of them abstaining, the default
 * {@link TaskHasher} is used, so a run without such a plugin hashes exactly as before.
 *
 * <p>A factory is expected to abstain (return {@code null}) whenever its cache is not the one in
 * use for the session: the plugin registry is process-wide, while the choice of hasher belongs to
 * the run.
 *
 * <p><b>Implementations must be stateless / thread-safe.</b> Extensions are resolved through pf4j's
 * {@code SingletonExtensionFactory}, so a single instance is shared by every processor of every run
 * in the JVM, and {@link #create} is called concurrently from every operator thread plus the retry
 * executor.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskHasherFactory extends ExtensionPoint {

    /**
     * Create the hasher for the given task.
     *
     * @param task The task to be hashed.
     * @return A {@link TaskHasher} for {@code task}, or {@code null} to abstain, in which case the
     *      next factory is asked and, failing all, the default {@link TaskHasher} is used.
     */
    TaskHasher create(TaskRun task)

}
