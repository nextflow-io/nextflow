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

import nextflow.processor.TaskRun
import org.pf4j.ExtensionPoint

/**
 * Extension point through which a plugin supplies the task hash spec for a task.
 *
 * Plugins contribute a spec rather than a hasher so that every hash goes through
 * BaseTaskHasher: per-key digests and spec fingerprints then hold for plugin-driven
 * runs too, and a plugin cannot duplicate the key list to vary two of its entries.
 */
interface TaskHashSpecFactory extends ExtensionPoint {

    /**
     * @return the spec for {@code task}, or {@code null} to abstain, in which case the
     *      next factory is asked and, failing all, the default spec is used.
     */
    TaskHashSpec create(TaskRun task)
}
