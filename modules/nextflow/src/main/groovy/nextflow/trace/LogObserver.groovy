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

package nextflow.trace

/**
 * Common interface for log observer implementations (ANSI and Agent modes).
 *
 * Only one log observer is active at a time. This interface defines the
 * shared contract used by {@link nextflow.Session} and
 * {@link nextflow.util.LoggerHelper}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
// declared as a trait to keep binary compatibility with plugins compiled with Groovy 4
trait LogObserver {

    abstract void appendInfo(String message)

    abstract void appendWarning(String message)

    abstract void appendError(String message)

    abstract void forceTermination()

    void appendSticky(String message) { }

    boolean getStarted() { true }

    boolean getStopped() { false }
}
