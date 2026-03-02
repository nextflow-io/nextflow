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

import nextflow.Session
import org.pf4j.ExtensionPoint
/**
 * Factory class creating {@link TraceObserver} instances
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Deprecated
interface TraceObserverFactory extends ExtensionPoint {

    /**
     * Register the observer on the current session object
     *
     * @param session The current {@link nextflow.Session} instance
     * @return One or more instances of {@link TraceObserver} objects
     */
    Collection<TraceObserver> create(Session session)

}
