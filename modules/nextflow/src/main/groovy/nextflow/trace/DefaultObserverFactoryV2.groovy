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

package nextflow.trace

import nextflow.Session

/**
 * Creates Nextflow observers
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultObserverFactoryV2 implements TraceObserverFactoryV2 {

    private Session session

    @Override
    Collection<TraceObserver> create(Session session) {
        this.session = session

        final result = new ArrayList(10)
        createAnsiLogObserver(result)
        return result
    }

    protected void createAnsiLogObserver(Collection<TraceObserver> result) {
        if( session.ansiLog ) {
            session.ansiLogObserver = new AnsiLogObserver()
            result << session.ansiLogObserver
        }
    }

}
