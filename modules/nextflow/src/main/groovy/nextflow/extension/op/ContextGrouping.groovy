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
 *
 */

package nextflow.extension.op

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.prov.OperatorRun

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ContextGrouping implements OpContext {

    private final Map<String,OperatorRun> holder = new ConcurrentHashMap<>(1);

    ContextGrouping(){
        holder.put('run', new OperatorRun())
    }

    @Override
    OperatorRun allocateRun() {
        final result = holder.get('run')
        log.trace "+ AllocateRun=$result"
        return result
    }
    
    @Override
    OperatorRun getOperatorRun() {
        final result = holder.get('run')
        log.trace "+ GetOperatorRun=$result"
        holder.put('run', new OperatorRun())
        return result
    }

}
