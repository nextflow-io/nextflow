/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.script.params

import nextflow.extension.CH
import nextflow.script.ProcessConfig
/**
 * Model a process default input parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
final class DefaultInParam extends ValueInParam {

    @Override
    String getTypeName() { 'default' }

    DefaultInParam(ProcessConfig config) {
        super(config)
        // This must be a dataflow queue channel to which
        // just a value is bound -- No STOP value has to be emitted
        // because this channel is used to control to process termination
        // See TaskProcessor.BaseProcessInterceptor#messageArrived
        final channel = CH.queue()
        channel.bind(Boolean.TRUE)
        setFrom(channel)
        bind('$')
    }
}
