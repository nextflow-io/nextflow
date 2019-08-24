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

import groovyx.gpars.dataflow.DataflowWriteChannel

/**
 * Model a process generic input parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

interface OutParam extends Cloneable {

    interface Mode {  }

    /**
     * @return The parameter name getter
     */
    String getName()

    /**
     * Defines the channel to which bind the output(s) in the script context
     *
     * @param value It can be a string representing a channel variable name in the script context. If
     *      the variable does not exist it creates a {@code DataflowVariable} in the script with that name.
     *      If the specified {@code value} is a {@code DataflowWriteChannel} object, use this object
     *      as the output channel
     * @return
     */
    OutParam into( def value )

    /**
     * @return The output channel instance
     */
    DataflowWriteChannel getOutChannel()

    @Deprecated List<DataflowWriteChannel> getOutChannels()

    short getIndex()

    @Deprecated Mode getMode()

    String getChannelEmitName()

}
