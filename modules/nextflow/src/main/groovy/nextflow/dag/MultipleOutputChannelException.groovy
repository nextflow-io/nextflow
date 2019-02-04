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

package nextflow.dag

import groovy.transform.PackageScope

/**
 * Exception raised then the same channel is declared as output more
 * than one time
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class MultipleOutputChannelException extends Exception {

    MultipleOutputChannelException( String name, DAG.Vertex duplicate, DAG.Vertex existing ) {
        super(message(name,duplicate,existing))
    }

    @PackageScope
    static message( String name, DAG.Vertex duplicate, DAG.Vertex existing ) {
        if( !name ) {
            return 'Channels cannot be used as output in more than one process or operator'
        }

        if( !duplicate || duplicate.type != DAG.Type.PROCESS ) {
            return "Channel `$name` has been used as an output by more than a process or an operator"
        }

        String message = "Channel `$name` has been used twice as an output by process `${duplicate.label}`"
        if( existing != duplicate )  {
            message += existing?.type == DAG.Type.PROCESS  ? " and process `${existing.label}`" : " and another operator"
        }

        return message
    }
}
