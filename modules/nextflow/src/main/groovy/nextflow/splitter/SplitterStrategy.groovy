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

package nextflow.splitter

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowQueue
/**
 * Interface that splitters object must implements
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface SplitterStrategy {

    /**
     * Abstract splitter method
     *
     * @param object The object to be splitted
     * @param params The map holding the splitting named parameters
     * @param closure An option closure applied to each split entry
     * @return
     */
    def SplitterStrategy target( targetObject )

    SplitterStrategy options( Map options )

    def split()

    long count()

    void each( Closure closure )

    List list()

    DataflowQueue channel()

}
