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

package nextflow.processor

/**
 * Defines the contract for {@link TaskHandler} classes that need
 * to aggregate multiple operation to optimise the interaction with
 * a remote execution system
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait BatchHandler<K,V> {

    /**
     * Defines a {@link BatchContext} for the this handler
     *
     * @param context The {@link BatchContext} instance
     */
    abstract void batch( BatchContext<K,V> context )

}
