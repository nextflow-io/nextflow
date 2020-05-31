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

package nextflow.util

import groovy.transform.CompileStatic

/**
 * Register a serializer class in the Kryo serializers list
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface SerializerRegistrant {

    /**
     * Serializer should implement this method adding the
     * serialized class and the serializer class to the
     * map passed as argument
     *
     * @param
     *      serializers The serializer map, where the key element
     *      represent the class to be serialized and the value
     *      the serialization object
     */
    void register(Map<Class,Object> serializers)

}
