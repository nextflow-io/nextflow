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

package nextflow.extension

import groovy.transform.CompileStatic

/**
 * Hack to capture the names of referenced variables into a closure
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CaptureProperties {

    List<String> names = new ArrayList<>(10)

    Object propertyMissing(String name) {
        if( names.contains(name) )
            throw new IllegalArgumentException("Duplicate channel definition: $name")
        names.add(name)
        return null
    }

    static List<String> capture(Closure holder) {
        final recorder = new CaptureProperties()
        holder.setResolveStrategy( Closure.DELEGATE_ONLY )
        holder.delegate = recorder
        holder.call()
        return recorder.names
    }

}
