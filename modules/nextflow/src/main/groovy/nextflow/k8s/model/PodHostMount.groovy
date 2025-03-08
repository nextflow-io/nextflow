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

package nextflow.k8s.model

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a K8s pod host mount definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
@ToString(includeNames = true)
@CompileStatic
class PodHostMount {

    private static Set<String> VALID_TYPES = [
        'DirectoryOrCreate', 'Directory', 'FileOrCreate', 'File', 'Socket', 'CharDevice', 'BlockDevice'
    ]

    String hostPath

    String mountPath

    boolean readOnly

    String type // must be empty or a value from VALID_TYPES

    PodHostMount(String host, String container, boolean readOnly=false, String type=null) {
        this.hostPath = host
        this.mountPath = container
        this.readOnly = readOnly
        this.type = type
        validateTypeEnum(type)
    }

    private static validateTypeEnum(String type) {
        if( type && !(type in VALID_TYPES) )
            throw new IllegalArgumentException("K8s invalid hostPath type: $type - must be empty or one of $VALID_TYPES")
    }

}
