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
 * Model a K8s pod emptyDir mount
 *
 * See also https://kubernetes.io/docs/concepts/storage/volumes/#emptydir
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class PodMountEmptyDir {

    String mountPath

    Map emptyDir

    PodMountEmptyDir( Map emptyDir, String mountPath ) {
        assert mountPath

        this.emptyDir = emptyDir
        this.mountPath = mountPath
    }

    PodMountEmptyDir( Map entry ) {
        this(entry.emptyDir as Map, entry.mountPath as String)
    }

}
