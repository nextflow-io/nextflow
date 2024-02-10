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

import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a K8s pod CSI ephemeral volume mount
 *
 * See also https://kubernetes.io/docs/concepts/storage/ephemeral-volumes/#csi-ephemeral-volumes
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class PodMountCsiEphemeral {

    String mountPath

    Map csi

    PodMountCsiEphemeral( Map csi, String mountPath ) {
        assert csi
        assert mountPath

        this.csi = csi
        this.mountPath = mountPath
    }

    PodMountCsiEphemeral( Map entry ) {
        this(entry.csi as Map, entry.mountPath as String)
    }

}
