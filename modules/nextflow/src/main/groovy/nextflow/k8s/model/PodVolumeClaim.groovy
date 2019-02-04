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

package nextflow.k8s.model

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a K8s pod persistent volume claim mount
 *
 * See https://kubernetes.io/docs/tasks/configure-pod-container/configure-persistent-volume-storage/#create-a-persistentvolumeclaim
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class PodVolumeClaim {

    String claimName

    String mountPath

    String subPath

    PodVolumeClaim(String name, String mount, String subPath=null) {
        assert name
        assert mount
        this.claimName = name
        this.mountPath = sanitize(mount)
        this.subPath = subPath
        validate(mountPath)
    }

    PodVolumeClaim(Map entry) {
        assert entry.volumeClaim
        assert entry.mountPath
        this.claimName = entry.volumeClaim
        this.mountPath = sanitize(entry.mountPath)
        this.subPath = entry.subPath
        validate(mountPath)
    }

    private static validate(String path) {
        if( !path.startsWith('/') )
            throw new IllegalArgumentException("K8s volume claim path must be an absolute path: $path")
    }

    static private String sanitize(path) {
        if( !path ) return null
        def result = path.toString().trim()
        while( result.endsWith('/') && result.size()>1 )
            result = result.substring(0,result.size()-1)
        return result
    }

}
