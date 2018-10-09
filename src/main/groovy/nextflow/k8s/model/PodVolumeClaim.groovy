/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
