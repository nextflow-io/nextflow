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

package nextflow.k8s

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class VolumeClaims {

    @Delegate
    private Map<String,Map<String,String>> target

    VolumeClaims(obj) {
        target = obj instanceof Map ? (Map)obj : Collections.emptyMap()
        // remove trailing slashes
        target.each { k, v -> v.mountPath = sanitize(v.mountPath) }
    }

    static VolumeClaims empty() { new VolumeClaims(Collections.emptyMap()) }

    String findVolumeByPath(String path) {
        target.find { k, v -> path.startsWith(v.mountPath?.toString()) }?.getKey()
    }

    Collection<String> getMountPaths() {
        target.values().collect { it.mountPath }
    }

    String getFirstMount() {
        target.values().iterator().next().mountPath
    }

    static private String sanitize(path) {
        if( !path ) return null
        def result = path.toString().trim()
        while( result.endsWith('/') && result.size()>1 )
            result = result.substring(0,result.size()-1)
        return result
    }
}
