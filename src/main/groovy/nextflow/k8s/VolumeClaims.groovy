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
 * Model a kubernetes volumes claims eg
 *
 *   [
 *      volumeClaimName1: [
 *                  mountPath: '/the/mounting/dir'
 *              ],
 *
 *      volumeClaimName2: [
 *                  mountPath: '/the/mounting/dir'
 *              ]
 *   ]
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class VolumeClaims {

    @Delegate
    private Map<String,Map<String,String>> target

    VolumeClaims(obj) {
        this.target = new LinkedHashMap<>()
        if( obj instanceof Map ) {
            addAllSkipExisting(obj)
        }
    }

    static VolumeClaims empty() { new VolumeClaims(Collections.emptyMap()) }

    String findVolumeByPath(String path) {
        target.find { k, v -> path.startsWith(v.mountPath?.toString()) }?.getKey()
    }

    /**
     * @return A collection of all mount paths
     */
    Collection<String> getMountPaths() {
        target.values().collect { it.mountPath }
    }

    /**
     * @return The mount path of the first volume claim entry
     */
    String getFirstMount() {
        target.values().iterator().next().mountPath
    }

    /**
     * @return Collection of all defined volume claim names
     */
    Collection<String> getClaimNames() {
        target.keySet()
    }

    /**
     * Add a volume claim entry
     *
     * @param name The volume claim name
     * @param mount The volume claim mount path
     * @return The object itself
     */
    VolumeClaims add( String name, String mount ) {
        target.put( name, [mountPath: sanitize(mount)] )
        return this
    }

    /**
     * Add a volume claim entry
     *
     * @param nameAndMount A string representing a claim name and mounth path separated by a `:` character
     * @return The object itself
     */
    VolumeClaims add( String nameAndMount ) {
        def (String name, String mount) = nameAndMount.tokenize(':')
        if( !name ) throw new IllegalArgumentException("Missing volume claim name")
        if( !mount ) throw new IllegalArgumentException("Missing volume claim mount path")
        add( name, mount )
    }

    /**
     * Add all volume entries skipping the claim names already defined
     *
     * @param volumes A map of volume claims eg {@code [claim-name1: [mountPath:'/this'], claim-name2: [mountPath: '/that']]}
     * @return The object itself
     */
    VolumeClaims addAllSkipExisting(Map<String,Map> volumes ) {
        for(Map.Entry<String,Map> entry : volumes ) {
            if( target.containsKey(entry.key))
                continue
            target.put(entry.key, entry.value)
            // remove trailing slashes
            entry.value.mountPath = sanitize(entry.value.mountPath)
        }
        return this
    }

    static private String sanitize(path) {
        if( !path ) return null
        def result = path.toString().trim()
        while( result.endsWith('/') && result.size()>1 )
            result = result.substring(0,result.size()-1)
        return result
    }
}
