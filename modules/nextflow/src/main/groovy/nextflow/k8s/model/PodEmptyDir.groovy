/*
 * Copyright 2019, WuxiNextcode
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
import nextflow.util.MemoryUnit

/**
 * Model a K8s pod emptyDir definition
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 */
@EqualsAndHashCode
@ToString(includeNames = true)
@CompileStatic
class PodEmptyDir {

    enum PodEmptyDirType {
        Disk, Memory
    }

    String name
    String mountPath
    PodEmptyDirType type
    MemoryUnit sizeLimit

    PodEmptyDir(String name, String mountPath,PodEmptyDirType type, String sizeLimit = null) {
        this.name = name
        this.type = type
        this.mountPath = mountPath
        if(sizeLimit) {
            this.sizeLimit = new MemoryUnit(sizeLimit)
        }
    }

    PodEmptyDir( Map entry ) {
            this(entry.emptyDir as String, entry.mountPath as String, PodEmptyDirType.valueOf(entry.type as String), entry.sizeLimit as String )
    }
}