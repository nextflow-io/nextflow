/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
 * Model a K8s pod ConfigMap mount
 *
 * See also https://kubernetes.io/docs/tasks/configure-pod-container/configure-pod-configmap/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class PodMountConfig {

    String mountPath

    String fileName

    String configName

    String configKey

    PodMountConfig( String config, String mount ) {
        assert config
        assert mount

        final path = Paths.get(mount)
        final tokens = config.tokenize('/')
        configName = tokens[0].trim()
        configKey = tokens.size()>1 ? tokens[1].trim() : null
        if( configKey ) {
            mountPath = path.parent.toString()
            fileName = path.fileName.toString()
        }
        else {
            mountPath = path.toString()
        }
    }

    PodMountConfig( Map entry ) {
        this(entry.config as String, entry.mountPath as String)
    }

}
