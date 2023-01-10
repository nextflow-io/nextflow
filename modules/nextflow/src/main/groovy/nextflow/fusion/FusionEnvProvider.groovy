/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

<<<<<<< HEAD:modules/nextflow/src/main/groovy/nextflow/k8s/model/PodMountEmptyDir.groovy
package nextflow.k8s.model

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
=======
package nextflow.fusion
>>>>>>> master:modules/nextflow/src/main/groovy/nextflow/fusion/FusionEnvProvider.groovy

import nextflow.Global
import nextflow.SysEnv
import nextflow.plugin.Plugins
/**
<<<<<<< HEAD:modules/nextflow/src/main/groovy/nextflow/k8s/model/PodMountEmptyDir.groovy
 * Model a K8s pod emptyDir mount
=======
 * Provider strategy for {@link FusionEnv}
>>>>>>> master:modules/nextflow/src/main/groovy/nextflow/fusion/FusionEnvProvider.groovy
 *
 * See also https://kubernetes.io/docs/concepts/storage/volumes/#emptydir
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
<<<<<<< HEAD:modules/nextflow/src/main/groovy/nextflow/k8s/model/PodMountEmptyDir.groovy
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
=======
class FusionEnvProvider {

    Map<String,String> getEnvironment(String scheme) {
        final config = new FusionConfig(Global.config?.fusion as Map ?: Collections.emptyMap(), SysEnv.get())
        final list = Plugins.getExtensions(FusionEnv)
        final result = new HashMap<String,String>()
        for( FusionEnv it : list ) {
            final env = it.getEnvironment(scheme,config)
            if( env ) result.putAll(env)
        }
        return result
>>>>>>> master:modules/nextflow/src/main/groovy/nextflow/fusion/FusionEnvProvider.groovy
    }

}