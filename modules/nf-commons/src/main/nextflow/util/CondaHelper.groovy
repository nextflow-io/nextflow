/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.util

/**
 * Helper methods to handle Conda input packages
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
@Slf4j
@CompileStatic
class CondaHelper {

    boolean containsPip(String str) {
        for( String it : str.tokenize() ) {
            if (it.startsWith("pip:"))
                return true
        }
        return false
    }

    Path condaPipPackagesToCondaFile(String packages, List<String> channels) {

    // DRAFT BACKBONE TO BE REWRITTEN

    //     def pipPackages = []
    //     def condaPackages = []
    //     for( String it : packages.tokenize() ) {
    //         if (it.startsWith("pip:")) {
    //             pipPackages.add(it.substring(4))
    //         } else {
    //             condaPackages.add(it)
    //         }
    //     }
    //     def condaFile = File.createTempFile("conda", ".yml")
    //     condaFile.deleteOnExit()
    //     def condaYml = new Yaml()
    //     def condaMap = [:]
    //     condaMap.channels = channels
    //     condaMap.dependencies = condaPackages
    //     condaMap.pip = pipPackages
    //     condaYml.dump(condaMap, condaFile.newWriter())
    //     return condaFile.toPath()
    // }

}
