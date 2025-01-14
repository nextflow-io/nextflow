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
 *
 */

package nextflow.processor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.MemoryUnit

/**
 * Models the offload configuration for publishing outputs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class PublishOffloadConfig {
    boolean enable = false
    boolean useFusion = false
    int batchSize = 10
    int maxParallel = 10

    PublishOffloadConfig( Map config){
        if( config.enable ){
            enable = config.enable as boolean
        }
        if( config.useFusion ){
            useFusion = config.useFusion as boolean
        }
        if( config.batchSize ){
            batchSize = config.batchSize as int
            if (batchSize < 0){
                log.warn("Publish offload 'batchSize' property must be bigger than 0. Setting 'batchSize' to 1")
                batchSize = 1
            }
        }
        if( config.maxParallel ){
            maxParallel = config.maxParallel as int
            if (maxParallel < 0){
                log.warn("Publish offload 'maxParallel' property must be bigger than 0. Setting 'maxParallel' to 1")
                maxParallel = 1
            }
        }else{
            // if no maxParallel specified set equal to batchSize
             maxParallel = batchSize
        }
    }




}
