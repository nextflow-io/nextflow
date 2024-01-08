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
 *
 */

package nextflow.configv2

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessOpts {

    String accelerator
    String afterScript
    String arch
    String beforeScript
    String cache
    String conda
    String cpus
    String container
    String containerOptions
    String cleanup
    String clusterOptions
    String debug
    String disk
    String errorStrategy
    String executor
    Map ext
    Boolean fair
    String machineType
    String queue
    String label
    String maxSubmitAwait
    String maxErrors
    String maxForks
    String maxRetries
    String memory
    String module
    String penv
    String pod
    String publishDir
    String scratch
    String shell
    String spack
    String tag
    String time
    String stageInMode
    String stageOutMode
    String resourceLabels
}
