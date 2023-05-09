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

package nextflow.k8s.model

import nextflow.executor.res.AcceleratorResource
import nextflow.util.MemoryUnit

/**
 * Interface for building specifications for compute resources
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
interface ResourceSpecBuilder {

    ResourceSpecBuilder withAccelerator(AcceleratorResource acc)

    ResourceSpecBuilder withActiveDeadline(int seconds)

    ResourceSpecBuilder withAnnotations(Map annotations)

    ResourceSpecBuilder withArgs(args)

    ResourceSpecBuilder withCommand(cmd)

    ResourceSpecBuilder withCpus(Integer cpus)

    ResourceSpecBuilder withDisk(MemoryUnit disk)

    ResourceSpecBuilder withEnv(PodEnv env)

    ResourceSpecBuilder withHostMount(String host, String mount)

    ResourceSpecBuilder withImageName(String name)

    ResourceSpecBuilder withLabels(Map labels)

    ResourceSpecBuilder withMemory(MemoryUnit mem)

    ResourceSpecBuilder withName(String name)

    ResourceSpecBuilder withNamespace(String name)

    ResourceSpecBuilder withPodOptions(PodOptions opts)

    ResourceSpecBuilder withPrivileged(boolean value)

    ResourceSpecBuilder withServiceAccount(String name)

    Map build()

}
