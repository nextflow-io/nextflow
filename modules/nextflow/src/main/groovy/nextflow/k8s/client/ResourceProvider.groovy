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

package nextflow.k8s.client

import java.nio.file.Path

import nextflow.k8s.client.K8sRequestDelegate
import nextflow.k8s.client.K8sResponseJson

/**
 * Models the API operations for a K8s compute resource (Pod, Job, etc).
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
interface ResourceProvider {

    ResourceProvider withDelegate(K8sRequestDelegate delegate)

    K8sResponseJson create(Map spec, Path saveYamlPath)

    K8sResponseJson delete(String name)

    K8sResponseJson list(boolean allNamespaces)

    InputStream fetchLog(Map params, String name)

    Map getState(String name)

}

