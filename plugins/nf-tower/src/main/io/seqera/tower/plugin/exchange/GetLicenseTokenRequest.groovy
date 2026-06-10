/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.tower.plugin.exchange


import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Models a REST request to obtain a license-scoped JWT token from Platform
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
@EqualsAndHashCode
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class GetLicenseTokenRequest {

    /** The product code */
    String product

    /** The product version */
    String version

    /**
     * The Platform workflow ID associated with this request
     */
    String workflowId

    /**
     * The Platform workspace ID associated with this request
     */
    String workspaceId

    /**
     * @return a Map representation of the request
     */
    Map<String, String> toMap() {
        final map = new HashMap<String, String>()
        map.product = this.product
        map.version = this.version
        map.workflowId = this.workflowId
        map.workspaceId = this.workspaceId
        return map
    }
}
