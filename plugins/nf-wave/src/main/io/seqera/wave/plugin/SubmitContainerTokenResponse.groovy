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

package io.seqera.wave.plugin

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Model a response for an augmented container
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class SubmitContainerTokenResponse {

    /**
     * Unique Id for this request
     */
    String requestId

    /**
     * A unique authorization token assigned to this request
     */
    String containerToken

    /**
     * The fully qualified wave container name to be used
     */
    String targetImage

    /**
     * The source container image that originated this request
     */
    String containerImage

    /**
     * The ID of the build associated with this request or null of the image already exists
     */
    String buildId

    /**
     * Whenever it's a cached build image. Only supported by API version v1alpha2
     */
    Boolean cached

    /**
     * When the result is a freeze container. Version v1alpha2 as later.
     */
    Boolean freeze;

    /**
     * When the result is a mirror container. Version v1alpha2 as later.
     */
    Boolean mirror

    /**
     * The id of the security scan associated with this container
     */
    String scanId

    /**
     * Whenever the container has been provisioned successfully or not. If false
     * the current status needs the be check via container status API
     */
    Boolean succeeded

}
