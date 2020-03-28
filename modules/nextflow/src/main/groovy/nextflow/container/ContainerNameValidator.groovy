/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.container

import java.util.regex.Pattern

/**
 * Helper class to validate a container image name
 *
 * Credits
 * https://github.com/jenkinsci/docker-plugin/blob/4733aee58fbcb11b050cfdd9f97fb0980e19bb28/src/main/java/com/nirima/jenkins/plugins/docker/builder/DockerBuilderPublisher.java
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerNameValidator {

    /**
     * The docker spec says "<i>A tag name may contain lowercase and uppercase
     * characters, digits, underscores, periods and dashes. A tag name may not
     * start with a period or a dash and may contain a maximum of 128
     * characters</i>". This is a simplified version of that specification.
     */
    private static final String TAG_REGEX = "[a-zA-Z0-9-_.]+";

    /**
     * The docker spec says "<i>Name components may contain lowercase
     * characters, digits and separators. A separator is defined as a period,
     * one or two underscores, or one or more dashes. A name component may not
     * start or end with a separator</i>". This is a simplified version of that
     * specification.
     */
    private static final String NAME_COMPONENT_REGEX = "[a-z0-9-_.]+";

    /**
     * The docker spec says "<i>The (registry) hostname must comply with
     * standard DNS rules, but may not contain underscores. If a hostname is
     * present, it may optionally be followed by a port number in the format
     * :8080</i>". This is a simplified version of that specification.
     */
    private static final String REGISTRY_HOSTNAME_REGEX = "[a-zA-Z0-9-.]+(:[0-9]+)?";

    /**
     * The docker spec says "<i>An image name is made up of slash-separated name
     * components, optionally prefixed by a registry hostname</i>".
     */
    private static final String IMAGE_NAME_REGEX = "(" + REGISTRY_HOSTNAME_REGEX + "/)?" + NAME_COMPONENT_REGEX + "(/" + NAME_COMPONENT_REGEX + ")*";

    /**
     * A regex matching IMAGE[:TAG] (from the "docker tag" command) where IMAGE
     * matches {@link #IMAGE_NAME_REGEX} and TAG matches {@link #TAG_REGEX}.
     */
    private static final String VALID_REPO_REGEX = "^" + IMAGE_NAME_REGEX + "(:" + TAG_REGEX + ")?\$"

    /** Compiled version of {@link #VALID_REPO_REGEX}. */
    private static final Pattern VALID_REPO_PATTERN = Pattern.compile(VALID_REPO_REGEX)


    static boolean isValidImageName(String image) {
        VALID_REPO_PATTERN.matcher(image).matches()
    }

}
