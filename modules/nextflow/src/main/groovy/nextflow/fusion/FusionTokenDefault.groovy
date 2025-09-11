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

package nextflow.fusion

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ReportWarningException
import org.pf4j.Extension

/**
 * Implements default for Fusion token environment provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Extension
@CompileStatic
class FusionTokenDefault implements FusionToken {

    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {
        final msg = "Unable to validate Fusion license - Make sure to have added in your config 'tower.accessToken' or have defined the variable TOWER_ACCESS_TOKEN in your launch environment"
        throw new ReportWarningException(msg, 'getFusionLicenseException')
    }
}
