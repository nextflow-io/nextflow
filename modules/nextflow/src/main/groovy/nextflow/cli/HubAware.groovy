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
 */

package nextflow.cli

import com.beust.jcommander.Parameter
import groovy.transform.CompileStatic
import nextflow.scm.HubOptions
/**
  * Command line parameters for commands that interact with a git repository provider i.e. GitHub or BitBucket
  *
  * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
  */
@CompileStatic
interface HubAware {

    @Parameter(names=['-hub'], description = "Service hub where the project is hosted")
    String hubProvider

    @Parameter(names='-user', description = 'Private repository user name')
    String hubUser

    default HubOptions toHubOptions() {
        return new HubOptions(hubProvider, hubUser)
    }
}
