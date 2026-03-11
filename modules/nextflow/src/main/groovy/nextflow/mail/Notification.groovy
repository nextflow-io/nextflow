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
package nextflow.mail

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

@ScopeName("notification")
@Description("""
    The `notification` scope controls the automatic sending of an email notification on workflow completion.
""")
@CompileStatic
@EqualsAndHashCode
class Notification implements ConfigScope {

    @ConfigOption
    @Description("""
        Map of variables that can be used in the template file.
    """)
    final Map attributes

    @ConfigOption
    @Description("""
        Send an email notification when the workflow execution completes (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        Sender address for the email notification.
    """)
    final String from

    @ConfigOption
    @Description("""
        Path of a template file containing the contents of the email notification.
    """)
    final Object template

    @ConfigOption
    @Description("""
        Recipient address for the email notification. Multiple addresses can be specified as a comma-separated list.
    """)
    final String to

    /* required by extension point -- do not remove */
    Notification() {}

    Notification(Map opts) {
        attributes = opts.attributes as Map
        enabled = opts.enabled as boolean
        from = opts.from
        template = opts.template
        to = opts.to
    }

}
