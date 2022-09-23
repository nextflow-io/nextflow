/*
 * Copyright (c) 2019-2022, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.plugin.BasePlugin
import nextflow.cli.PluginExecAware
import org.pf4j.PluginWrapper
/**
 * Nextflow Tower plugin
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerPlugin extends BasePlugin implements PluginExecAware {

    @Delegate private CacheCommand delegate

    TowerPlugin(PluginWrapper wrapper) {
        super(wrapper)
        this.delegate = new CacheCommand()
    }

}
