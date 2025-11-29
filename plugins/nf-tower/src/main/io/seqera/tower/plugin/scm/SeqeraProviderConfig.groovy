/*
 * Copyright 2013-2025, Seqera Labs
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

package io.seqera.tower.plugin.scm

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.BaseCommandImpl
import io.seqera.tower.plugin.TowerConfig
import io.seqera.tower.plugin.TowerRetryPolicy
import nextflow.Const
import nextflow.Global
import nextflow.SysEnv
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.platform.PlatformHelper
import nextflow.scm.ProviderConfig

/**
 * Implements a provider config for Seqera Platform data-links git-remote repositories
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class SeqeraProviderConfig extends ProviderConfig {

    private TowerConfig config

    SeqeraProviderConfig(String name, Map values) {
        super(name, [server: "seqera://$name"] + values)
        setValues(values)
    }

    SeqeraProviderConfig(String name) {
        super(name, [platform: 'seqera', server: "seqera://$name"])
        setValues()
    }

    private void setValues(Map values = Map.of()) {
        // Get tower config from session if exists
        def towerConfig = getDefaultValues()


        // Merge with SCM values
        config = new TowerConfig(new HashMap(towerConfig + values), SysEnv.get())

        if (!config.accessToken) {
            throw new AbortOperationException("Seqera Platform access token not configured. Set TOWER_ACCESS_TOKEN environment variable or configure it in nextflow.config")
        }
    }

    String getEndpoint() {
        this.config.endpoint
    }

    String getAccessToken() {
        this.config.accessToken
    }

    String getWorkspaceId() {
        this.config.workspaceId
    }

    TowerRetryPolicy getRetryPolicy() {
        this.config.retryPolicy
    }

    @Override
    protected String resolveProjectName(String path) {
        log.debug("Resolving project name from $path")
        if (!server.startsWith('seqera://')) {
            throw new AbortOperationException("Seqera project server doesn't start with seqera://")
        }
        return "${server.substring('seqera://'.size())}/$path"
    }

    /**
     * Get the default values from the session or the global config
     */
    private Map getDefaultValues() {
        Map defaultValues = Global.session?.config?.tower as Map
        if(!defaultValues)
            defaultValues = new ConfigBuilder().setHomeDir(Const.APP_HOME_DIR).setCurrentDir(Const.APP_HOME_DIR).buildConfigObject()?.tower as Map
        return defaultValues ?: Map.of()
    }
}
