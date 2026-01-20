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

package nextflow.cli

import org.pf4j.ExtensionPoint

/**
 * Extension point interface for the `auth` command.
 *
 * @see io.seqera.tower.plugin.auth.AuthCommandImpl
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
interface AuthCommand extends ExtensionPoint {
    /**
     * Authenticates with Seqera Platform and saves credentials to config.
     *
     * @param url the Seqera Platform API endpoint URL (null for default)
     */
    void login(String url)

    /**
     * Revokes access token and removes authentication from local config.
     */
    void logout()

    /**
     * Configures Seqera Platform settings (workspace, monitoring, compute environment).
     */
    void config()

    /**
     * Displays current authentication status and configuration sources.
     */
    void status()
}
