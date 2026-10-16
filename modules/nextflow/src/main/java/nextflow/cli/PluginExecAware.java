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

package nextflow.cli;

import java.util.List;
import java.util.Map;

/**
 * Define the interface for plugin commands
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public interface PluginExecAware {

    String CMD_SEP = ":";

    int exec(String pluginId, String cmd, List<String> args);

    /**
     * Execute a plugin command with the config resolved by the CLI, so that
     * options such as {@code -c} and {@code -config-ignore-includes}
     * are honoured.
     *
     * <p>Defined as a Java default method so that plugins compiled against a previous
     * version fall back to {@link #exec(String, String, List)}.
     */
    default int exec(String pluginId, String cmd, List<String> args, Map config) {
        return exec(pluginId, cmd, args);
    }

}
