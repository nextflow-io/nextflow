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

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.plugin.cli.PluginAbstractExec
/**
 * Implements nextflow cache and restore commands
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CacheCommand implements PluginAbstractExec {

    List<String> getCommands() { ['cache-backup', 'cache-restore'] }

    @Override
    int exec(String cmd, List<String> args) {

        if( cmd == 'cache-backup') {
            cacheBackup()
        }
        if( cmd == 'cache-restore' )
            cacheRestore()
        return 0
    }

    protected void cacheBackup() {
        log.debug "Running Nextflow cache backup"
        final env = System.getenv()
        if( !env.get('NXF_CLOUDCACHE_PATH') ) {
            // legacy cache manager
            new CacheManager(env).saveCacheFiles()
        }
        else {
            new LogsHandler(getSession(), env).saveFiles()
        }
    }

    protected void cacheRestore() {
        final env = System.getenv()
        if( !env.get('NXF_CLOUDCACHE_PATH') ) {
            log.debug "Running Nextflow cache restore"
            // legacy cache manager
            new CacheManager(env).restoreCacheFiles()
        }
        else {
            // this command is only kept for backward compatibility
            log.debug "Running Nextflow cache restore - DO NOTHING"
        }
    }

}
