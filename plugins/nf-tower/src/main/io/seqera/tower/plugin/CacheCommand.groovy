/*
 * Copyright 2020-2022, Seqera Labs
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

import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cli.Launcher
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.cli.PluginExecAware

/**
 * Implements nextflow cache and restore commands
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CacheCommand implements PluginExecAware {

    @Override
    int exec(Launcher launcher, List<String> args) {
        if( !args )
            throw new AbortOperationException("Missing plugin command")
        final cmd = args.pop()
        if( cmd !in ['cache-backup','cache-restore'])
            throw new AbortOperationException("Unknown plugin command: ${cmd}")

        // create the config
        final config = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(Paths.get('.'))
                .build()
        // create the session object
        final sess = new Session(config)
        try {
            if( cmd == 'cache-backup')
                cacheBackup()
            if( cmd == 'cache-restore' )
                cacheRestore()
            return 0
        }
        catch (Throwable e) {
            log.error("Unable to perform plugin command: $cmd - Cause: ${e.message}", e)
            return 1
        }
        finally {
            sess.destroy()
        }
    }

    protected void cacheBackup() {
        log.debug "Running Nextflow cache backup"
        new CacheManager(System.getenv()).saveCacheFiles()
    }

    protected void cacheRestore() {
        log.debug "Running Nextflow cache restore"
        new CacheManager(System.getenv()).restoreCacheFiles()
    }

}
