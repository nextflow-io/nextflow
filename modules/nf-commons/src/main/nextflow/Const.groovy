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

package nextflow

import java.nio.file.Path
import java.nio.file.Paths
/**
 * Application main constants
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Const {

    static final public String ISO_8601_DATETIME_FORMAT = "yyyy-MM-dd'T'HH:mm:ss'Z'"

    static final public transient BOOL_YES = ['true','yes','on']

    static final public transient BOOL_NO = ['false','no','off']

    /**
     * The application main package name
     */
    static public final String MAIN_PACKAGE = Const.name.split('\\.')[0]

    /**
     * The application main name
     */
    static public final String APP_NAME = MAIN_PACKAGE

    /**
     * The application home folder
     */
    static public final Path APP_HOME_DIR = getHomeDir(APP_NAME)

    static Path sysHome() {
        def home = System.getProperty("user.home")
        if( !home || home=='?' )
            home = System.getenv('HOME')
        if( !home )
            throw new IllegalStateException("Unable to detect system home path - Make sure the variable HOME or NXF_HOME is defined in your environment")
        return Path.of(home)
    }

    private static Path getHomeDir(String appname) {
        final home = System.getenv('NXF_HOME')
        final result = home ? Paths.get(home) : sysHome().resolve(".$appname")

        if( !result.exists() && !result.mkdir() ) {
            throw new IllegalStateException("Cannot create path '${result}' -- check file system access permission")
        }

        return result
    }

    static final Path getAppCacheDir() {
        return Path.of(SysEnv.get('NXF_CACHE_DIR', '.nextflow'))
    }

    static public final String ROLE_WORKER = 'worker'

    static public final String ROLE_MASTER = 'master'

    static public final String MANIFEST_FILE_NAME = 'nextflow.config'

    static public final String DEFAULT_MAIN_FILE_NAME = 'main.nf'

    static public final String DEFAULT_ORGANIZATION = System.getenv('NXF_ORG') ?: 'nextflow-io'

    static public final String DEFAULT_HUB = System.getenv('NXF_HUB') ?: 'github'

    static public final File DEFAULT_ROOT = System.getenv('NXF_ASSETS') ? new File(System.getenv('NXF_ASSETS')) : APP_HOME_DIR.resolve('assets').toFile()

    static public final String DEFAULT_BRANCH = 'master'

    static public final String SCOPE_SEP = ':'

}
