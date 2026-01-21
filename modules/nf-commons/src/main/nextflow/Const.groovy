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
class Const extends AppConst{

    static final public String ISO_8601_DATETIME_FORMAT = "yyyy-MM-dd'T'HH:mm:ss'Z'"

    static final public transient BOOL_YES = ['true','yes','on']

    static final public transient BOOL_NO = ['false','no','off']

    static final Path getAppCacheDir() {
        return Path.of(SysEnv.get('NXF_CACHE_DIR', '.nextflow'))
    }

    static public final String ROLE_WORKER = 'worker'

    static public final String ROLE_MASTER = 'master'

    static public final String SCOPE_SEP = ':'

}
