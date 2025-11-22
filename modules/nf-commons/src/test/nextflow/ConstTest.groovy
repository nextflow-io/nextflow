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

package nextflow

import java.nio.file.Path

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConstTest extends Specification {

    def 'should get app default cache dir' () {
        expect:
        Const.appCacheDir == Path.of('.nextflow')
    }

    def 'should get app custom cache dir' () {
        given:
        SysEnv.push(NXF_CACHE_DIR: '/some/path')
        expect:
        Const.appCacheDir == Path.of('/some/path')
        cleanup:
        SysEnv.pop()
    }

}
