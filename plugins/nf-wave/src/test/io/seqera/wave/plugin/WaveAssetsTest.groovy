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

package io.seqera.wave.plugin

import nextflow.script.bundle.ResourcesBundle
import nextflow.util.CacheHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveAssetsTest extends Specification {

    def 'should compute hash key' () {
        given:
        def IMAGE = 'foo:latest'
        def BUNDLE = Mock(ResourcesBundle) { fingerprint() >> '12345' }
        
        expect:
        new WaveAssets(IMAGE).fingerprint() == CacheHelper.hasher([IMAGE]).hash().toString()
//        new WaveAssets(IMAGE,BUNDLE).hashKey() == CacheHelper.hasher([IMAGE, BUNDLE]).hash().toString()

    }

}
