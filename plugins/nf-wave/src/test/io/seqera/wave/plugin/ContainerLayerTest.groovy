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

package io.seqera.wave.plugin

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerLayerTest extends Specification {

    def 'should convert to string' () {
        when:
        def l1 = new ContainerLayer( 'data:ABC1234567890', 'sha256:12345', 100, 'sha256:67890' )
        then:
        l1.toString() == 'ContainerLayer[location=data:ABC1234567890; tarDigest=sha256:67890; gzipDigest=sha256:12345; gzipSize=100]'

        when:
        def l2 = new ContainerLayer( 'data:12345678901234567890', 'sha256:12345', 100, 'sha256:67890' )
        then:
        l2.toString() == 'ContainerLayer[location=data:12345678901234567890; tarDigest=sha256:67890; gzipDigest=sha256:12345; gzipSize=100]'

        when:
        def l3= new ContainerLayer( 'data:12345678901234567890x', 'sha256:12345', 100, 'sha256:67890' )
        then:
        l3.toString() == 'ContainerLayer[location=data:12345678901234567890...; tarDigest=sha256:67890; gzipDigest=sha256:12345; gzipSize=100]'

    }

    
}
