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

package nextflow.module

import nextflow.module.spi.RemoteModuleResolverProvider
import spock.lang.Specification

/**
 * Test for DefaultRemoteModuleResolver SPI implementation
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class DefaultRemoteModuleResolverTest extends Specification {

    def 'should load resolver via SPI'() {
        when:
        def resolver = RemoteModuleResolverProvider.getInstance()

        then:
        resolver != null
        resolver.class.name == 'nextflow.module.DefaultRemoteModuleResolver'
        resolver.priority == 0
    }

    def 'should return default priority'() {
        given:
        def resolver = new DefaultRemoteModuleResolver()

        expect:
        resolver.getPriority() == 0
    }
}
