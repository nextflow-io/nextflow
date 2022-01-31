/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.k8s.model

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PodMountConfigTest extends Specification {

    def 'should create mount for configmap' () {

        when:
        def opt = new PodMountConfig(mountPath: '/etc/some/name', config: 'here' )
        then:
        opt.mountPath == '/etc/some/name'
        opt.fileName == null
        opt.configName == 'here'
        opt.configKey == null

        when:
        opt = new PodMountConfig(mountPath: '/etc/some/name', config: 'here/there.txt' )
        then:
        opt.mountPath == '/etc/some'
        opt.fileName == 'name'
        opt.configName == 'here'
        opt.configKey == 'there.txt'

    }

}
