/*
 * Copyright 2019, WuxiNextcode
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

import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Ólafur Haukur Flyenring <olafurh@wuxinextcode.com>
 */
class PodEmptyDirTest extends Specification {

    def 'should create emptyDirs' () {
        when:
        def opt = new PodEmptyDir('empty-1', '/etc/some/name', PodEmptyDir.PodEmptyDirType.Memory, '2 GB' )
        then:
        opt.name == 'empty-1'
        opt.mountPath == '/etc/some/name'
        opt.type == PodEmptyDir.PodEmptyDirType.Memory
        opt.sizeLimit == new MemoryUnit('2 GB')

        when:
        opt = new PodEmptyDir('empty-2','/etc/some',PodEmptyDir.PodEmptyDirType.Disk)
        then:
        opt.name == 'empty-2'
        opt.mountPath == '/etc/some'
        opt.type == PodEmptyDir.PodEmptyDirType.Disk

        when:
        def entry = [emptyDir: 'empty-3', mountPath: '/etc/some', type: 'Memory']
        opt = new PodEmptyDir(entry)
        then:
        opt.name == 'empty-3'
        opt.mountPath == '/etc/some'
        opt.type == PodEmptyDir.PodEmptyDirType.Memory

        when:
        entry = [emptyDir: 'empty-4', mountPath: '/etc/some', type: 'Disk', sizeLimit: '3 GB']
        opt = new PodEmptyDir(entry)
        then:
        opt.name == 'empty-4'
        opt.mountPath == '/etc/some'
        opt.type == PodEmptyDir.PodEmptyDirType.Disk
        opt.sizeLimit == new MemoryUnit('3 GB')
    }
}