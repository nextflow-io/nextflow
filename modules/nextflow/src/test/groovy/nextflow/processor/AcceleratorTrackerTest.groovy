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

package nextflow.processor

import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class AcceleratorTrackerTest extends Specification {

    @Unroll
    def 'should create tracker from device environment variables'() {
        given:
        def tracker = AcceleratorTracker.create(ENV)
        
        expect:
        tracker.name() == NAME
        tracker.total() == TOTAL
        tracker.available() == TOTAL
        
        where:
        ENV                                                         | NAME                     | TOTAL
        [:]                                                         | null                     | 0
        ['CUDA_VISIBLE_DEVICES': '0,1,2']                           | 'CUDA_VISIBLE_DEVICES'   | 3
        ['HIP_VISIBLE_DEVICES': 'gpu0,gpu1']                        | 'HIP_VISIBLE_DEVICES'    | 2
        ['ROCR_VISIBLE_DEVICES': 'device1,device2,device3,device4'] | 'ROCR_VISIBLE_DEVICES'   | 4
    }

    def 'should handle single device ID'() {
        given:
        def env = ['CUDA_VISIBLE_DEVICES': '0']
        
        when:
        def tracker = AcceleratorTracker.create(env)
        
        then:
        tracker.name() == 'CUDA_VISIBLE_DEVICES'
        tracker.total() == 1
        tracker.available() == 1
    }

    def 'should handle UUID device IDs'() {
        given:
        def env = ['CUDA_VISIBLE_DEVICES': 'GPU-12345678-1234-1234-1234-123456789abc,GPU-87654321-4321-4321-4321-cba987654321']
        
        when:
        def tracker = AcceleratorTracker.create(env)
        
        then:
        tracker.name() == 'CUDA_VISIBLE_DEVICES'
        tracker.total() == 2
        tracker.available() == 2
    }

    def 'should acquire and release permits correctly'() {
        given:
        def env = ['CUDA_VISIBLE_DEVICES': '0,1,2']
        def tracker = AcceleratorTracker.create(env)
        
        when:
        def acquired1 = tracker.acquire(2)
        
        then:
        acquired1.size() == 2
        acquired1.containsAll(['0', '1'])
        tracker.available() == 1
        
        when:
        def acquired2 = tracker.acquire(1)
        
        then:
        acquired2.size() == 1
        acquired2.contains('2')
        tracker.available() == 0
        
        when:
        tracker.release(acquired1)
        
        then:
        tracker.available() == 2
        
        when:
        tracker.release(acquired2)
        
        then:
        tracker.available() == 3
    }

    def 'should handle acquiring all permits'() {
        given:
        def env = ['CUDA_VISIBLE_DEVICES': '0,1,2']
        def tracker = AcceleratorTracker.create(env)
        
        when:
        def acquired = tracker.acquire(3)
        
        then:
        acquired.size() == 3
        acquired.containsAll(['0', '1', '2'])
        tracker.available() == 0
        
        when:
        tracker.release(acquired)
        
        then:
        tracker.available() == 3
    }

    def 'should handle empty device list'() {
        given:
        def env = ['CUDA_VISIBLE_DEVICES': '']
        
        when:
        def tracker = AcceleratorTracker.create(env)
        
        then:
        tracker.name() == 'CUDA_VISIBLE_DEVICES'
        tracker.total() == 0
        tracker.available() == 0
    }

}
