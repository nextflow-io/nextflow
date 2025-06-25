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

package nextflow.util

import java.util.concurrent.Semaphore

import groovy.transform.CompileStatic

/**
 * Specialized semaphore that keeps track of which slots
 * are being used.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class TrackingSemaphore {
    private final Semaphore semaphore
    private final Map<Integer,Boolean> availIds

    TrackingSemaphore(int permits) {
        semaphore = new Semaphore(permits)
        availIds = new HashMap<>(permits)
        if ( SysEnv.containsKey('CUDA_VISIBLE_DEVICES') ) {
            final cudaDevices = SysEnv.get('CUDA_VISIBLE_DEVICES').split(',').collect { it.trim() as int }
            if ( cudaDevices.size() < permits) {
                throw new IllegalConfigException("The environment variable CUDA_VISIBLE_DEVICES contains ${cudaDevices.size()} device(s), but ${permits} gpus were requested. Please update either `CUDA_VISIBLE_DEVICES` or `excutor.gpus`.")
                }
            else {
                for ( int i=0; i< cudaDevices.size(); i++ ) {
                availIds.put(cudaDevices[i], true)
                }
            }
        }
        else {
        for( int i=0; i<permits; i++ )
            availIds.put(i, true)
        }
    }

    int availablePermits() {
        return semaphore.availablePermits()
    }

    List<Integer> acquire(int permits) {
        semaphore.acquire(permits)
        final result = new ArrayList<Integer>(permits)
        for( final entry : availIds.entrySet() ) {
            if( entry.getValue() ) {
                entry.setValue(false)
                result.add(entry.getKey())
            }
            if( result.size() == permits )
                break
        }
        return result
    }

    void release(List<Integer> ids) {
        semaphore.release(ids.size())
        for( id in ids )
            availIds.put(id, true)
    }

}
