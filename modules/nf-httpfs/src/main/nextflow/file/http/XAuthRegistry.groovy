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

package nextflow.file.http

import groovy.transform.CompileStatic

/**
 * Provides a pluggable authentication registry for {@link XFileSystemProvider}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Singleton(strict = false)
@CompileStatic
class XAuthRegistry {

    private List<XAuthProvider> providers = new ArrayList<>()

    protected XAuthRegistry() {}

    void register(XAuthProvider provider) {
        providers.add(provider)
    }

    void unregister(XAuthProvider provider) {
        final p = providers.indexOf(provider)
        if( p!=-1 )
            providers.remove(p)
    }

    boolean authorize(URLConnection connection) {
        for( XAuthProvider it : providers ) {
            if( it.authorize(connection) )
                return true
        }
        return false
    }

    boolean refreshToken(URLConnection connection) {
        for( XAuthProvider it : providers ) {
            if( it.refreshToken(connection) )
                return true
        }
        return false
    }
}
