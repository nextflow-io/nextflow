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

package test

import java.lang.reflect.Method

import groovy.transform.CompileDynamic

/**
 * Test helper to modify the process environment as seen by {@link System#getenv()}.
 *
 * This is needed to test code paths resolving their settings directly from the
 * environment e.g. proxy configuration. It relies on reflection over
 * {@code java.lang.ProcessEnvironment}, therefore it requires the test JVM to be
 * launched with {@code --add-opens=java.base/java.lang=ALL-UNNAMED} and
 * {@code --add-opens=java.base/java.util=ALL-UNNAMED} (already set by the
 * project build).
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@CompileDynamic
class EnvHelper {

    private static final Map ENV_MAP
    private static final Method VAR_OF
    private static final Method VAL_OF

    static {
        final clazz = Class.forName('java.lang.ProcessEnvironment')
        final field = clazz.getDeclaredField('theEnvironment')
        field.setAccessible(true)
        ENV_MAP = (Map) field.get(null)
        VAR_OF = Class.forName('java.lang.ProcessEnvironment$Variable').getDeclaredMethod('valueOf', String)
        VAL_OF = Class.forName('java.lang.ProcessEnvironment$Value').getDeclaredMethod('valueOf', String)
        VAR_OF.setAccessible(true)
        VAL_OF.setAccessible(true)
    }

    /**
     * Run the given action with the specified variables applied to the process
     * environment, restoring the original values on completion. A {@code null}
     * value removes the variable.
     *
     * @param vars The environment variables to apply
     * @param action The action to run
     */
    static void withEnv(Map<String,String> vars, Closure action) {
        final Map<String,String> saved = new HashMap<>()
        for( final entry : vars.entrySet() ) {
            saved.put(entry.key, System.getenv(entry.key))
            setEnv(entry.key, entry.value)
        }
        try {
            action.call()
        }
        finally {
            for( final entry : saved.entrySet() )
                setEnv(entry.key, entry.value)
        }
    }

    /**
     * Set or remove ({@code value == null}) a variable in the process environment
     */
    static void setEnv(String name, String value) {
        if( value != null )
            ENV_MAP.put(VAR_OF.invoke(null, name), VAL_OF.invoke(null, value))
        else
            ENV_MAP.remove(VAR_OF.invoke(null, name))
    }
}
