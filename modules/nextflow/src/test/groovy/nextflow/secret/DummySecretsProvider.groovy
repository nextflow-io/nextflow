/*
 * Copyright 2021, Sage-Bionetworks
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

package nextflow.secret

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DummySecretsProvider implements SecretsProvider {

    private Map<String,Secret> target = [:]

    DummySecretsProvider(Map<String,String> secrets) {
        secrets.each {
            target.put(it.key, new SecretImpl(it.key,it.value))
        }
    }

    @Override
    boolean activable() { true }

    @Override
    SecretsProvider load() {
        return this
    }

    @Override
    Secret getSecret(String name) {
        return target.get(name)
    }

    @Override
    void putSecret(String name, String value) {
        target.put(name, new SecretImpl(name,value))
    }

    @Override
    void removeSecret(String name) {
       target.remove(name)
    }

    @Override
    Set<String> listSecretsNames() {
        return target.keySet()
    }

    @Override
    void close() throws IOException { }

    @Deprecated
    String getSecretsEnv() {
        return null
    }

    @Override
    String getSecretsEnv(List<String> names) {
        String result = ''
        target.each { k,v -> result += "export $k=$v\n" }
        return result
    }
}
