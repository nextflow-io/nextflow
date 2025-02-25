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

package nextflow.secret

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.errors.NotSupportedException

/**
 * A secret provider only used to report an error when secrets are disabled
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@PackageScope
@CompileStatic
class NullProvider implements SecretsProvider {
    @Override
    boolean activable() {
        return false
    }

    @Override
    SecretsProvider load() {
        return this
    }

    @Override
    Secret getSecret(String name) {
        throw new AbortOperationException("Unable to access 'secrets.$name' because secrets feature is disabled - Enable it setting the variable NXF_ENABLE_SECRETS=true in your environment")
    }

    @Override
    void putSecret(String name, String value) {
        throw new NotSupportedException("Operation 'putSecret' is not supported by ${this.class.name}")
    }

    @Override
    void removeSecret(String name) {
        throw new NotSupportedException("Operation 'removeSecret' is not supported by ${this.class.name}")
    }

    @Override
    Set<String> listSecretsNames() {
        throw new NotSupportedException("Operation 'listSecretsNames' is not supported by ${this.class.name}")
    }

    @Override
    String getSecretsEnv(List<String> secretNames) {
        throw new NotSupportedException("Operation 'getSecretsEnv' is not supported by ${this.class.name}")
    }

    @Override
    String getSecretsEnv() {
        throw new NotSupportedException("Operation 'getSecretsEnv' is not supported by ${this.class.name}")
    }

    @Override
    void close() throws IOException {

    }
}
