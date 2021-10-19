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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SecretsContextTest extends Specification {

    def 'should resolve secrets' () {
        given:
        def provider = new DummySecretsProvider([foo:'one', bar:'two'])
        def resolver = new SecretsContext()
        and:
        def shell = new GroovyShell(new Binding([secrets:resolver]))
        
        when:
        def secret = shell.evaluate('secrets.foo')
        then:
        secret instanceof SecretHolder

        when:
        def holder = (SecretHolder) secret
        then:
        holder.getSecretName() == 'foo'
        holder.call() == 'secrets.foo'
        
        when:
        holder.getSecretValue()
        then:
        thrown(IllegalStateException)

        when:
        holder.bind(provider)
        then:
        holder.getSecretValue() == 'one'
        holder.call() == 'one'

    }

}
