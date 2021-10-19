/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.config

import nextflow.secret.DummySecretsProvider
import nextflow.secret.SecretHolder
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigMapTest extends Specification {

    def 'should allow cfg attributes' () {

        given:
        def cfg = new ConfigMap([foo:1, bar:'two'])
        expect:
        cfg.foo == 1
        cfg.bar == 'two'
        cfg.unknown == null
        and:
        cfg.get('foo') == 1
        cfg.get('bar') == 'two'
        cfg.get('unknown') == null

    }

    def 'should return secret values' () {
        given:
        def SECRET1 = "I'm the secret value"
        def SECRET2 = "The other secret"
        and:
        def provider = new DummySecretsProvider([foo:SECRET1, bar:SECRET2])
        def holder1 = new SecretHolder('foo')
        def holder2 = new SecretHolder('bar')
        def cfg = new ConfigMap()

        when:
        cfg.put('alpha', holder1)
        cfg.put('delta', new ConfigMap(gamma: holder2))
        then:
        cfg.alpha == holder1
        cfg.delta.gamma == holder2 

        when:
        cfg.withSecretProvider(provider)
        then:
        cfg.alpha == SECRET1
        cfg.delta.gamma == SECRET2
    }

}
