/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.config


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

}
