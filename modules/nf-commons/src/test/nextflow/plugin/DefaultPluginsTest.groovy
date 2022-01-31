/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.plugin

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultPluginsTest extends Specification {

    def 'should validate get plugins' () {
        given:
        def AMAZON = new PluginSpec('nf-amazon', '0.1.0')
        def GOOGLE = new PluginSpec('nf-google', '0.2.0')
        def AZURE = new PluginSpec('nf-azure', '0.3.0')
        and:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': AMAZON,
                'nf-google': GOOGLE,
                'nf-azure': AZURE  ])

        expect:
        defaults.getPlugin('nf-amazon') == AMAZON
        defaults.hasPlugin('nf-amazon')
        and:
        !defaults.hasPlugin('nf-foo')

        when:
        defaults.getPlugin('nf-foo')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Unknown Nextflow plugin 'nf-foo'"
    }

    def 'should get sorted string' () {
        given:
        def DELTA = new PluginSpec('nf-delta', '1.0.0')
        def OMEGA = new PluginSpec('nf-omega', '2.0.0')
        def ALPHA = new PluginSpec('nf-alpha', '3.0.0')
        and:
        def defaults = new DefaultPlugins(plugins: [
                'nf-delta': DELTA,
                'nf-omega': OMEGA,
                'nf-alpha': ALPHA  ])

        expect:
        defaults.toSortedString() == 'nf-alpha@3.0.0,nf-delta@1.0.0,nf-omega@2.0.0'
    }

}
