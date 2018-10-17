package nextflow.config

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ManifestTest extends Specification {

    def 'should check manifest object' () {

        given:
        def MAN = [author: 'pablo', nextflowVersion: '1.2.3', name: 'foo']

        when:
        def manifest = new Manifest(MAN)
        then:
        manifest.with {
            author == 'pablo'
            nextflowVersion == '1.2.3'
            name == 'foo'
        }

    }

    def 'should check empty manifest' () {

        // check empty manifest
        when:
        def manifest = new Manifest(new ConfigObject())
        then:
        manifest.with {
            homePage == null
            defaultBranch == 'master'
            description == null
            author == null
            mainScript == 'main.nf'
            gitmodules == null
            nextflowVersion == null
            version == null
            name == null
        }

    }


}
