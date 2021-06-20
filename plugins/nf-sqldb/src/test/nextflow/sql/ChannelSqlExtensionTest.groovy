package nextflow.sql

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelSqlExtensionTest extends Specification {

    def 'should normalise query' () {
        given:
        def ext = new ChannelSqlExtension()

        expect:
        ext.normalize('select * from x; ') == 'select * from x;'
        ext.normalize('select * from x ')  == 'select * from x;'
    }

    def 'should create db params' () {
        given:
        def ext = new ChannelSqlExtension()
        when:
        def props = ext.dbProps([:])
        then:
        props.driver == 'org.h2.Driver'
        props.url == 'jdbc:h2:mem:'

    }

    def 'should connect db' () {
        given:
        def ext = new ChannelSqlExtension()
        when:
        def conn = ext.connect([:])
        then:
        conn != null
        cleanup:
        conn.close()
    }
}
