package nextflow.sql

import java.sql.Connection
import java.sql.ResultSet
import java.sql.Statement
import java.util.concurrent.CompletableFuture

import groovy.sql.Sql
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Global
import nextflow.NF
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.ChannelExtensionPoint
/**
 * Provide a channel factory extension that allows the execution of Sql queries
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelSqlExtension implements ChannelExtensionPoint {

    private static String DEFAULT_URL = 'jdbc:h2:mem:'
    private static String DEFAULT_DRIVER = 'org.h2.Driver'
    private static String DEFAULT_USER = 'sa'

    private Session session
    private Map config

    void init(Session session) {
        this.session = session
        this.config = (Map) session.config.dataSources
    }

    @Override
    String getScope() {
        return 'sql'
    }

    protected Map<String,Object> dbProps(Map<String,Object> opts) {
        Map result = new HashMap()
        result.driver = opts.driver ?: DEFAULT_DRIVER
        result.url = opts.url ?: DEFAULT_URL
        result.user = opts.user ?: DEFAULT_USER
        result.password = opts.password
        return result
    }

    @CompileDynamic
    protected Connection connect(Map<String,Object> opts) {
        def props = dbProps(opts)
        Sql.newInstance(props).getConnection()
    }

    protected String normalize(String q) {
        if( !q )
            throw new IllegalArgumentException("Missing query argument")
        def result = q.trim()
        if( !result.endsWith(';') )
            result += ';'
        return result
    }

    DataflowWriteChannel fromQuery(String query) {
        fromQuery([:], query)
    }

    DataflowWriteChannel fromQuery(Map opts, String query) {
        return queryToChannel(query, opts)
    }

    protected DataflowWriteChannel queryToChannel(String query, Map opts) {
        final channel = CH.create()
        if(NF.dsl2) {
            session.addIgniter {-> final conn = connect(opts); asyncQuery(query, conn, channel) }
        }
        else {
            final conn = connect(opts)
            asyncQuery(query, conn, channel)
        }
        return channel
    }

    protected asyncQuery(String query, Connection conn, DataflowWriteChannel channel) {
        def future = CompletableFuture.runAsync ({ query0(query, conn, channel) })
        future.exceptionally(this.&handlerException)
    }

    static private void handlerException(Throwable e) {
        final error = e.cause ?: e
        log.error(error.message, error)
        final session = Global.session as Session
        session?.abort(error)
    }

    protected void query0(String query, Connection conn, DataflowWriteChannel channel) {
        try {
            try (Statement stm = conn.createStatement()) {
                try( def rs = stm.executeQuery(normalize(query)) ) {
                    emitRows0(rs, channel)
                }
            }
        }
        finally {
            conn.close()
        }
    }

    protected emitRows0(ResultSet rs, DataflowWriteChannel channel) {
        final cols = rs.getMetaData().getColumnCount()
        try {
            while( rs.next() ) {
                def item = new ArrayList(cols)
                for( int i=0; i<cols; i++) {
                    item[i] = rs.getString(i+1)
                }
                // emit the value
                channel.bind(item)
            }
        }
        finally {
            // close the channel
            channel.bind(Channel.STOP)
        }
    }

}
