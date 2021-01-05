package nextflow.sql

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.dag.NodeMarker
import nextflow.Channel

import java.util.concurrent.CompletableFuture

Channel.metaClass.static.propertyMissing = { String name ->
    return new Channel(name) as DataflowWriteChannel
}


Channel.metaClass.static.create << { String title -> new Book(title: title) }

class SQLChannel extends Channel {
    static DataflowWriteChannel fromSql(query) {
        fromSql(Collections.emptyMap(), query)
    }

    static DataflowWriteChannel fromSql(Map opts, query) {
        fromSql(opts, query, null)
    }

    static DataflowWriteChannel fromSql(Map opts, query, Closure closure) {
        //CheckHelper.checkParams('fromSRA', opts, SraExplorer.PARAMS)
        def String queryString = query instanceof GString ? query.toString() : query
        //place something to avoid SQL injection.
        def target = new DataflowQueue()
        def DBSql sql = new DBSql(target, opts).setQuery(queryString)
        if( NF.isDsl2() ) {
            session.addIgniter { fromSql0(sql, closure) }
        } else {
            fromSql0(sql, closure)
        }

        NodeMarker.addSourceNode('Channel.fromSql', target)
        return target
    }

    static private void fromSql0(DBSql sql, Closure closure) {
        def future = CompletableFuture.runAsync({ sql.apply(closure == null ? { row -> row } : closure) } as Runnable)
        fromPath0Future = future.exceptionally(Channel.&handlerException)
    }

}
