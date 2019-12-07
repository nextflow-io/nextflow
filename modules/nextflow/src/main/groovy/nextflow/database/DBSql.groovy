package nextflow.database

import groovy.util.logging.Slf4j
import groovy.sql.Sql
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.extension.CH

@Slf4j
class DBSql {
    private DataflowWriteChannel target
    private String dbUrl
    private String dbUser
    private String dbPassword
    private String dbDriver
    private Sql connection
    private boolean delayConnection = false
    private query
    private columnCount = 0
    private List<String> columns = [] 
    private Closure metaclosure = { meta ->
	columns = []
	columnCount = meta.getColumnCount()
        for (int i = 1; i <= columnCount; i++) {
	    columns << meta.getColumnLabel(i)
	}
    }

    protected void init(Map opts) {
        if( opts.dbUrl )
            dbUrl = opts.dbUrl
        if( opts.dbUser != null )
            dbUser = opts.dbUser
        if( opts.dbPassword )
            dbPassword = opts.dbPassword
        if (opts.dbDriver)
            dbDriver = opts.dbDriver
        if (opts.dbDelayConnecton != null)
            delayConnection = opts.dbDelayConnection as boolean
    }

    DBSql() {

    }

    DBSql(DataflowWriteChannel target, Map opts) {
        this.target = target
        init(opts)
        if (!delayConnection)
            connect()
    }

    DBSql setQuery(query) {
        this.query = query
        return this
    }

    void connect() {
        this.connection=Sql.newInstance(dbUrl, dbUser, dbPassword, dbDriver)
    }


    DataflowWriteChannel apply(Closure closure) {
        if( target == null )
            target = new DataflowQueue()
	    
        if( connection == null )
            connect()
	    
        query0(query, closure)
        
        target.bind(Channel.STOP)
        connection.close()
        return target
    }

    protected void query0(Object query, Closure closure) {
        connection.eachRow(query, metaclosure) { row ->
	    Map<String, Object> lhm = new LinkedHashMap<String, Object>(columnCount, 1);
	    for (int i = 1; i <= columnCount; i++) {
		lhm.put(columns.get(i-1), row.getObject(i));
	    }
            target.bind(closure.call(lhm))
        }
    }
}
