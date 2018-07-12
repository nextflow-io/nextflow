package nextflow.prov

import nextflow.CacheDB
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.trace.TraceRecord
import nextflow.util.HistoryFile
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdProv extends CmdBase {

    @Override
    String getName() {
        return null
    }

    @Override
    void run() {

        def history = HistoryFile.DEFAULT

        if( !history.exists() || history.empty() )
            throw new AbortOperationException("It looks no pipeline was executed in this folder (or execution history is empty)")

        def db = new CacheDB(history.getLast())
        db
            .openForRead()
            .eachRecord(this.&renderProv)
            .close()
    }

    protected void renderProv(TraceRecord record) {

    }
}
