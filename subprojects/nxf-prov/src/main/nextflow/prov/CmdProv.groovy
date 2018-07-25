package nextflow.prov

import nextflow.CacheDB
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.trace.TraceRecord
import nextflow.util.HistoryFile
import org.apache.taverna.robundle.Bundle
import org.openprovenance.prov.model.Document

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdProv extends CmdBase {

    private Map provMap=[:]

    ProvenanceGenerator provGenerator = new ProvenanceGenerator()
    ResearchObjectGenerator roGenerator = new ResearchObjectGenerator()


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

        provMap = this.readProvMetadataFile("${db.dataDir}/prov.${db.runName}")
        roGenerator.setProvValues(provMap)

        Document provDocument = provGenerator.getProvDocument()
        /**
         * Introduce all the provenance object to the PROV file
         */
        provGenerator.setElementsToProvFile(provDocument)

        /**
         * Generate the RO BUNDLE
         */
        Bundle bundle = roGenerator.generateROBundle()

        roGenerator.generateFileStructure(bundle)

        /**
         * Modify the bundle manifest
         */
        roGenerator.setManifest(bundle)

        /**
         * Generate the PROV provDocument
         */
        provGenerator.generateProvFile(provDocument)
        roGenerator.addProvenanceFile(bundle)

        /**
         * Generate the snapshot content
         */
        roGenerator.generateSnapshot(bundle)
        /**
         * Generate Workflow folder
         */
        //TODO check random ERROR
        //  Exception in thread "Thread-1" groovy.lang.GroovyRuntimeException: exception while reading process stream
        //  at org.codehaus.groovy.runtime.ProcessGroovyMethods$TextDumper.run(ProcessGroovyMethods.java:496)
        //
        // at java.lang.Thread.run(Thread.java:748)
        // Caused by: java.io.IOException: Stream closed
        // 	    at java.io.BufferedInputStream.getBufIfOpen(BufferedInputStream.java:170)
        //  	at java.io.BufferedInputStream.read(BufferedInputStream.java:336)
        //	    at sun.nio.cs.StreamDecoder.readBytes(StreamDecoder.java:284)
        //	    at sun.nio.cs.StreamDecoder.implRead(StreamDecoder.java:326)
        //	    at sun.nio.cs.StreamDecoder.read(StreamDecoder.java:178)
        //	    at java.io.InputStreamReader.read(InputStreamReader.java:184)
        //	    at java.io.BufferedReader.fill(BufferedReader.java:161)
        //	    at java.io.BufferedReader.readLine(BufferedReader.java:324)
        //	    at java.io.BufferedReader.readLine(BufferedReader.java:389)
        //	    at org.codehaus.groovy.runtime.ProcessGroovyMethods$TextDumper.run(ProcessGroovyMethods.java:489)
        //	    ... 1 more
        roGenerator.generateWorkflowFolder(bundle)
        /**
         * Generate the Data folder with th input files
         */
        roGenerator.generateDataFolder(bundle)
        /**
         * Generate the Output folder
         */
        //TODO check random ERROR
        //  Exception in thread "Thread-1" groovy.lang.GroovyRuntimeException: exception while reading process stream
        //  at org.codehaus.groovy.runtime.ProcessGroovyMethods$TextDumper.run(ProcessGroovyMethods.java:496)
        //
        roGenerator.getCleanOutputFiles()
        roGenerator.generateOutputFolder(bundle)
        /**
         * save data into METADATA file
         */
        roGenerator.generateMetadataFolder(bundle)

        /**
         * save data into LOG file
         */
        roGenerator.generateLogFile(bundle)

        /**
         * save the entire ZIP bundle
         */
        roGenerator.saveBundle(bundle)
    }

    protected void renderProv(TraceRecord record) {
        provGenerator.generateProvenance(record)
        /**
         * Add the input files of the process to the main list
         */
        roGenerator.getCleanInputFiles(record)
    }
    private Map readProvMetadataFile(String path){
        String fileContents=new File(path).text

        fileContents = fileContents.replace("\n",',') //put everithing on a single line
        fileContents = fileContents.replace(" = ",'=') // remove space on the equal symbol

        def newMap = [:]
        fileContents.tokenize(',').each {
            def kvTuple = it.tokenize('=')
            kvTuple[1]=removeComma(kvTuple[1],"'","'")
            newMap[kvTuple[0]] = kvTuple[1]
        }
        return newMap
    }
    private String removeComma(String variable, String ini, String fin){
        if (variable[0].equals(ini) && variable[variable.length()-1].equals(fin)){
            return variable.substring(1, variable.length() - 1);
        } else{
            return variable
        }
    }
}
