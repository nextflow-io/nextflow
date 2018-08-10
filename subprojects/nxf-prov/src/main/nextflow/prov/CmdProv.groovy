/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */
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
        roGenerator.generateWorkflowFolder(bundle)
        /**
         * Generate the Data folder with th input files
         */
        roGenerator.generateDataFolder(bundle)
        /**
         * Generate the Output folder
         */
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

        fileContents = fileContents.replace("\n",',') //put everything on a single line
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
