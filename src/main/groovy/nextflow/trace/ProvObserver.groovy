package nextflow.trace

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j;
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.prov.ProvenanceGenerator
import nextflow.prov.*


import org.openprovenance.prov.model.*
import org.apache.taverna.robundle.Bundle;


/**
 * Created by edgar on 10/05/18.
 */

@Slf4j
@CompileStatic
public class ProvObserver implements TraceObserver {

    private ProvenanceGenerator provGenerator = new ProvenanceGenerator()
    private ResearchObjectGenerator roGenerator = new ResearchObjectGenerator()

    private Map configManifest
    private Map configMap

    private String author

    @Override
    public void onFlowStart(Session session) {
        /**
         * Load Manifest from config file to a local variable
         */
        configManifest= getManifestConfig(session)
        configMap = getConfigMap(session)
        author = getManifestConfig(session).author
    }

    @Override
    public void onFlowComplete() {

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
        roGenerator.setManifest(bundle, configManifest)

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
         * save data into METADATA file
         */
        roGenerator.generateMetadataFile(bundle,configMap)

        /**
         * save data into LOG file
         */
        roGenerator.generateLogFile(bundle,author)

        /**
        * save the entire ZIP bundle
        */
        roGenerator.saveBundle(bundle)
    }

    @Override
    public void onProcessCreate(TaskProcessor process) {

    }

    @Override
    public void onProcessSubmit(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    public void onProcessStart(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    public void onProcessComplete(TaskHandler handler, TraceRecord trace) {

        provGenerator.generateProvenance(trace)
    }

    @Override
    public void onProcessCached(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    public boolean enableMetrics() {
        return false;
    }

    private Map getConfigMap(Session session){
        return session.config instanceof Map ? (Map)session.config : Collections.emptyMap()
    }

    private Map getManifestConfig(Session session) {
        return session.config.manifest instanceof Map ? (Map)session.config.manifest : Collections.emptyMap()
    }

}