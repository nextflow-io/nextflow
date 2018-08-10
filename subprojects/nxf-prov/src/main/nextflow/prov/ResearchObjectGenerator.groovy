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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.trace.TraceRecord
import org.apache.taverna.robundle.Bundle
import org.apache.taverna.robundle.Bundles
import org.apache.taverna.robundle.manifest.Agent
import org.apache.taverna.robundle.manifest.Manifest
import org.apache.taverna.robundle.manifest.PathMetadata

import java.nio.file.*
import java.text.SimpleDateFormat
import groovy.io.FileType

/**
 *
 * @author Edgar Garriga <edgano@gmail.com>
 */
@Slf4j
@CompileStatic
public class ResearchObjectGenerator {
    /**
     * Filename and file structure of the RO zip (bundle) file
     */
    //RO structure names
    private String roBundleName="XXXXbundle"
    private String metadataFolder = "metadata"
        private String metadataFileName = "metadata.xml"
        private String provenanceFileName = "provenance.json"
    private String logFileName = "log.txt"
    private String snapshotFolderName = "snapshot"
    private String workflowFolderName = "workflow"
    private String dataFolderName = "data"
    private String outputFolderName= "outputs"
    // .ro
        //bundle manifest.xml
    private String commandLineFileName = "commandLine.txt"
    private String enviromentFileName = "enviroment.txt"


    private String author
    private String authorORCID
    private String outDirFolder
    private String baseDir
    private boolean dockerImage
    private boolean singularityImage
    private String containerTech
    private String containerName
    private String containerSha
    private String dockerSHAPrefix ="sha256:"
    private String singularitySHAPrefix=""
    private String commandLine
    private String uuid
    private String nextflowVersionPrefix = "nextflow version "
    private String nextflowVersion
    private String orcidERROR= "**ORCID_not_provided**"

    private List<String> inputFiles =[]
    private List<String> outputFiles =[]

    ResearchObjectGenerator(){}

    /**
     * Creation of the RO zip bundle object
     * @return Bundle
     */
    public Bundle generateROBundle(){
        return Bundles.createBundle()
    }

    /**
     * Generation of all the file system and structure of the RO bundle
     * @param bundle
     */
    public void generateFileStructure(Bundle bundle){
        Path metaFolderPath = bundle.getRoot().resolve(metadataFolder);
        createROdirectory(metaFolderPath);

        Path snapFolderPath = bundle.getRoot().resolve(snapshotFolderName);
        createROdirectory(snapFolderPath);

        Path workflowFolderPath = bundle.getRoot().resolve(workflowFolderName);
        createROdirectory(workflowFolderPath);

        Path dataFolderPath = bundle.getRoot().resolve(dataFolderName);
        createROdirectory(dataFolderPath);

        Path outputFolderPath = bundle.getRoot().resolve(outputFolderName);
        createROdirectory(outputFolderPath);
        log.debug "RO structure generated"
    }

    /**
     * Set :
     *      -Author
     *      -ORCID
     *      -ID
     *      -Agregation Files
     * into the RO manifest.
     * @param bundle
     */
    public void setManifest(Bundle bundle){
        // https://github.com/apache/incubator-taverna-language/blob/master/taverna-robundle/src/test/java/org/apache/taverna/robundle/manifest/TestManifestJSON.java
        Manifest manifest = bundle.getManifest();

        setAuthorInformation(manifest)

        manifest.setId(URI.create("/"))

        setAggregationManifest(manifest)

        log.debug "set RO manifest"
    }

    /**
     * Create/update the log file
     * @param bundle
     */
    public void generateLogFile(Bundle bundle){
        File logFile = setLogInfo()
        fileToBundle(bundle,Paths.get(logFile.path),logFileName, metadataFolder)
        log.debug "Generate Log File DONE"
    }

    /**
     * Add the provenance file to the RO bundle
     * @param bundle
     */
    public void addProvenanceFile(Bundle bundle){
        File file = new File(provenanceFileName)
        fileToBundle(bundle, Paths.get(file.getPath()), provenanceFileName,metadataFolder)
        boolean result = removeFile(provenanceFileName)
        log.debug "Add prov to bundle: ${result} DONE"
    }

    /**
     * Insert all the input files into the data folder in the RO bundle
     * @param bundle
     */
    public void generateDataFolder(Bundle bundle){
        inputFiles.unique()
        for (file in inputFiles){
            Path auxPath = Paths.get(file)
            fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), dataFolderName)
        }
        log.debug "Generate Data folder DONE"
    }

    /**
     * Insert all the output files into the output folder in the RO bundle
     * @param bundle
     *
     * It needs to be declared into the config file
     */
    public Path generateOutputFolder(Bundle bundle){
        outputFiles.unique()
        for (file in outputFiles){
            Path auxPath = Paths.get(file)
            fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), outputFolderName)
        }
        log.debug "Generate Output folder DONE"
    }

    /**
     * Copy the pipeline's work directory into the RO bundle
     * @param bundle
     */
    public void generateWorkflowFolder(Bundle bundle){
        def result = getFilesFromDir(baseDir)
        log.debug "Get Files from workflow dir DONE"
        for (element in result) {
            Path auxPath = Paths.get(element.toString())
            if (auxPath.toString().contains(".git")||auxPath.toString().contains(".DS_Store")) {}   //TODO FIX encoded problem with git files
            else {
                String trimmedPath = auxPath.toString().replace(baseDir, ""); //remove the baseDir from the filePath -> relativePath for bundle
                trimmedPath = trimmedPath.replace(auxPath.getFileName().toString(), "");// remove fileName from the relativePath

                if (Paths.get(trimmedPath).getNameCount() == 0) {   //if the file is in the "relative" root -> file2bundle
                    fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), "${workflowFolderName}${trimmedPath}")
                    log.debug "File ${auxPath.getFileName().toString()} workdir into bundle DONE"
                } else {    //the file is not in the root -> there are intermediate folder
                    directoryToBundle(bundle, Paths.get(trimmedPath),workflowFolderName)
                    fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), "${workflowFolderName}${trimmedPath}")
                }
            }
        }
        log.debug "Generate Workflow folder DONE"
    }

    /**
     * Fill the Snapshot folder with the commandLine file
     * @param bundle
     */
    public void generateSnapshot(Bundle bundle){
        File scriptFile = generateScript()
        fileToBundle(bundle,Paths.get(scriptFile.path),commandLineFileName, snapshotFolderName)
        boolean result = removeFile(commandLineFileName)
        log.debug "Generate Snapshot file: ${result}"
    }

    /**
     * Fill the Metadata folder with the files needed
     * @param bundle
     */
    public void generateMetadataFolder(Bundle bundle){
        Path metadataFile = generateMetadataFile()
        fileToBundle(bundle, metadataFile,metadataFileName,metadataFolder)
        log.debug "Generate Metadata folder DONE"
    }

    /**
     * Get the files path of the given directory
     * @param directory
     * @return A List with the Path of all the files in the directory
     */
    protected List getFilesFromDir(String directory){
        //https://stackoverflow.com/questions/3953965/get-a-list-of-all-the-files-in-a-directory-recursive
        def list = []
        def dire = new File(directory)
        dire.eachFileRecurse (FileType.FILES) { file ->
            list << file
        }
        return list
    }

    /**
     * Create the metadata file with teh required information
     * @return Path of the metadata file
     */
    protected Path generateMetadataFile(){
        Path metadataFile = Files.createTempFile("","");
        metadataFile.write"Command Line: ${commandLine}\n"
        metadataFile << "UUID: ${uuid}\n"
        metadataFile << "Nextflow version: ${nextflowVersion}\n"

        return metadataFile
    }

    /**
     * Fill inputFiles variable with the trace.input file's path
     * @param trace
     */
    public void getCleanInputFiles(TraceRecord trace){
        def inputFiles = trace.getFmtStr("input")
        inputFiles= inputFiles.replaceAll("\\s","")         //remove whitespace from the string
        List<String> inputList = Arrays.asList(inputFiles.split(';'));

        for (inputElem in inputList){
            if (!inputElem.contains(/work/)){
                this.inputFiles.add(inputElem.trim())
            }
        }
        log.debug "Get Clean Input Files (data folder)"
    }

    /**
     * Fill outputFiles variable with the file's path from the output folder
     * It is needed to have declared the output dir into the config file
     */
    public void getCleanOutputFiles(){
        if (outDirFolder!=null){ // because we miss the "concept" null when we add it on the map>>file>>map
            def outDirFiles = getOutDirFiles(outDirFolder)
            File fileAux = new File ("${outDirFolder}/${outDirFiles[0].toString()}")
            String filePath = fileAux.absolutePath.substring(0,fileAux.absolutePath.lastIndexOf(File.separator));
            for (element in outDirFiles){
                String aux = "${filePath}/${element}"
                outputFiles.add(aux)
            }
        }else{
            // NOT NEEDED, its controled on the CMDPROV --> log.warn("You need to specify the output directory inside nextflow.config/params for the RO zip file")
        }
        log.debug "Get Clean Output Files (output folder)"
    }

    /**
     * Function to insert a file into the RO bundle
     * @param bundle
     * @param filePath
     * @param fileName
     * @param folderName
     */
    protected void fileToBundle(Bundle bundle, Path filePath, String fileName, String folderName){
        if (filePath.isFile()){
            Path folderPath = bundle.getRoot().resolve(folderName);
            Path bundleFilePath = folderPath.resolve(fileName);
            Bundles.setStringValue(bundleFilePath, filePath.text.toString());
            Files.copy(bundleFilePath, filePath, StandardCopyOption.REPLACE_EXISTING);
        }else if (filePath.isDirectory()){
            log.warn("The element: \"${fileName}\" is a directory")
            //directoryToBundle(bundle,filePath,fileName,folderName)
        }
        log.debug "File to bundle: ${fileName} to ${folderName}"
    }

    /**
     * Function to create a directory/subdirectory into the RO bundle
     * @param bundle
     * @param filePath
     * @param folderName
     */
    private void directoryToBundle(Bundle bundle, Path filePath, String folderName){
        //TODO --> commented in fileToBundle -> not needed to control there??
        Path auxPath = bundle.getRoot().resolve("${folderName}");
        //**
        // generate intermeditate directories
        //**
        for (int i=0; i<filePath.getNameCount(); i++){
            auxPath = Paths.get("${auxPath}/${filePath.getName(i)}")
            String stringAuxPath = auxPath.toString().substring(1) // remove first "/"
            Path intermediateFolderPath = bundle.getRoot().resolve(stringAuxPath);

            if (!Files.exists(intermediateFolderPath)) {
                Files.createDirectory(intermediateFolderPath);
                log.debug "Directory ${intermediateFolderPath.toString()} Created!"

            }
        }
    }

    /**
     * Save the bundle into a Zip file
     * @param bundle
     */
    public void saveBundle(Bundle bundle) {
        Path ro = Paths.get("${System.getProperty("user.dir")}/${roBundleName}.zip")

        //Remove the older zip if it exist
        boolean pathExists =Files.exists(ro, LinkOption.NOFOLLOW_LINKS);
        if (pathExists) {
            try {
                Files.delete(ro);
            } catch (IOException e) {
                //deleting file failed
                e.printStackTrace();
            }
        }
        Path zip = Files.createFile(ro)
        Bundles.closeAndSaveBundle(bundle, zip);
        log.info "RO Bundle saved to ${zip}"
    }

    /**
     * Get the files from the output folder if it exists
     * @param outDirFolder
     * @return List with the output file's path
     */
    private List getOutDirFiles(String outDirFolder){
        if (outDirFolder!=null){
            log.debug "OutDir folder exist"
            return getFilesFromDir(outDirFolder)
        }
        log.warn "OutDir folder NOT exist"
    }

    /**
     * Set Author's name, author's ORCID into the RO manifest
     * @param manifest
     */
    protected void setAuthorInformation(Manifest manifest) {
        if (author!=null){
            Agent createdBy = new Agent(author);

            if (authorORCID!=null) {
                createdBy.setOrcid(URI.create(authorORCID));
            }else{      //>> Do we want to save an "ORCID ERROR" or just not create the instance?
                createdBy.setOrcid(URI.create(orcidERROR));
            }
            manifest.setCreatedBy(createdBy);
        }else{
            // what to do if we dont have author on the manifest
            // We dont "generate" the entity --> not ORCID allow neither
        }
        log.debug "set Author info to RO manifest"
    }

    /**
     * Auxiliar function to recover each value from the map if it's exist
     * @param map
     * @param field
     * @return String with the value of the map
     */
    private String getFieldMap(Map map, String field){
        if(map.get(field) != "not_defined"){ //from map to string to file -> to map again... we miss the "concept" null, and it stay as a string value
            return map.get(field)
        }else{
            log.warn("The field ${field} is not configured on the project.")
        }
    }

    /**
     * Set all the agrgegates files into the RO manifest
     * @param manifest
     */
    protected void setAggregationManifest(Manifest manifest){
        LinkedList<PathMetadata> aggregationList = new LinkedList<PathMetadata>()
        manifest.setAggregates(aggregationList)

        log.debug "Set Aggregation to RO manifest"
    }

    /**
     * Generate the scriptFile with the commandLine for the Snapshot folder
     * @return File with the commandLine used
     */
    protected File generateScript(){
        File scriptFile= new File(commandLineFileName)
        scriptFile.write(commandLine)

        return scriptFile
    }

    /**
     * Remove a file given the file's path
     * @param fileName
     * @return boolean if it was possible to remove the file or not
     */
    private boolean removeFile(String fileName){
        return Files.deleteIfExists(Paths.get(fileName))
    }

    /**
     * Update/generate the log file
     * @return Log File
     */
    protected File setLogInfo(){
        def today = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date())
        File logFile = new File(logFileName)
        logFile.append("${today} -- ${author} -- ${uuid} -- ${commandLine}\n")

        return logFile
    }

    /**
     * Update the prov class variables using a given Map.
     * @param provMap
     */
    public void setProvValues(Map provMap){
        author = getFieldMap(provMap, 'author')
        authorORCID = getFieldMap(provMap,'orcid')
        commandLine= getFieldMap(provMap, 'commandLine')
        uuid= getFieldMap(provMap, 'uuid')
        nextflowVersion= getFieldMap(provMap, 'nfVersion').substring(nextflowVersionPrefix.size())
        outDirFolder = getFieldMap(provMap, 'outDir')
        baseDir = getFieldMap(provMap, 'baseDir')
    }

    /**
     * Print a given map
     * @param map
     */
    private void printMap(Map map){
        for (Map.Entry<String, String> element : map.entrySet()){
            print   "-->>value: ..${element.getKey()}..  -- ${element.getValue()} \n"
        }
    }

    /**
     * Create the folders of the RO structure
     * @param path
     */
    protected void createROdirectory(Path path){
        if (!Files.exists(path)) {
            try {
                Files.createDirectories(path);
            } catch (IOException e) {
                //fail to create directory
                e.printStackTrace();
            }
        }
    }
}