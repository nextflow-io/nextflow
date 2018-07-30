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
 * Created by edgar on 19/06/18.
 */
@Slf4j
@CompileStatic
public class ResearchObjectGenerator {
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

    public Bundle generateROBundle(){
        return Bundles.createBundle()
    }

    public void generateFileStructure(Bundle bundle){
        Path metaFolderPath = bundle.getRoot().resolve(metadataFolder);
        Files.createDirectory(metaFolderPath);

        Path snapFolderPath = bundle.getRoot().resolve(snapshotFolderName);
        Files.createDirectory(snapFolderPath);

        Path workflowFolderPath = bundle.getRoot().resolve(workflowFolderName);
        Files.createDirectory(workflowFolderPath);

        Path dataFolderPath = bundle.getRoot().resolve(dataFolderName);
        Files.createDirectory(dataFolderPath);

        Path outputFolderPath = bundle.getRoot().resolve(outputFolderName);
        Files.createDirectory(outputFolderPath);
        log.info "RO structure generated"

    }

    public void setManifest(Bundle bundle){
        // https://github.com/apache/incubator-taverna-language/blob/master/taverna-robundle/src/test/java/org/apache/taverna/robundle/manifest/TestManifestJSON.java
        Manifest manifest = bundle.getManifest();

        this.setAuthorInformation(manifest)

        manifest.setId(URI.create("/"))

        setAggregationManifest(manifest)

        log.info "set RO manifest"
    }

    public void generateLogFile(Bundle bundle){
        File logFile = getLogInfo()
        fileToBundle(bundle,Paths.get(logFile.path),logFileName, metadataFolder)
        log.info "Generate Log File DONE"
    }

    public void addProvenanceFile(Bundle bundle){
        File file = new File(provenanceFileName)
        fileToBundle(bundle, Paths.get(file.getPath()), provenanceFileName,metadataFolder)
        boolean result = removeFile(provenanceFileName)
        log.info "Add prov to bundle: ${result} DONE"
    }

    public void generateDataFolder(Bundle bundle){
        inputFiles.unique()
        for (file in inputFiles){
            Path auxPath = Paths.get(file)
            fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), dataFolderName)
        }
        log.info "Generate Data folder DONE"

    }
    public Path generateOutputFolder(Bundle bundle){
        outputFiles.unique()
        for (file in outputFiles){
            Path auxPath = Paths.get(file)
            fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), outputFolderName)
        }
        log.info "Generate Output folder DONE"

    }

    public void generateWorkflowFolder(Bundle bundle){
        def result = getFilesFromDir(baseDir)
        log.info "Get Files from workflow dir DONE"
        for (element in result){
            Path auxPath = Paths.get("${baseDir}/${element}")
            log.info "Get Path from .${element}. workflow dir path: ${auxPath} DONE"
            fileToBundle(bundle, auxPath,auxPath.getFileName().toString(),workflowFolderName)
            log.info "File ${auxPath.getFileName().toString()} workdir into bundle DONE"
        }
        log.info "Generate Workflow folder DONE"

    }

    public void generateSnapshot(Bundle bundle){
        File scriptFile = generateScript()
        fileToBundle(bundle,Paths.get(scriptFile.path),commandLineFileName, snapshotFolderName)
        boolean result = removeFile(commandLineFileName)
        log.info "Generate Snapshot file: ${result}"

    }

    public void generateMetadataFolder(Bundle bundle){
        Path metadataFile = generateMetadataFile()
        fileToBundle(bundle, metadataFile,metadataFileName,metadataFolder)
        log.info "Generate Metadata folder DONE"

    }

    private List getFilesFromDir(String dir){
        //https://stackoverflow.com/questions/3953965/get-a-list-of-all-the-files-in-a-directory-recursive
        def list = []

        def dire = new File(dir)
        dire.eachFileRecurse (FileType.FILES) { file ->
            list << file
        }
        return list
    }

    private Path generateMetadataFile(){
        Path metadataFile = Files.createTempFile("","");
        metadataFile.write"Command Line: ${commandLine}\n"
        /*       "Container technology: ${containerTech}\n"
        metadataFile << "Container name: ${containerName}\n"
        metadataFile << "Container SHA256: ${containerSha}\n"
        metadataFile << */
        metadataFile << "UUID: ${uuid}\n"
        metadataFile << "Nextflow version: ${nextflowVersion}\n"

        return metadataFile
    }

    public void getCleanInputFiles(TraceRecord trace){
        def inputFiles = trace.getFmtStr("input")
        inputFiles= inputFiles.replaceAll("\\s","")//remove whitespace from the string
        List<String> inputList = Arrays.asList(inputFiles.split(';'));

        for (inputElem in inputList){
            if (!inputElem.contains(/work/)){
                this.inputFiles.add(inputElem.trim())
            }
        }
        log.info "Get Clean Input Files (data folder)"

    }

    public void getCleanOutputFiles(){
        if (!outDirFolder.equals("null")){ // because we miss the "concept" null when we add it on the map>>file>>map
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
        log.info "Get Clean Output Files (output folder)"
    }

    private void fileToBundle(Bundle bundle, Path filePath, String fileName, String folderName){
        if (filePath.isFile()){
            Path folderPath = bundle.getRoot().resolve(folderName);
            Path bundleFilePath = folderPath.resolve(fileName);
            Bundles.setStringValue(bundleFilePath, filePath.text);
            Files.copy(bundleFilePath, filePath, StandardCopyOption.REPLACE_EXISTING);
        }else if (filePath.isDirectory()){
            log.warn("The element: \"${fileName}\" is a directory")
            directoryToBundle(bundle,filePath,fileName,folderName)
        }
        log.info "File to bundle: ${fileName} to ${folderName}"
    }
    private void directoryToBundle(Bundle bundle, Path filePath, String fileName, String folderName){
        //TODO NEED TO BE IMPLEMENT!!!
    }
    public void saveBundle(Bundle bundle) {
        Path ro = Paths.get("${System.getProperty("user.dir")}/${roBundleName}.zip")

        /**
         * Remove the older zip if it exist
         */
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

    private List getOutDirFiles(String outDirFolder){
        if (outDirFolder!=null){
            log.info "OutDir folder exist"
            return getFilesFromDir(outDirFolder)
        }
        log.warn "OutDir folder NOT exist"
    }
    private void setAuthorInformation(Manifest manifest) {
        if (author!=null){
            Agent createdBy = new Agent(author);

            if (authorORCID!=null) {
                createdBy.setOrcid(URI.create(authorORCID));
            }else{
                createdBy.setOrcid(URI.create(orcidERROR));
            }
            manifest.setCreatedBy(createdBy);
        }else{
            // what to do if we dont have author on the manifest?? -> WARNING?
            //Now we dont "generate" the entity --> not ORCID allow neither
        }
        log.info "set Author info to RO manifest"

    }

    private String getFieldMap(Map map, String field){
        if(map.get(field) != "not_defined"){ //from map to string to file -> to map again... we miss the "concept" null, and it stay as a string value
            return map.get(field)
        }else{
            log.warn("The field ${field} is not configured on the project.")
        }
    }

    private void setAggregationManifest(Manifest manifest){
        LinkedList<PathMetadata> aggregationList = new LinkedList<PathMetadata>()
        manifest.setAggregates(aggregationList)

        log.info "Set Aggregation to RO manifest"
    }

    /*
    // TODO IT WILL BE PART OF THE AGENT
    private String getContainerTechnology(Map config) {
        if (config.containsKey("singularity") && config.get("singularity").getAt("enabled").toString().equals("true")) {
            singularityImage = true
            return "Singularity"
        } else if (config.containsKey("docker") && config.get("docker").getAt("enabled").toString().equals("true")) {
            dockerImage = true
            return "Docker"
        }
        //TODO capture if user uses CONDA ??
    }

    private String getContainerName(Map config){
        return config.getAt('process').getAt('container')
    }

    private String getContainerSHA256() {
        def cmd
        def duration = '10min'

         //decide the command line to get the SHA

        if (dockerImage == true) {
            //https://stackoverflow.com/questions/32046334/where-can-i-find-the-sha256-code-of-a-docker-image
            //give a diff sha256 the .repoDigest and the Id value --> TO CHECK
            cmd = "docker inspect --format='{{index .Id}}' ${containerName}"
        } else if (singularityImage == true) {
            //TODO does it work on macOX ??
            //--> DONT like it! it takes a while to digest the sha256
            //*-*- get the value when we pull the image
            //**** singularity pull --hash shub://vsoch/hello-world
            //**** Progress |===================================| 100.0%
            //**** Done. Container is at: /home/vanessa/ed9755a0871f04db3e14971bec56a33f.simg                      <-- file hash

            // use RGEGISTRY (singularity global client) : sregistry inspect $IMAGE -> version
            cmd = "sha256sum ./work/singularity/cbcrg-regressive-msa-v0.2.6.img | cut -d' ' -f1" //use container PATH
        }

        final max = Duration.of(duration).toMillis()
        final builder = new ProcessBuilder(['bash', '-c', cmd.toString()])
        final proc = builder.start()
        final err = new StringBuilder()
        //https://stackoverflow.com/questions/8149828/read-the-output-from-java-exec
        builder.redirectErrorStream(true)
        BufferedReader ine;
        ine = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        String line
        String containerAux
        while ((line = ine.readLine()) != null) {
            containerAux = line
            // remove the prefix for the sha256 information
            if (dockerImage) {
                containerAux = containerAux.substring(dockerSHAPrefix.size())
            } else {
                containerAux = containerAux.substring(singularitySHAPrefix.size())
            }
        }
        ine.close();
        proc.consumeProcessErrorStream(err)
        proc.waitForOrKill(max)
        def status = proc.exitValue()
        if (status != 0) {
            def msg = "Failed to get image info\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
        return containerAux
    }

    private void getMetadataInfo(){
        containerTech = getContainerTechnology(configMap)
        if (dockerImage || singularityImage) {
            containerName = getContainerName(configMap)
            containerSha = getContainerSHA256()
        }
    }*/
    private File generateScript(){
        File scriptFile= new File(commandLineFileName)
        scriptFile.write(commandLine)

        return scriptFile
    }
    private boolean removeFile(String fileName){
        return Files.deleteIfExists(Paths.get(fileName))
    }
    private File getLogInfo(){
        def today = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date())
        File logFile = new File(logFileName)
        logFile.append("${today} -- ${author} -- ${uuid} -- ${commandLine}\n")

        return logFile
    }
    public void setProvValues(Map provMap){
        //TODO modify variable to "null" concept?

        author = getFieldMap(provMap, 'author')
        authorORCID = getFieldMap(provMap,'orcid')
        commandLine= getFieldMap(provMap, 'commandLine')
        uuid= getFieldMap(provMap, 'uuid')
        nextflowVersion= getFieldMap(provMap, 'nfVersion').substring(nextflowVersionPrefix.size())
        outDirFolder = getFieldMap(provMap, 'outDir')
        baseDir = getFieldMap(provMap, 'baseDir')

    }
    public void printMap(Map map){
     for (Map.Entry<String, String> element : map.entrySet()){
         print   "-->>value: ..${element.getKey()}..  -- ${element.getValue()} \n"
     }
    }


}
