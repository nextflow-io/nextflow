package nextflow.prov

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import nextflow.cli.Launcher
import nextflow.trace.TraceRecord
import nextflow.util.Duration
import org.apache.taverna.robundle.Bundle
import org.apache.taverna.robundle.Bundles
import org.apache.taverna.robundle.manifest.Agent
import org.apache.taverna.robundle.manifest.Manifest
import org.apache.taverna.robundle.manifest.PathMetadata

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.text.SimpleDateFormat

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

    private boolean dockerImage
    private boolean singularityImage
    private String containerTech
    private String containerName
    private String containerSha
    private String dockerSHAPrefix ="sha256:"
    private String singularitySHAPrefix=""
    private String commandLine=System.getenv('NXF_CLI')
    private String uuid
    private String nextflowVersionPrefix = "nextflow version "
    private String nextflowVersion

    private String commandLineFileName = "commandLine.txt"
    private String enviromentFileName = "enviroment.txt"

    private String orcidERROR= "**ORCID_not_provided**"

    private List<String> inputFiles =[]
    private List<String> outputFiles =[]

    ResearchObjectGenerator(){

    }
    public Bundle generateROBundle(){
        return Bundles.createBundle();
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
    }

    public void setManifest(Bundle bundle, Map configManifest){
        // https://github.com/apache/incubator-taverna-language/blob/master/taverna-robundle/src/test/java/org/apache/taverna/robundle/manifest/TestManifestJSON.java
        Manifest manifest = bundle.getManifest();

        this.setAuthorInformation(manifest,configManifest)

        manifest.setId(URI.create("/"))

        setAggregationManifest(manifest)
    }

    public void generateLogFile(Bundle bundle, String author ){
        File logFile = getLogInfo(author)
        fileToBundle(bundle,Paths.get(logFile.path),logFileName, metadataFolder)
    }

    public void addProvenanceFile(Bundle bundle){
        File file = provDocumentToFile(provenanceFileName)
        fileToBundle(bundle, Paths.get(file.getPath()), provenanceFileName,metadataFolder)
        boolean result = Files.deleteIfExists(Paths.get(provenanceFileName))
    }

    public void generateDataFolder(Bundle bundle){
        inputFiles.unique()
        for (file in inputFiles){
            Path auxPath = Paths.get(file)
            fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), dataFolderName)
        }
        //TODO generate sha256 HERE ??
    }
    public Path generateOutputFolder(Bundle bundle){
        outputFiles.unique()
        for (file in outputFiles){
            Path auxPath = Paths.get(file)
            fileToBundle(bundle, auxPath, auxPath.getFileName().toString(), outputFolderName)
        }
        //TODO generate sha256 HERE??
    }

    public void generateWorkflowFolder(Bundle bundle, String baseDir){
        def result = getFilesFromDir(baseDir).tokenize('\n')
        for (element in result){
            Path auxPath = Paths.get("${baseDir}/${element}")
            fileToBundle(bundle, auxPath,auxPath.getFileName().toString(),workflowFolderName)
        }
    }

    public void generateSnapshot(Bundle bundle){

    }

    public void generateMetadataFolder(Bundle bundle, Map configMap){
        Path metadataFile = generateMetadataFile(configMap)
        fileToBundle(bundle, metadataFile,metadataFileName,metadataFolder)
    }

    private String getFilesFromDir(String dir){
        def proc = "ls ${dir}".execute()
        def b = new StringBuffer()
        proc.consumeProcessErrorStream(b)
        def resultado = proc.text
        return resultado
    }

    private Path generateMetadataFile(Map configMap){
        getMetadataInfo(configMap)
        return writeMetadataIntoFile()
    }

    private File provDocumentToFile(String provenanceFileName){
        File file = new File(provenanceFileName)
        return file
    }

    public void getCleanInputFiles(TraceRecord trace){
        def inputFiles = trace.getFmtStr("input")
        List<String> inputList = Arrays.asList(inputFiles.split(";"));

        for (inputElem in inputList){
            if (!inputElem.contains(/work/)){
                this.inputFiles.add(inputElem.trim())
            }
        }
    }

    public void getCleanOutputFiles(Map manifest){
        def outDirFolder = getOutdirFolder(manifest)
        if (outDirFolder !=null){
            def outDirFiles = getOutDirFiles(outDirFolder)
            File fileAux = new File ("${outDirFolder}/${outDirFiles[0].toString()}")
            String filePath = fileAux.absolutePath.substring(0,fileAux.absolutePath.lastIndexOf(File.separator));
            for (element in outDirFiles){
                String aux = "${filePath}/${element}"
                outputFiles.add(aux)
            }
        }else{
            log.warn("You need to specify the output directory inside nextflow.config/params for the RO zip file")
        }
    }

    private void fileToBundle(Bundle bundle, Path filePath, String fileName, String folderName){
        if (filePath.isFile()){
            Path folderPath = bundle.getRoot().resolve(folderName);
            Path bundleFilePath = folderPath.resolve(fileName);
            Bundles.setStringValue(bundleFilePath, filePath.text);
            Files.copy(bundleFilePath, filePath, StandardCopyOption.REPLACE_EXISTING);
        }else if (filePath.isDirectory()){
            log.warn("The element: \"${fileName}\" is a directory")
        }
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

    private String getOutdirFolder(Map manifest){
        def paramsMap = manifest.getAt('params')
        return paramsMap.getAt('outdir')
    }
    private List getOutDirFiles(String outDirFolder){
        if (outDirFolder!=null){
            return getFilesFromDir(outDirFolder).tokenize('\n')
        }
    }
    private void setAuthorInformation(Manifest manifest, Map configManifest) {
        def author = getAuthor(configManifest)
        if (author!=null){
            Agent createdBy = new Agent(author);

            def authorORCID = getAuthorORCID(configManifest)
            if (authorORCID!=null) {                   //TODO allow author ORCID on nextflow.config.manifest ?
                createdBy.setOrcid(URI.create(authorORCID));
            }else{
                createdBy.setOrcid(URI.create(orcidERROR));
            }
            manifest.setCreatedBy(createdBy);
        }else{
            // what to do if we dont have author on the manifest?? -> WARNING?
            //Now we dont "generate" the entity --> not ORCID allow neither
        }
    }

    private String getAuthor(Map manifestConfig) {
        if(manifestConfig.author!=null){
            return manifestConfig.author
        }else {
            log.warn("You need to specify the Author's name  inside nextflow.config/manifest for the RO zip file")
        }
    }

    private String getAuthorORCID(Map manifestConfig) {
        if(manifestConfig.ORCID!=null){
            return manifestConfig.ORCID
        }else {
            log.warn("You need to specify the Author's ORCID  inside nextflow.config/manifest for the RO zip file")
        }
    }

    private void setAggregationManifest(Manifest manifest){
        //TODO fill list
        LinkedList<PathMetadata> aggregationList = new LinkedList<PathMetadata>()
        PathMetadata metaPath = new PathMetadata("URI_Example")        // set URI

        aggregationList.add(metaPath)
        manifest.setAggregates(aggregationList)
    }

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
        /**
         * decide the command line to get the SHA
         */
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

    private void getMetadataInfo(Map configMap){
        containerTech = getContainerTechnology(configMap)
        if (dockerImage || singularityImage) {
            containerName = getContainerName(configMap)
            containerSha = getContainerSHA256()
        }
        uuid= "" //TODO GET UUID --> "${session.getUniqueId().toString()}"
        nextflowVersion= Launcher.getVersion().substring(nextflowVersionPrefix.size())

    }

    private Path  writeMetadataIntoFile(){
        Path metadataFile = Files.createTempFile("","");
        metadataFile.write "Container technology: ${containerTech}\n"
        metadataFile << "Container name: ${containerName}\n"
        metadataFile << "Container SHA256: ${containerSha}\n"
        metadataFile << "Command Line: ${commandLine}\n"
        metadataFile << "UUID: ${uuid}\n"
        metadataFile << "Nextflow version: ${nextflowVersion}\n"

        return metadataFile
    }

    private File getLogInfo(String author){
        def today = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date())
        File logFile = new File(logFileName)
        logFile.append("${today} -- ${author} -- ${uuid} -- ${containerSha}\n\t${commandLine}\n")

        return logFile
    }

    public void printMap(Map map){
     for (Map.Entry<String, String> element : map.entrySet()){
         print   "-->>value: ${element.getKey()}  -- ${element.getValue()} \n"
     }
    }
}
