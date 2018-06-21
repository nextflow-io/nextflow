package nextflow.prov

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cli.Launcher
import nextflow.util.Duration
import org.apache.taverna.robundle.Bundle
import org.apache.taverna.robundle.Bundles
import org.apache.taverna.robundle.manifest.Agent
import org.apache.taverna.robundle.manifest.Manifest
import org.apache.taverna.robundle.manifest.PathMetadata
import org.openprovenance.prov.model.Document

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.nio.file.StandardOpenOption
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
    // .ro
        //bundle manifest.xml

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

    private String commandLineFileName = "commandLine.txt"
    private String enviromentFileName = "enviroment.txt"


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
        fileToBundle(bundle,Paths.get(logFile.path),logFileName, "")
    }

    public void addProvenanceFile(Bundle bundle){
        File file = provDocumentToFile(provenanceFileName)
        fileToBundle(bundle, Paths.get(file.getPath()), provenanceFileName,metadataFolder)
        boolean result = Files.deleteIfExists(Paths.get(provenanceFileName))
    }

    public void generateMetadataFile(Bundle bundle, Map configMap){
        getMetadataInfo(configMap)
        Path metadataFile = writeMetadataIntoFile()
        fileToBundle(bundle, metadataFile,metadataFileName,metadataFolder)
    }

    public void generateSnapshot(Bundle bundle){

    }

    private File provDocumentToFile(String provenanceFileName){
        File file = new File(provenanceFileName)
        return file
    }

    private void fileToBundle(Bundle bundle, Path filePath, String fileName, String folderName){
        // Get the inputs
        Path folderPath = bundle.getRoot().resolve(folderName);

        // Get an input port:
        Path bundleFilePath = folderPath.resolve(fileName);

        // Setting a string value for the input port:
        Bundles.setStringValue(bundleFilePath, filePath.text);

        // Or Java 7 style
        Files.copy(bundleFilePath, filePath, StandardCopyOption.REPLACE_EXISTING);
    }

    public void saveBundle(Bundle bundle){
        Path zip = Files.createTempFile(roBundleName, ".zip");
        Bundles.closeAndSaveBundle(bundle, zip);
        log.info "RO Bundle saved to ${zip}"
    }

    private void setAuthorInformation(Manifest manifest, Map configManifest) {
        def author = getAuthor(configManifest)
        Agent createdBy = new Agent(author);

        def authorORCID = getAuthorORCID(configManifest)
        if (!authorORCID.isEmpty()) {                   //TODO allow author ORCID on nextflow.config.manifest ?
            createdBy.setOrcid(URI.create(authorORCID));
        }
        manifest.setCreatedBy(createdBy);
    }

    private String getAuthor(Map manifestConfig) {
        return manifestConfig.author
    }

    private String getAuthorORCID(Map manifestConfig) {
        return manifestConfig.ORCID
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
    }

    //TODO why it need to be public
    public String getContainerName(Map config){
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
            cmd = "sha256sum ./work/singularity/cbcrg-regressive-msa-v0.2.6.img | cut -d' ' -f1" //${containerName}
        }

        final max = Duration.of(duration).toMillis()
        final builder = new ProcessBuilder(['bash', '-c', cmd.toString()])
        //builder.directory(storePath.toFile())
        //builder.environment().remove('SINGULARITY_PULLFOLDER')
        final proc = builder.start()
        final err = new StringBuilder()
        //*********START read proc stdout       https://stackoverflow.com/questions/8149828/read-the-output-from-java-exec
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
        //*******END read proc stdout
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
        commandLine= System.getenv('NXF_CLI')
        uuid= "" //TODO GET UUID --> "${session.getUniqueId().toString()}"
        nextflowVersion= Launcher.getVersion().substring(nextflowVersionPrefix.size())

        /*
        println "**** METADATA info ****\n"+
                "containerTech: ${containerTech}\n"+
                "dockerImage: ${dockerImage}\tsingularityImage: ${singularityImage}\n"+
                "containerName: ${containerName}\n"+
                "containerSha:  ${containerSha}\n"+
                "commandLine:  ${commandLine}\n"+
                "uuid:  ${uuid}\n"+
                "nextflowVersion:  ${nextflowVersion}\n"+
                "****        ****"
       */
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
