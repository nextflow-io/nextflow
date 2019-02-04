/*
 * Copyright 2018, WuxiNextcode
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.cloud.google

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.compute.Compute
import com.google.api.services.compute.model.*
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovyjarjarcommonscli.MissingArgumentException
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper

import java.nio.file.Path
import java.security.GeneralSecurityException

/**
 * Helper class for Google Compute Engine.
 *
 * @author Vilmundur Pálmason <vilmundur@wuxinextcode.com>
 */
@CompileStatic
class GceApiHelper {
    private static final String PROJECT_PREFIX = "https://www.googleapis.com/compute/v1/projects/"
    public static final String GAC_ENV = "GOOGLE_APPLICATION_CREDENTIALS"

    final String project
    final String zone
    final Compute compute

    Random random = new Random()


    /**
     * only for testing purpose -- do not use
     */
    GceApiHelper(String project, String zone, Compute compute) {
        this.project = project
        this.zone = zone
        this.compute = compute
    }

    GceApiHelper(String project, String zone) throws IOException, GeneralSecurityException {
        this.project = project ?: readProject()
        this.zone = zone ?: readZone()
        this.compute = createComputeService(GoogleCredential.getApplicationDefault())
    }

    static Compute createComputeService(GoogleCredential credential) throws IOException, GeneralSecurityException {
        def httpTransport = GoogleNetHttpTransport.newTrustedTransport()
        def jsonFactory = JacksonFactory.getDefaultInstance()

        if (credential.createScopedRequired()) {
            credential =
                    credential.createScoped(Arrays.asList("https://www.googleapis.com/auth/cloud-platform"))
        }

        new Compute.Builder(httpTransport, jsonFactory, credential)
                .setApplicationName("nextflow/${Const.APP_VER}")
                .build()
    }

    String projectZonePrefix() {
        "${PROJECT_PREFIX}$project/zones/$zone/"
    }

    /**
     * Full name of machine type
     * @param shortName Short name such as "n1-standard-1"
     * @return Fully qualified machine type
     */
    String instanceType(String shortName) {
        "${projectZonePrefix()}machineTypes/$shortName"
    }

    /**
     * Full name of image.
     * @param imagePath including image project (e.g. "debian-cloud/global/images/debian-7-wheezy-v20150710" )
     * @return Fully qualified image name
     */
    static String imageName(String imagePath) {
        PROJECT_PREFIX + imagePath
    }

    AttachedDisk createBootDisk(String name, String imagePath) {
        def disk = new AttachedDisk()
        disk.setBoot(true)
        disk.setAutoDelete(true)
        disk.setType("PERSISTENT")
        def params = new AttachedDiskInitializeParams()
        // Assign the Persistent Disk the same name as the VM Instance.
        if (name != null) {
            params.setDiskName(name)
        }
        // Specify the source operating system machine image to be used by the VM Instance.
        params.setSourceImage(imageName(imagePath))
        // Specify the disk type as Standard Persistent Disk
        params.setDiskType(projectZonePrefix() + "diskTypes/pd-standard")
        disk.setInitializeParams(params)
        disk
    }

    NetworkInterface createNetworkInterface() {
        def ifc = new NetworkInterface()
        ifc.setNetwork("${PROJECT_PREFIX}${project}/global/networks/default")
        List<AccessConfig> configs = []
        def config = new AccessConfig()
        config.setType("ONE_TO_ONE_NAT")
        config.setName("External NAT")
        configs.add(config)
        ifc.setAccessConfigs(configs)
        ifc
    }

    String randomName(String baseName) {
        baseName + randomName()
    }

    String randomName() {
        def bytes = new byte[5]
        random.nextBytes(bytes)
        new BigInteger(bytes).abs().toString(16)
    }

    List<Instance> getInstanceList(String filter) {
        def listRequest = compute.instances().list(project,zone)
        listRequest.setFilter(filter)
        listRequest.execute().getItems()
    }

    /**
     * Block until all operations are complete or if any results in an error.
     */
    Operation.Error blockUntilComplete(Iterable<Operation> ops, long timeoutMs = 20000, long pollingIntervalMs = 5000) {
        long start = System.currentTimeMillis()
        for (Operation op : ops) {
            Operation.Error result = blockUntilComplete(op, timeoutMs - (System.currentTimeMillis() - start),pollingIntervalMs)
            if (result != null) return result
        }
        return null
    }

    Operation.Error blockUntilComplete(Operation operation, long timeoutMs = 10000, long pollingIntervalMs = 5000) {
        def start = System.currentTimeMillis()

        def opZone = operation.getZone()  // null for global/regional operations

        if (opZone != null) {
            opZone = opZone.split("/").last()
        }
        def status = operation.getStatus()
        def opId = operation.getName()
        while (operation != null && status != "DONE") {
            Thread.sleep(pollingIntervalMs)
            long elapsed = System.currentTimeMillis() - start
            if (elapsed >= timeoutMs) {
                throw new InterruptedException("Timed out waiting for operation to complete")
            }
            if (opZone != null) {
                Compute.ZoneOperations.Get get = compute.zoneOperations().get(project, opZone, opId)
                operation = get.execute()
            } else {
                Compute.GlobalOperations.Get get = compute.globalOperations().get(project, opId)
                operation = get.execute()
            }
            status = operation?.getStatus()
        }
        operation?.getError()
    }

    Image lookupImage(String imageId) {

        def img = imageId.split("/")[0]

        def list =compute.images().list(img).execute()

        list.getItems().find {it.getSelfLink().endsWith(imageId)}
    }

    MachineType lookupMachineType(String machineType) {
        try {
            compute.machineTypes().get(project, zone, machineType).execute()
        } catch(GoogleJsonResponseException re) {
            if(re.details.getCode() == 404)
                return null
            else throw re
        }
    }

    String instanceIdToPrivateDNS(String instanceId) {
        "${instanceId}.c.${project}.internal"
    }

    String publicIpToDns(String ip) {
        if (ip == null) return null
        String[] parts = ip.split("\\.")
        if (parts.length != 4) throw new IllegalArgumentException("Expected IPv4 Public IP address instead of '" + ip + "'")

        parts.reverse().join(".") + ".bc.googleusercontent.com"
    }

    /**
     * Check if value is valid as a label value as specified here: https://cloud.google.com/compute/docs/labeling-resources
     * @return null if valid or error message
     */
    String validateLabelValue(String value) {
        if (value == null)
            return null

        if (!value.matches("[a-z0-9-_]*")) {
            return "Value must consist of lowercase letters, numbers, underscores and dashes only"
        }
        if (value.length() > 63) {
            return "Value exceeds maximum length of 63"
        }
        return null
    }

    def setStartupScript(Instance instance, String script) {
        addMetadataItem(instance, "startup-script", script)
    }

    def setShutdownScript(Instance instance, String script) {
        addMetadataItem(instance, "shutdown-script", script)
    }

    def addMetadataItem(Instance instance, String key, String value) {
        def metadata = instance.getMetadata() ?: new Metadata()
        List<Metadata.Items> items = metadata.getItems() ?: []

        Metadata.Items item = new Metadata.Items()
        item.setKey(key)
        item.setValue(value)
        items.add(item)
        metadata.setItems(items)
        instance.setMetadata(metadata)
    }

    String readGoogleMetadata(String meta) {
        try {
            def payload = "http://metadata/computeMetadata/v1/$meta".toURL().getText(requestProperties: ['Metadata-Flavor': 'Google'])
            if(!payload || payload?.isEmpty()) {
                throw new AbortOperationException("Could not find google metadata: $meta")
            } else {
                return payload
            }
        } catch (Exception e) {
            throw new AbortOperationException("Cannot read Google metadata $meta: ${e.getClass()}(${e.getMessage()})", e)
        }
    }

    String readProject() {
        if(isCredentialLocationDefined()) {
            def cred = new JsonSlurper().parseText(getCredentialsFile())
            def project = cred['project_id']
            if(!project) throw new AbortOperationException("Could not read project from Credential file")
            return project
        }
        readGoogleMetadata('project/project-id')
    }

    String readZone() {
        readGoogleMetadata('instance/zone').split("/").last()
    }

    String readInstanceId() {
        readGoogleMetadata('instance/name')
    }

    Scheduling createScheduling(boolean preemptible) {
        new Scheduling().setPreemptible(preemptible)
    }

    String getCredentialsFile() {
        String credFileLocation = credentialFileLocation()

        if(!credFileLocation)
            throw new MissingArgumentException("$GAC_ENV is not defined in your environment" )

        Path credFile = FileHelper.asPath(credFileLocation)
        if (credFile && credFile.exists() ) {
            credFile.toFile().text
        } else {
            throw new FileNotFoundException("Could not find Google credentials file '$credFileLocation'")
        }
    }

    boolean isCredentialLocationDefined() {
        credentialFileLocation()
    }

    String credentialFileLocation() {
        System.getenv(GAC_ENV)
    }
}
