package nextflow.cloud.gce;

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential;
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport;
import com.google.api.client.http.HttpTransport;
import com.google.api.client.json.JsonFactory;
import com.google.api.client.json.jackson2.JacksonFactory;
import com.google.api.services.compute.Compute;
import com.google.api.services.compute.model.*;

import java.io.IOException;
import java.math.BigInteger;
import java.security.GeneralSecurityException;
import java.util.*;

/**
 * Helper class for Google Compute Engine
 *
 * @author Vilmundur Pálmason <vilmundur@wuxinextcode.com>
 */
public class GceApiHelper {
    private static final String PROJECT_PREFIX = "https://www.googleapis.com/compute/v1/projects/";
    private final String project;
    private final String zone;
    private Compute compute;
    Random random = new Random();

    public GceApiHelper(String project, String zone) throws IOException, GeneralSecurityException {
        this.project = project;
        this.zone = zone;
        this.compute = createComputeService();
    }

    public Compute compute() {
        return compute;
    }

    public Compute createComputeService() throws IOException, GeneralSecurityException {
        HttpTransport httpTransport = GoogleNetHttpTransport.newTrustedTransport();
        JsonFactory jsonFactory = JacksonFactory.getDefaultInstance();

        GoogleCredential credential = GoogleCredential.getApplicationDefault();

        if (credential.createScopedRequired()) {
            credential =
                    credential.createScoped(Arrays.asList("https://www.googleapis.com/auth/cloud-platform"));
        }

        return new Compute.Builder(httpTransport, jsonFactory, credential)
                .setApplicationName("NextCode-Experiments/0.1")
                .build();
    }

    public String projectZonePrefix() {
        return PROJECT_PREFIX + project + "/zones/" + zone +"/";
    }

    /**
     * Full name of machine type
     * @param shortName  Short name such as "n1-standard-1"
     * @return Fully qualifie machine type
     */
    public String instanceType(String shortName) {
        return  projectZonePrefix() + "machineTypes/"+shortName;
    }

    /**
     * Full name of image.
     * @param imagePath including image project (e.g. "debian-cloud/global/images/debian-7-wheezy-v20150710" )
     * @return Fully qualified image name
     */
    public String imageName(String imagePath) {
        return PROJECT_PREFIX + imagePath;
    }

    public AttachedDisk setBootDisk(Instance instance, String imagePath) {
        AttachedDisk disk = new AttachedDisk();
        disk.setBoot(true);
        disk.setAutoDelete(true);
        disk.setType("PERSISTENT");
        AttachedDiskInitializeParams params = new AttachedDiskInitializeParams();
        // Assign the Persistent Disk the same name as the VM Instance.
        if (instance.getName() != null) {
            params.setDiskName(instance.getName());
        }
        // Specify the source operating system machine image to be used by the VM Instance.
        params.setSourceImage(imageName(imagePath));
        // Specify the disk type as Standard Persistent Disk
        params.setDiskType(projectZonePrefix()+ "diskTypes/pd-standard");
        disk.setInitializeParams(params);
        instance.setDisks(Collections.singletonList(disk));
        return disk;
    }

    public NetworkInterface setNetworkInterface(Instance inst) {
        NetworkInterface ifc = new NetworkInterface();
        ifc.setNetwork(PROJECT_PREFIX + project + "/global/networks/default");
        List<AccessConfig> configs = new ArrayList<>();
        AccessConfig config = new AccessConfig();
        config.setType("ONE_TO_ONE_NAT");
        config.setName("External NAT");
        configs.add(config);
        ifc.setAccessConfigs(configs);
        inst.setNetworkInterfaces(Collections.singletonList(ifc));
        return ifc;
    }

    public String randomName(String baseName) {
        return baseName + randomName();
    }

    public String randomName() {
        byte[] bytes = new byte[5];
        random.nextBytes(bytes);
        return new BigInteger(bytes).abs().toString(16);
    }

    public Metadata createMetadata(String ... tagVal) {
        Metadata metadata = new Metadata();

        List<Metadata.Items> items = new ArrayList<>();
        for (int i=0;i<tagVal.length -1; i+=2) {
            Metadata.Items it = new Metadata.Items();
            it.set(tagVal[i],tagVal[i+1]);
            items.add(it);
        }
        metadata.setItems(items);
        System.out.println(metadata);
        return metadata;
    }

    public Operation.Error blockUntilComplete(Iterable<Operation> ops, long timeoutMs) throws InterruptedException, IOException {
        long start = System.currentTimeMillis();
        for (Operation op: ops) {
            Operation.Error result = blockUntilComplete(op,timeoutMs - (System.currentTimeMillis() - start));
            if (result != null) return null;
        }
        return null;
    }

    public Operation.Error blockUntilComplete(Operation operation, long timeoutMs) throws InterruptedException, IOException {
        long start = System.currentTimeMillis();
        final long pollInterval = 5 * 1000;
        String zone = operation.getZone();  // null for global/regional operations
        if (zone != null) {
            String[] bits = zone.split("/");
            zone = bits[bits.length - 1];
        }
        String status = operation.getStatus();
        String opId = operation.getName();
        while (operation != null && !status.equals("DONE")) {
            Thread.sleep(pollInterval);
            long elapsed = System.currentTimeMillis() - start;
            if (elapsed >= timeoutMs) {
                throw new InterruptedException("Timed out waiting for operation to complete");
            }
            if (zone != null) {
                Compute.ZoneOperations.Get get = compute.zoneOperations().get(project, zone, opId);
                operation = get.execute();
            } else {
                Compute.GlobalOperations.Get get = compute.globalOperations().get(project, opId);
                operation = get.execute();
            }
            if (operation != null) {
                status = operation.getStatus();
            }
        }
        return operation == null ? null : operation.getError();
    }

    public Image lookupImage(String imagePath) throws IOException {
        return compute.images().get(project,imageName(imagePath)).execute();
    }

    public MachineType lookupMachineType(String machineType) throws IOException {
        return compute.machineTypes().get(project,zone,machineType).execute();
    }

    public String instanceIdToPrivateDNS(String instanceId) {
        return instanceId+".c."+project+".internal";
    }

    public String publicIpToDns(String ip) {
        if (ip == null) return null;
        String[] parts = ip.split("\\.");
        if (parts.length != 4) throw new IllegalArgumentException("Expected IPv4 Public IP address instead of '"+ip+"'");

        // TODO: Is this domain name stable ?
        return parts[3]+"."+parts[2]+"."+parts[1]+"."+parts[0]+".bc.googleusercontent.com";
    }

    /**
     * Check if value is valid as a label value as specified here: https://cloud.google.com/compute/docs/labeling-resources
     * @return null if valid or error message
     */
    public String validateLabelValue(String value) {
        if (value == null) return null;

        if (!value.matches("[a-z0-9-_]*")) {
            return "Value must consist of lowercase letters, numbers, underscores and dashes only";
        }
        if (value.length() > 63) {
            return "Value exceeds maximum length of 63";
        }
        return null;
    }


}
