package nextflow.cloud.gce.pipelines

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.v2alpha1.Genomics
import com.google.api.services.genomics.v2alpha1.model.Action
import com.google.api.services.genomics.v2alpha1.model.Disk
import com.google.api.services.genomics.v2alpha1.model.Mount
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import com.google.api.services.genomics.v2alpha1.model.Resources
import com.google.api.services.genomics.v2alpha1.model.ServiceAccount
import com.google.api.services.genomics.v2alpha1.model.VirtualMachine
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

@Slf4j
@CompileStatic
class GooglePipelinesHelper {

    static final String SCOPE_CLOUD_PLATFORM = "https://www.googleapis.com/auth/cloud-platform"
    static final List<String> ENV_VAR_TO_INCLUDE = ["NXF_DEBUG"]

    enum ActionFlags {
        FLAG_UNSPECIFIED,
        IGNORE_EXIT_STATUS,
        RUN_IN_BACKGROUND,
        ALWAYS_RUN,
        ENABLE_FUSE,
        PUBLISH_EXPOSED_PORTS,
        DISABLE_IMAGE_PREFETCH,
        DISABLE_STANDARD_ERROR_CAPTURE
    }

    static String sanitizeName(String name) {
        name.replaceAll(/[^a-zA-Z0-9\-_]+/, '-').take(63)
    }

    static Genomics createGenomicClient() {

        def credentials = GoogleCredential.applicationDefault

        if (credentials.createScopedRequired()) {
            credentials = credentials.createScoped([SCOPE_CLOUD_PLATFORM])
        }

        new Genomics.Builder(GoogleNetHttpTransport.newTrustedTransport(), JacksonFactory.defaultInstance, credentials)
                .setApplicationName("Nextflow GooglePipelinesExecutor")
                .build()
    }

    static Map<String,String> getEnvironment() {
        def ret = [:]
        def env = System.getenv()

        env.each { kv ->
            if(ENV_VAR_TO_INCLUDE.contains(kv.key))
                ret << kv
        }
        ret
    }

    static Action createAction(String name, String imageUri, List<String> commands, List<Mount> mounts, List<ActionFlags> flags = [], String entrypoint = null) {
        new Action()
                .setName(name)
                .setImageUri(imageUri)
                .setCommands(commands)
                .setMounts(mounts)
                .setFlags(flags.collect{flag -> flag.toString()})
                .setEntrypoint(entrypoint)
                .setEnvironment(getEnvironment())
    }

    static Pipeline createPipeline(List<Action> actions, Resources resources) {
        new Pipeline().setActions(actions).setResources(resources)
    }

    //TODO: Do we want to configure this via nextflow config?
    static Resources configureResources(String instanceType, String projectId, String zone, String diskName, List<String> scopes = null,boolean preEmptible = false) {

        def disk = new Disk()
        disk.setName(diskName)

        def serviceAccount = new ServiceAccount()
        if (scopes)
            serviceAccount.setScopes(scopes)

        def vm = new VirtualMachine()
                .setMachineType(instanceType)
                .setDisks([disk])
                .setServiceAccount(serviceAccount)
                .setPreemptible(preEmptible)


        new Resources()
                .setProjectId(projectId)
                .setZones([zone])
                .setVirtualMachine(vm)
    }

    static Mount configureMount(String diskName, String mountPath, boolean readOnly = false) {
        new Mount().setDisk(diskName).setPath(mountPath).setReadOnly(readOnly)
    }
}