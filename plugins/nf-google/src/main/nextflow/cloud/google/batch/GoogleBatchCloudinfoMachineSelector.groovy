package nextflow.cloud.google.batch

import groovy.json.JsonSlurper
import groovy.transform.Immutable
import groovy.transform.Memoized

class GoogleBatchCloudinfoMachineSelector {

    private static final CLOUD_INFO_API = "https://cloudinfo.seqera.io/api/v1"

    // Some families CPUs are faster so this is a cost correction factor
    // for processes that request more than 2 CPUs or 2GB
    // https://cloud.google.com/compute/docs/cpu-platforms
    private static final familyCostCorrection = [
            'e2' : 1.0,   // Mix of processors, tend to be similar in performance to N1

            // INTEL
            'n1' : 1.0,   // Skylake, Broadwell, Haswell, Sandy Bridge, and Ivy Bridge ~2.7 Ghz
            'n2' : 0.85,  // Intel Xeon Gold 6268CL ~3.4 Ghz
            'c2' : 0.8,   // Intel Xeon Gold 6253CL ~3.8 Ghz
            'm1' : 1.0,   // Intel Xeon E7-8880V4 ~2.7 Ghz
            'm2' : 0.85,  // Intel Xeon Platinum 8280L ~3.4 Ghz
            'm3' : 0.85,  // Intel Xeon Platinum 8373C ~3.4 Ghz
            'a2' : 0.9,   // Intel Xeon Platinum 8273CL ~2.9 Ghz

            // AMD
            't2d': 1.0,   // AMD EPYC Milan ~2.7 Ghz
            'n2d': 1.0,   // AMD EPYC Milan ~2.7 Ghz
    ]

    // Families that will be use as default if Fusion is enabled but no list is provided
    // https://cloud.google.com/compute/docs/disks#local_ssd_machine_type_restrictions
    // LAST UPDATE 2023-02-17
    private static final defaultFamiliesWithSSD = ['n1', 'n2', 'n2d', 'c2', 'c2d', 'm3']

    @Immutable
    static class MachineType {
        String type
        String family
        float spotPrice
        float onDemandPrice
        int cpusPerVm
        int memPerVm
    }

    private static final jsonParser = new JsonSlurper()

    static String bestMachineType(int cpus, int memoryMB, String region, boolean spot, boolean localSSD, List<String> families) {
        final machineTypes = getAvailableMachineTypes(region)
        final memoryGB = Math.ceil(memoryMB / 1024) as int
        final cores = cpus > 1 ? cpus : 2

        // Use only families that can have a local SSD
        if (!families && localSSD)
            families = defaultFamiliesWithSSD

        // All types are valid if no families are defined, otherwise at least it has to start with one of the given values
        final matchMachineType = { t -> !families || families.find{t.startsWith(it)} }

        // find machines with enough resources and SSD local disk
        final validMachineTypes = machineTypes.findAll {
                    it.cpusPerVm >= cores &&
                    it.memPerVm >= memoryGB &&
                    matchMachineType(it.type)
        }.collect()

        final sortedByCost = validMachineTypes.sort {
            (it.cpusPerVm > 2 || it.memPerVm > 2 ? familyCostCorrection.get(it.family, 1.0) : 1.0) * (spot ? it.spotPrice : it.onDemandPrice)
        }

        return sortedByCost.first().type
    }

    @Memoized
    static List<MachineType> getAvailableMachineTypes(String region) {
        final json = "${CLOUD_INFO_API}/providers/google/services/compute/regions/${region}/products".toURL().text
        final data = jsonParser.parseText(json)
        data['products'].collect {
            new MachineType(
                    type: it.type,
                    family: it.type.split('-')[0],
                    spotPrice: it.spotPrice*.price.with { sum { it } / size() } as float,
                    onDemandPrice: it.onDemandPrice as float,
                    cpusPerVm: it.cpusPerVm as int,
                    memPerVm: it.memPerVm as int
            )
        }
    }
}
