package nextflow.cloud.google.batch

import groovy.json.JsonSlurper
import groovy.transform.Immutable
import groovy.transform.Memoized

class GoogleBatchCloudinfoMachineSelector {

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

    // Some families CPUs are faster so this is
    // a cost correction factor
    private static final familyCostCorrection = [
            'n1': 1.0,   // Skylake, Broadwell, Haswell, Sandy Bridge, and Ivy Bridge
            'n2': 0.85,  // Cascade Lake and Ice Lake
            'c2': 0.9,   // Cascade Lake
            'm3': 0.85   // Ice Lake
    ]

    // Using only Intel with local SSD
    private static final validFamilies = ['n1', 'n2', 'c2', 'm3']

    static String bestMachineType(int cpus, int memoryMB, String region) {
        final machineTypes = getAvailableMachineTypes(region)
        final memoryGB = ((memoryMB / 1024) as int) + 1

        // find machines with enough resources and SSD local disk
        final validMachineTypes = machineTypes.findAll {
                    it.cpusPerVm >= cpus &&
                    it.memPerVm >= memoryGB &&
                    it.family in validFamilies }.collect()

        // TODO spot vs on-demand
        final sortedByCost = validMachineTypes.sort {
            familyCostCorrection[it.family] * it.spotPrice
        }

        return sortedByCost.first().type
    }

    @Memoized
    static List<MachineType> getAvailableMachineTypes(String region) {
        final json = "https://cloudinfo.seqera.io/api/v1/providers/google/services/compute/regions/${region}/products".toURL().text
        final data = jsonParser.parseText(json)
        data['products'].collect { new MachineType(
                type: it.type,
                family: it.type.split('-')[0],
                spotPrice: it.spotPrice*.price.with { sum{it} / size() } as float,
                onDemandPrice: it.onDemandPrice as float,
                cpusPerVm:  it.cpusPerVm as int,
                memPerVm: it.memPerVm as int
        )}
    }
}
