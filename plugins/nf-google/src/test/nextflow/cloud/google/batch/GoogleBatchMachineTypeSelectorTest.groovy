package nextflow.cloud.google.batch

import nextflow.cloud.google.batch.GoogleBatchMachineTypeSelector.MachineType
import spock.lang.Specification

class GoogleBatchMachineTypeSelectorTest extends Specification {

    static final MACHINE_TYPES = [
            new MachineType(type: 'e2-type01', family: 'e2', 'spotPrice': 0.001, 'onDemandPrice': 0.01, 'cpusPerVm': 1, 'memPerVm': 1),
            new MachineType(type: 'n1-type02', family: 'n1', 'spotPrice': 0.002, 'onDemandPrice': 0.02, 'cpusPerVm': 2, 'memPerVm': 2),
            new MachineType(type: 'e2-type03', family: 'e2', 'spotPrice': 0.010, 'onDemandPrice': 0.05, 'cpusPerVm': 4, 'memPerVm': 4),
            new MachineType(type: 'n2-type04', family: 'n2', 'spotPrice': 0.011, 'onDemandPrice': 0.15, 'cpusPerVm': 4, 'memPerVm': 4),
            new MachineType(type: 'e2-type05', family: 'e2', 'spotPrice': 0.020, 'onDemandPrice': 0.20, 'cpusPerVm': 6, 'memPerVm': 6),
            new MachineType(type: 'n1-type06', family: 'n1', 'spotPrice': 0.025, 'onDemandPrice': 0.25, 'cpusPerVm': 6, 'memPerVm': 6),
            new MachineType(type: 'm1-type07', family: 'm1', 'spotPrice': 0.030, 'onDemandPrice': 0.30, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'm2-type08', family: 'm2', 'spotPrice': 0.036, 'onDemandPrice': 0.35, 'cpusPerVm': 8, 'memPerVm': 8),
            new MachineType(type: 'n2-type09', family: 'n2', 'spotPrice': 0.040, 'onDemandPrice': 0.40, 'cpusPerVm': 10, 'memPerVm': 10),
            new MachineType(type: 'c2-type10', family: 'c2', 'spotPrice': 0.045, 'onDemandPrice': 0.45, 'cpusPerVm': 10, 'memPerVm': 10),
    ]

    def 'should select best machine type'() {
        given:
        final selector = Spy(GoogleBatchMachineTypeSelector) {
            getAvailableMachineTypes(REGION) >> MACHINE_TYPES
        }
        expect:
        selector.bestMachineType(CPUS, MEM, REGION, SPOT, SSD, FAMILIES) == EXPECTED

        where:
        CPUS | MEM  | REGION | SPOT  | SSD   | FAMILIES                   | EXPECTED
        1    | 1000 | 'reg'  | true  | false | null                       | 'e2-type01'
        1    | 1000 | 'reg'  | false | true  | null                       | 'n1-type02'
        4    | 4000 | 'reg'  | false | false | []                         | 'e2-type03'
        4    | 4000 | 'reg'  | true  | false | []                         | 'n2-type04'
        4    | 4000 | 'reg'  | false | true  | []                         | 'n2-type04'
        6    | 6000 | 'reg'  | true  | false | null                       | 'e2-type05'
        6    | 6000 | 'reg'  | true  | true  | null                       | 'n1-type06'
        6    | 6000 | 'reg'  | true  | false | ['n1-*', 'm1-*']           | 'n1-type06'
        8    | 8000 | 'reg'  | true  | false | null                       | 'm1-type07'
        8    | 8000 | 'reg'  | false | false | ['m?-*', 'c2-*']           | 'm2-type08'
        8    | 8000 | 'reg'  | false | false | ['m1-type07', 'm2-type66'] | 'm1-type07'


    }

    def 'should not select a machine type'() {
        given:
        final selector = Spy(GoogleBatchMachineTypeSelector) {
            getAvailableMachineTypes(REGION) >> MACHINE_TYPES
        }
        when:
        selector.bestMachineType(CPUS, MEM, REGION, SPOT, SSD, FAMILIES)

        then:
        thrown(NoSuchElementException)

        where:
        CPUS | MEM  | REGION | SPOT | SSD   | FAMILIES
        8    | 9000 | 'reg'  | true | true  | ['s2-*', 'e2-*']
        12   | 1000 | 'reg'  | true | false | null

    }

    def 'should parse Seqera cloud info API'() {
        when:
        GoogleBatchMachineTypeSelector.INSTANCE.getAvailableMachineTypes("europe-west2")

        then:
        noExceptionThrown()
    }
}
