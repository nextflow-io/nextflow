package nextflow.processor

import spock.lang.Specification

class PublishOffloadConfigTest extends Specification {

    def "should initialize with default values when config is empty"() {
        given:
        Map config = [:]

        when:
        def publishOffloadConfig = new PublishOffloadConfig(config)

        then:
        publishOffloadConfig.enable == false
        publishOffloadConfig.useFusion == false
        publishOffloadConfig.batchSize == 10
        publishOffloadConfig.maxParallel == 10
    }

    def "should set values correctly based on the config map"() {
        given:
        Map config = [
            enable     : true,
            useFusion  : true,
            batchSize  : 20,
            maxParallel: 15
        ]

        when:
        def publishOffloadConfig = new PublishOffloadConfig(config)

        then:
        publishOffloadConfig.enable == true
        publishOffloadConfig.useFusion == true
        publishOffloadConfig.batchSize == 20
        publishOffloadConfig.maxParallel == 15
    }

    def "should enforce batchSize and maxParallel minimum value of 1"() {
        given:
        Map config = [
            batchSize  : -10,
            maxParallel: -5
        ]

        when:
        def publishOffloadConfig = new PublishOffloadConfig(config)

        then:
        publishOffloadConfig.batchSize == 1
        publishOffloadConfig.maxParallel == 1
    }

    def "should set maxParallel to batchSize if maxParallel is not specified"() {
        given:
        Map config = [
            batchSize: 25
        ]

        when:
        def publishOffloadConfig = new PublishOffloadConfig(config)

        then: "maxParallel is set equal to batchSize"
        publishOffloadConfig.batchSize == 25
        publishOffloadConfig.maxParallel == 25
    }

    def "should respect values for optional parameters when explicitly set in the config map"() {
        given:
        Map config = [
            enable   : false,
            useFusion: true
        ]

        when:
        def publishOffloadConfig = new PublishOffloadConfig(config)

        then:
        publishOffloadConfig.enable == false
        publishOffloadConfig.useFusion == true
        publishOffloadConfig.batchSize == 10 // default value
        publishOffloadConfig.maxParallel == 10 // default value
    }
}
