package io.seqera.tower.plugin.datalink

import groovy.transform.CompileStatic

@CompileStatic
class DataLink {
    private String name
    private String id
    private String credentialsId
    private String provider

    DataLink(String name, String id, String credentialsId, String provider) {
        this.name = name
        this.id = id
        this.credentialsId = credentialsId
        this.provider = provider
    }
    String getName(){
        return name
    }

    String getId() {
        return id
    }

    String getCredentialsId() {
        return credentialsId
    }

    String getProvider() {
        return provider
    }

}
