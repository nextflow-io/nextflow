package io.seqera.tower.plugin.exchange


import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Models a REST request to obtain a license-scoped JWT token from Platform
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
@EqualsAndHashCode
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class GetLicenseTokenRequest {

    /** The product code */
    String product

    /** The product version */
    String version

    /**
     * The Platform workflow ID associated with this request
     */
    String workflowId

    /**
     * The Platform workspace ID associated with this request
     */
    String workspaceId

    /**
     * @return a Map representation of the request
     */
    Map<String, String> toMap() {
        final map = new HashMap<String, String>()
        map.product = this.product
        map.version = this.version
        map.workflowId = this.workflowId
        map.workspaceId = this.workspaceId
        return map
    }
}
