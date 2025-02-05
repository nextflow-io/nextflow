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
}
