package io.seqera.tower.plugin.exchange

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Models a REST response containing a license-scoped JWT token from Platform
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
@CompileStatic
@ToString(includeNames = true, includePackage = false)
class LicenseTokenResponse {
    /**
     * The signed JWT token
     */
    String signedToken

    /**
     * The expiration date of the token
     */
    Date expirationDate
}
