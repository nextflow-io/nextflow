package io.seqera.tower.plugin.exchange

import java.time.Instant

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Models a REST response containing a license-scoped JWT token from Platform
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
@CompileStatic
@ToString(includeNames = true, includePackage = false)
class GetLicenseTokenResponse {
    /**
     * The signed JWT token
     */
    String signedToken

    /**
     * The expiration timestamp of the token
     */
    Instant expiresAt

    /**
     * The Exception returned while trying to access the token
     */
    Throwable error
}
