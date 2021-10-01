package nextflow.cloud.azure.config

import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper

/**
 * Model the settings to access to an Azure File Share.
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AzFileShareOpts implements CacheFunnel {

	static public final String DEFAULT_MOUNT_OPTIONS = '-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp'

	String mountPath
	String mountOptions

	AzFileShareOpts(Map opts) {
		assert opts!=null
		this.mountPath = opts.mountPath ?: ''
		this.mountOptions = opts.mountOptions ?: DEFAULT_MOUNT_OPTIONS
	}

	AzFileShareOpts() {
		this(Collections.emptyMap())
	}

	@Override
	Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
		hasher.putUnencodedChars(mountPath)
		hasher.putUnencodedChars(mountOptions)
		return hasher
	}
}