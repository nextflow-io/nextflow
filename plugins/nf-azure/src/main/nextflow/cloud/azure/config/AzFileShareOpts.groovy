/*
 * Copyright 2021, Microsoft Corp
 * Copyright 2022, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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
		assert opts != null
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
