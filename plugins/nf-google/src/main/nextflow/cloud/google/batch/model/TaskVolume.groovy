/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Model a Google Batch task volume mount
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class TaskVolume {

    // A Google Cloud Storage source for the volume.
    Map gcs

    // An NFS source for the volume (could be a Filestore, for example).
    Map nfs

    // A persistent disk source for the volume.
    Map pd

    // Mount path for the volume, e.g. /mnt/share
    String mountPath

    // Mount options
    // For Google Cloud Storage, mount options are the global options supported by
    // gcsfuse tool. Batch will use them to mount the volume with the following
    // command:
    // "gcsfuse [global options] bucket mountpoint".
    // For PD, NFS, mount options are these supported by /etc/fstab. Batch will
    // use Fstab to mount such volumes.
    // https://help.ubuntu.com/community/Fstab
    List<String> mountOptions

    TaskVolume withMountOptions(List<String> opts) {
        mountOptions = opts
        return this
    }

    TaskVolume withMountPath(String path) {
        this.mountPath = path
        return this
    }
}
