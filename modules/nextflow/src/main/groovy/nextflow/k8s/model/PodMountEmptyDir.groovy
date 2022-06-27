/**
 * Model a K8s pod emptyDir mount
 *
 * See also https://kubernetes.io/docs/concepts/storage/volumes/#emptydir
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class PodMountEmptyDir {

    String mountPath

    Map emptyDir

    PodMountEmptyDir( Map emptyDir, String mountPath ) {
        assert emptyDir
        assert mountPath

        this.emptyDir = emptyDir
        this.mountPath = mountPath
    }

    PodMountEmptyDir( Map entry ) {
        this(entry.emptyDir as Map, entry.mountPath as String)
    }

}