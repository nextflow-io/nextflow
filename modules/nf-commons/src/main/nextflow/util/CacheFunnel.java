package nextflow.util;

import com.google.common.hash.Hasher;
/**
 * Interface to delegate cache hashing to
 * a the implementing object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface CacheFunnel {

    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode);

}
