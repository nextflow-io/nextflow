package nextflow.util

import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Helper class to hide sensitive data 
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SecretHelper {

    static public final Pattern SECRET_KEYS = ~/(?im)^AWS.+|.*TOKEN.*|.*SECRET.*|.*accessKey.*/

    // note: ?i stands for ignore case - ?m stands for multiline
    static public final Pattern SECRET_REGEX = ~/(?im)(^AWS[^=]*|.*TOKEN[^=]*|.*SECRET[^=]*)=(.*)$/

    static String secureEnvString( String str ) {
        str.replaceAll(SECRET_REGEX, '$1=[secure]')
    }

    static Object hideSecrets( obj ) {
        if( obj == null )
            return 
        
        if( obj instanceof Map ) {
            final names = obj.keySet()
            for( String n : names )  {
                if( SECRET_KEYS.matcher(n).find() ) {
                    obj.put(n, '[secret]')
                }
                else {
                    hideSecrets(obj.get(n))
                }
            }
        }
        else if( obj instanceof Collection ) {
            for( Object item : ((Collection)obj) ) {
                hideSecrets(item)
            }
        }
        else if( obj.getClass().isArray() ) {
            for( Object item : ((Object[])obj) ) {
                hideSecrets(item)
            }    
        }

        return obj
    }
    
}
