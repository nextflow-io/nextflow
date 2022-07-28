package nextflow.hello.functions

import nextflow.extension.FunctionExtensionPoint


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class HelloFunctions extends FunctionExtensionPoint{

    def sayHello(String lang='en'){
        switch( lang ){
            case 'es':
                return 'hola'
            case 'en':
                return 'hi'
            case 'it':
                return 'ciao'
            default:
                return '???'
        }
    }

}
