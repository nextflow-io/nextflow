package nextflow.hello.functions
/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class HelloFunctions{

    String sayHello(String lang='en'){
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
