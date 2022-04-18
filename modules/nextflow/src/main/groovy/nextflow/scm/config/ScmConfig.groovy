package nextflow.scm.config


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
interface ScmConfig {

    String getName()

    Map enrichConfiguration(Map map)

}