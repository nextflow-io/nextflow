package nextflow.scm

import nextflow.AppConst

class ScmConst extends AppConst{

    static public final String DEFAULT_BRANCH = 'master'

    static public final String DEFAULT_ORGANIZATION = System.getenv('NXF_ORG') ?: 'nextflow-io'

    static public final String DEFAULT_HUB = System.getenv('NXF_HUB') ?: 'github'

    static public final File DEFAULT_ROOT = System.getenv('NXF_ASSETS') ? new File(System.getenv('NXF_ASSETS')) : APP_HOME_DIR.resolve('assets').toFile()

}
