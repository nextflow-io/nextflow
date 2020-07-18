

import nextflow.file.FileHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestFs {


    static void main(String...ars) {

        def uri = "azb://my-data:/nf-DE8YJPrG.txt?account=nfstore"
        def file = FileHelper.asPath(uri)
        println file.text
    }
}
