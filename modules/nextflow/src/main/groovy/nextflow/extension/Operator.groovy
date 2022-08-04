package nextflow.extension

/**
 * An annotation interface for operators that the plugin want to expose
 * Nextflow will search for all methods annotated with @Operators in the ExtensionPoint and allow to the user imported them
 *
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
import java.lang.annotation.ElementType
import java.lang.annotation.Retention
import java.lang.annotation.RetentionPolicy
import java.lang.annotation.Target

@Retention(RetentionPolicy.RUNTIME)
@Target([ElementType.METHOD])
@interface Operator {

}