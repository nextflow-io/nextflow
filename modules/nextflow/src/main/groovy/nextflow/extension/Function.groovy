package nextflow.extension

/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
import java.lang.annotation.ElementType
import java.lang.annotation.Retention
import java.lang.annotation.RetentionPolicy
import java.lang.annotation.Target

@Retention(RetentionPolicy.RUNTIME)
@Target([ElementType.METHOD, ElementType.TYPE])
@interface Function {

}