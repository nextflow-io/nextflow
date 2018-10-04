package nextflow.ast;

import java.nio.file.Path;

import groovy.transform.CompileStatic;
import groovy.transform.PackageScope;
import nextflow.processor.TaskPath;
import nextflow.util.Duration;
import nextflow.util.MemoryUnit;
import org.codehaus.groovy.runtime.ScriptBytecodeAdapter;
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation;

/**
 * Implements Nextflow DSL helper methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@PackageScope
@CompileStatic
class LangHelpers {

    /**
     * This method is used in place of `==` operator to allow
     * the comparisons of `Path` objects which by default are not supported
     * because it implements the Comparator interface
     *
     * See
     *  {@link ScriptBytecodeAdapter#compareEqual(Object,Object)}
     *  https://stackoverflow.com/questions/28355773/in-groovy-why-does-the-behaviour-of-change-for-interfaces-extending-compar#comment45123447_28387391
     *
     * @param left Left equals operand
     * @param right Right equals operand
     * @return 
     */
    @PackageScope
    static boolean compareEqual( Object left, Object right )  {
        // -- legacy
        if (left==right) return true;
        Class<?> leftClass = left==null?null:left.getClass();
        Class<?> rightClass = right==null?null:right.getClass();
        if (leftClass ==Integer.class && rightClass==Integer.class) {
            return left.equals(right);
        }
        if (leftClass ==Double.class && rightClass==Double.class) {
            return left.equals(right);
        }
        if (leftClass ==Long.class && rightClass==Long.class) {
            return left.equals(right);
        }

        // -- compare Paths
        if( left instanceof Path && right instanceof Path ) {
            return TaskPath.equals((Path)left, (Path)right);
        }
        // -- compare memory unit
        if( left instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)left, right)==0;
        }
        if( right instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)right, left)==0;
        }
        // -- compare duration
        if( left instanceof Duration) {
            return Duration.compareTo((Duration)left, right)==0;
        }
        if( right instanceof Duration ) {
            return Duration.compareTo((Duration)right, left)==0;
        }

        // -- fallback on default
        return DefaultTypeTransformation.compareEqual(left, right);
    }

    /**
     * Extends default less than comparison provided by {@link ScriptBytecodeAdapter#compareLessThan(Object, Object)}
     *
     * @param left Left operand
     * @param right Right operand
     * @return The comparison result
     */
    static boolean compareLessThan( Object left, Object right ) {
        Class<?> leftClass = left==null?null:left.getClass();
        Class<?> rightClass = right==null?null:right.getClass();
        if (leftClass ==Integer.class && rightClass==Integer.class) {
            return (Integer) left < (Integer) right;
        }
        if (leftClass ==Double.class && rightClass==Double.class) {
            return (Double) left < (Double) right;
        }
        if (leftClass ==Long.class && rightClass==Long.class) {
            return (Long) left < (Long) right;
        }

        // -- compare memory unit
        if( left instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)left, right)<0;
        }
        if( right instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)right, left)>0;
        }
        // -- compare duration
        if( left instanceof Duration ) {
            return Duration.compareTo((Duration)left, right)<0;
        }
        if( right instanceof Duration ) {
            return Duration.compareTo((Duration)right, left)>0;
        }

        // -- fallback on default
        return ScriptBytecodeAdapter.compareTo(left, right) < 0;
    }

    /**
     * Extends default less than comparison provided by {@link ScriptBytecodeAdapter#compareLessThanEqual(Object, Object)}
     *
     * @param left Left operand
     * @param right Right operand
     * @return The comparison result
     */
    static boolean compareLessThanEqual( Object left, Object right ) {
        Class<?> leftClass = left==null?null:left.getClass();
        Class<?> rightClass = right==null?null:right.getClass();
        if (leftClass ==Integer.class && rightClass==Integer.class) {
            return (Integer) left <= (Integer) right;
        }
        if (leftClass ==Double.class && rightClass==Double.class) {
            return (Double) left <= (Double) right;
        }
        if (leftClass ==Long.class && rightClass==Long.class) {
            return (Long) left <= (Long) right;
        }

        // -- compare memory unit
        if( left instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)left, right)<=0;
        }
        if( right instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)right, left)>=0;
        }
        // -- compare duration
        if( left instanceof Duration ) {
            return Duration.compareTo((Duration)left, right)<=0;
        }
        if( right instanceof Duration ) {
            return Duration.compareTo((Duration)right, left)>=0;
        }

        // -- fallback on default
        return ScriptBytecodeAdapter.compareTo(left, right) <= 0;
    }

    /**
     * Extends default less than comparison provided by {@link ScriptBytecodeAdapter#compareGreaterThan(Object, Object)}
     *
     * @param left Left operand
     * @param right Right operand
     * @return The comparison result
     */
    static boolean compareGreaterThan( Object left, Object right ) {
        Class<?> leftClass = left==null?null:left.getClass();
        Class<?> rightClass = right==null?null:right.getClass();
        if (leftClass ==Integer.class && rightClass==Integer.class) {
            return (Integer) left > (Integer) right;
        }
        if (leftClass ==Double.class && rightClass==Double.class) {
            return (Double) left > (Double) right;
        }
        if (leftClass ==Long.class && rightClass==Long.class) {
            return (Long) left > (Long) right;
        }

        // -- compare memory unit
        if( left instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)left, right)>0;
        }
        if( right instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)right, left)<0;
        }
        // -- compare duration
        if( left instanceof Duration ) {
            return Duration.compareTo((Duration)left, right)>0;
        }
        if( right instanceof Duration ) {
            return Duration.compareTo((Duration)right, left)<0;
        }
        // -- fallback on default
        return ScriptBytecodeAdapter.compareTo(left, right) > 0;
    }

    /**
     * Extends implementation of {@link ScriptBytecodeAdapter#compareGreaterThanEqual(Object, Object)}
     *
     * @param left Left operand
     * @param right Right operand
     * @return The comparison result
     */
    static boolean compareGreaterThanEqual( Object left, Object right ) {
        Class<?> leftClass = left==null?null:left.getClass();
        Class<?> rightClass = right==null?null:right.getClass();
        if (leftClass ==Integer.class && rightClass==Integer.class) {
            return (Integer) left >= (Integer) right;
        }
        if (leftClass ==Double.class && rightClass==Double.class) {
            return (Double) left >= (Double) right;
        }
        if (leftClass ==Long.class && rightClass==Long.class) {
            return (Long) left >= (Long) right;
        }

        // -- compare memory unit
        if( left instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)left, right)>=0;
        }
        else if( right instanceof MemoryUnit ) {
            return MemoryUnit.compareTo((MemoryUnit)right, left)<=0;
        }
        // -- compare duration
        if( left instanceof Duration ) {
            return Duration.compareTo((Duration)left, right)>=0;
        }
        else if( right instanceof Duration ) {
            return Duration.compareTo((Duration)right, left)<=0;
        }

        // -- fallback on default 
        return ScriptBytecodeAdapter.compareTo(left, right) >= 0;
    }

}
