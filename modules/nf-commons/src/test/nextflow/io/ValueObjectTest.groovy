package nextflow.io

import groovy.test.GroovyAssert
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ValueObjectTest extends Specification {

    void 'should validate ValueObjectBase annotation'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.SerializableObject
            
            @SerializableObject
            class Foo {
                String foo
            }
            
            def x = new Foo(foo: 'hello')
            assert x instanceof Serializable 
            '''
    }


    void 'should validate ValueObject annotation'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.ValueObject
            
            @ValueObject
            class Foo {
                String foo
            }
            
            def x = new Foo(foo: 'hello')
            def y = new Foo(foo: 'hello')
            assert x instanceof Serializable 
            assert x == y 
            '''
    }

    void 'should be cloneable'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.ValueObject
            
            @ValueObject
            class Obj {
                String foo
                String bar
            }
            
            def obj = new Obj(foo: 'hello', bar:'world')
            assert obj == obj.clone()
            '''
    }

    void 'should be copyable'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.ValueObject
            
            @ValueObject
            class Obj {
                String foo
                String bar
            }
            
            def obj = new Obj(foo: 'hello', bar:'world')
            def that = obj.copyWith(foo:'hola')
            
            assert that.foo == 'hola'
            assert that.bar == 'world' 
            '''
    }

}
