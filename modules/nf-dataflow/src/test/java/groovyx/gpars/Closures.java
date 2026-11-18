/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package groovyx.gpars;

import java.util.Arrays;

import groovy.lang.Closure;

/**
 * Builds Groovy closures from Java lambdas, declaring the arity that GPars inspects.
 */
public final class Closures {

    private Closures() {}

    @FunctionalInterface
    public interface Action0 { void call() throws Throwable; }

    @FunctionalInterface
    public interface Action1 { void call(Object a) throws Throwable; }

    @FunctionalInterface
    public interface Action2 { void call(Object a, Object b) throws Throwable; }

    @FunctionalInterface
    public interface Action3 { void call(Object a, Object b, Object c) throws Throwable; }

    @FunctionalInterface
    public interface Function0 { Object call() throws Throwable; }

    @FunctionalInterface
    public interface Function1 { Object call(Object a) throws Throwable; }

    /**
     * A closure body that also receives the closure delegate, e.g. the operator or the agent.
     */
    @FunctionalInterface
    public interface Body { Object call(Object delegate, Object[] args) throws Throwable; }

    public static Closure<Object> action(Action0 fn) {
        return closure(0, (d, args) -> { fn.call(); return null; });
    }

    public static Closure<Object> action(Action1 fn) {
        return closure(1, (d, args) -> { fn.call(args[0]); return null; });
    }

    public static Closure<Object> action(Action2 fn) {
        return closure(2, (d, args) -> { fn.call(args[0], args[1]); return null; });
    }

    public static Closure<Object> action(Action3 fn) {
        return closure(3, (d, args) -> { fn.call(args[0], args[1], args[2]); return null; });
    }

    public static Closure<Object> function(Function0 fn) {
        return closure(0, (d, args) -> fn.call());
    }

    public static Closure<Object> function(Function1 fn) {
        return closure(1, (d, args) -> fn.call(args[0]));
    }

    public static Closure<Object> closure(int arity, Body body) {
        return new LambdaClosure(arity, body);
    }

    private static final class LambdaClosure extends Closure<Object> {

        private final Body body;

        LambdaClosure(int arity, Body body) {
            super(null, null);
            this.body = body;
            this.maximumNumberOfParameters = arity;
            this.parameterTypes = new Class<?>[arity];
            Arrays.fill(this.parameterTypes, Object.class);
        }

        @Override
        public Object call(Object... args) {
            final Object[] actual = args != null ? args : new Object[] { null };
            try {
                return body.call(getDelegate(), actual);
            }
            catch( Throwable e ) {
                throw Closures.<RuntimeException>sneakyThrow(e);
            }
        }

        public Object doCall(Object... args) {
            return call(args);
        }
    }

    @SuppressWarnings("unchecked")
    private static <E extends Throwable> E sneakyThrow(Throwable e) throws E {
        throw (E) e;
    }
}
