/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package com.upplication.s3fs.experiment;

import java.math.BigInteger;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicReference;

import com.google.common.annotations.Beta;

/**
 * Implement an BitInteger with *atomic* behavior
 *
 * DON'T USE, IT ADDS NOT NEGLIGIBLE EXECUTION OVERHEAD
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Beta
public final class AtomicBigInteger {

    private final AtomicReference<BigInteger> current;

    public AtomicBigInteger() {
        this(new BigInteger("0"));
    }

    public AtomicBigInteger(final BigInteger value) {
        this.current = new AtomicReference<>(Objects.requireNonNull(value));
    }

    public BigInteger incrementAndGet() {
        return current.accumulateAndGet(BigInteger.ONE, BigInteger::add);
    }

    public BigInteger incrementAndGet(BigInteger value) {
        return current.accumulateAndGet(value, BigInteger::add);
    }

    public BigInteger incrementAndGet(long value) {
        return incrementAndGet(BigInteger.valueOf(value));
    }

    public BigInteger getAndIncrement() {
        return getAndIncrement(BigInteger.ONE);
    }

    public BigInteger getAndIncrement(long value) {
        return getAndIncrement(BigInteger.valueOf(value));
    }

    public BigInteger getAndIncrement(BigInteger value) {
        return current.getAndAccumulate(value, BigInteger::add);
    }

    public BigInteger get() {
        return current.get();
    }
}
