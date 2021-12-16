/*
 * MIT License
 *
 * Copyright (c) 2020 Jacob Glickman
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package com.upplication.s3fs.ng;

import java.nio.Buffer;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.NavigableMap;
import java.util.Objects;
import java.util.Optional;
import java.util.TreeMap;

/**
 * A pool of {@link T buffers}.
 * <br><br>
 * All buffers dispatched from this pool will be reused, resulting in significant performance improvements from not
 * having to constantly allocate new buffers.
 *
 * @author Jacob G.
 * @param <T> The type of {@link Buffer} that this pool holds.
 * @since February 23, 2019
 */
public abstract class AbstractBufferPool<T extends Buffer> implements AutoCloseable {

    /**
     * A heuristic used for the initial capacity of the backing {@link Deque deques}.
     */
    private static final int HEURISTIC = 3;

    /**
     * The data structure that holds all pooled {@link T buffers}.
     */
    private final NavigableMap<Integer, Deque<T>> buffers = new TreeMap<>();

    /**
     * An {@code abstract} method that allocates a new {@link T buffer} with the specified capacity.
     *
     * @param capacity the capacity of the buffer to create.
     * @return         a newly-allocated buffer.
     */
    protected abstract T allocate(int capacity);

    /**
     * Attempts to take a {@link T buffer} from this {@link AbstractBufferPool<T> pool}.
     * <br><br>
     * If no buffer can be found (with a capacity of at-least {@code capacity}) within this pool, then a new one is
     * created.
     *
     * @param capacity the capacity of the buffer requested, which will be interpreted differently for each
     *                 implementation of this pool. A {@link com.github.pbbl.heap.ByteBufferPool} measures
     *                 {@code capacity} in bytes, a {@link com.github.pbbl.heap.CharBufferPool} measures
     *                 {@code capacity} in chars, etc.
     * @return         a buffer with a capacity greater than or equal to {@code capacity}, a limit set to
     *                 {@code capacity}, and position set to {@code 0}.
     */
    public T take(int capacity) {
        Optional<T> maybeBuffer;

        synchronized (buffers) {
            maybeBuffer = buffers.tailMap(capacity, true)
                    .values()
                    .stream()
                    .map(Deque::poll)
                    .filter(Objects::nonNull)
                    .findAny();
        }

        return maybeBuffer.map((T buffer) -> {
            buffer.clear().limit(capacity);
            return buffer;
        }).orElseGet(() -> allocate(capacity));
    }

    /**
     * Gives the specified {@link T buffer} to this {@link AbstractBufferPool<T> pool}.
     * <br><br>
     * This method should only be invoked with an argument that was returned from a call to {@link #take(int)}.
     *
     * @param buffer the buffer to return to this pool.
     */
    public void give(T buffer) {
        synchronized (buffers) {
            buffers.computeIfAbsent(buffer.capacity(), capacity -> new ArrayDeque<>(HEURISTIC)).offer(buffer);
        }
    }

    /**
     * Closes this {@link AbstractBufferPool} by clearing {@link #buffers}.
     */
    @Override
    public void close() {
        synchronized (buffers) {
            buffers.clear();
        }
    }
}

