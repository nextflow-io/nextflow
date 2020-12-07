// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

/**
 * RESERVED FOR INTERNAL USE.
 *
 * An enum to indicate the status of a directory.
 */
enum DirectoryStatus {
    EMPTY, // The directory at least weakly exists and is empty.

    NOT_EMPTY, // The directory at least weakly exists and has one or more children.

    DOES_NOT_EXIST, // There is no resource at this path.

    NOT_A_DIRECTORY; // A resource exists at this path, but it is not a directory.

    static boolean isDirectory(DirectoryStatus status) {
        return EMPTY.equals(status) || NOT_EMPTY.equals(status);
    }
}
