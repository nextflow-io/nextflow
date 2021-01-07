// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;

/**
 * Only a minimal Utility class to get around a shortcoming in Core's logging.
 */
final class LoggingUtility {
    public static <T extends Exception> T logError(ClientLogger logger, T e) {
        logger.error(e.getMessage());
        return e;
    }
}
