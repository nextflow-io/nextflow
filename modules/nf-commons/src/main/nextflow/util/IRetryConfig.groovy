package nextflow.util

interface IRetryConfig {

    Duration getDelay();

    Duration getMaxDelay();

    int getMaxAttempts();

    double getJitter();

}