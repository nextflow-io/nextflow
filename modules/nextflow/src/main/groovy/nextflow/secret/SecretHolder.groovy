/*
 * Copyright 2021, Sage-Bionetworks
 *
 * All Rights reserved
 *
 */

package nextflow.secret

/**
 * Hold a config secret
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SecretHolder extends Closure {

    private String secretName

    private SecretsProvider provider

    SecretHolder(String name) {
        super(null, null);
        setResolveStrategy(Closure.TO_SELF)
        this.secretName = name
    }

    SecretHolder bind(SecretsProvider provider) {
        this.provider = provider
        return this
    }

    boolean isBound() {
        return provider != null
    }

    String getSecretName() {
        return secretName
    }

    Object getSecretValue() {
        if( provider == null )
            throw new IllegalStateException("Missing secret provider — unable to resolve secret '$secretName'")
        final secret = provider.getSecret(secretName)
        if( secret==null )
            throw new MissingSecretException("Unknown config secret '$secretName'")
        return secret.value
    }

    Object call(final Object... args) {
        return resolve0()
    }

    @Override
    Object call() {
        return resolve0()
    }

    private Object resolve0() {
        isBound() ? getSecretValue() : "secrets.$secretName"
    }
}
