# Testing

## S3 integration tests

S3 tests run only when `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` are set; otherwise they are silently skipped.

- `make test` exports them from the `nextflow-oss` SSO profile (override with `aws_profile=<name>`) and prints which case applies.
- With `./gradlew` directly, check the session with `aws sts get-caller-identity --profile nextflow-oss`, then run
  `eval "$(aws configure export-credentials --profile nextflow-oss --format env)"; unset AWS_PROFILE` in the same shell command.
- If the session is expired, ask the user to run `aws sso login --profile nextflow-oss`.
- Changing only env vars does not invalidate Gradle's test cache: pass `--no-build-cache`
  (e.g. `make test gradle="./gradlew --console=plain --no-build-cache"`).

## Other suites

- `tests/`: integration workflows with expected outputs, run with `tests/checks/run.sh`.
- `validation/`: cloud end-to-end tests; require provider credentials.
- `docs/snippets/`: verified documentation examples; update them alongside the docs that include them.
