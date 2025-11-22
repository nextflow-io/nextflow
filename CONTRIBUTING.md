# Contributing to Nextflow

This guide documents the various ways to contribute to Nextflow, including what is required before submitting a code change.

Contributing to Nextflow doesn't just mean writing code. Helping new users in the community, testing releases and bug fixes, and improving documentation are all essential and valuable contributions. Helping in these ways is an excellent way to become an effective contributor and gain credibility within the community, which makes it easier to make larger contributions like code changes and new features.

## Helping Other Users

A great way to contribute to Nextflow is to answer user questions on the [community forum](https://community.seqera.io) and the [Nextflow Slack](https://www.nextflow.io/slack-invite.html). Contributors should ideally be active members here and keep up with the latest developments in the Nextflow community. There are always many new Nextflow users, so taking a few minutes to help answer a question is a valuable community service and a great way to demonstrate your expertise.

## Documentation Changes

Propose changes to the [Nextflow documentation](https://nextflow.io/docs/latest/) by editing the source files in the [docs](https://github.com/nextflow-io/nextflow/tree/master/docs) directory. The `README.md` in that directory describes how to build and preview the docs locally. Finally, open a pull request with the proposed changes.

## Bug Reports

Submitting a bug report is one of the simplest and most useful ways to contribute, as it helps us to quickly identify and fix issues and thereby make Nextflow more stable.

Report a bug using the **New issue** button on the [issues page](https://github.com/nextflow-io/nextflow/issues). A good bug report should include a minimal test case that can replicate the reported bug. Please follow the instructions in the issue template when submitting a bug report.

## Bug Fixes

Contributing bug fixes is the best way to gain experience with the Nextflow codebase and credibility within the community as a project contributor.

If you are new to the Nextflow codebase and want to get involved, check out issues marked as [`help wanted`](https://github.com/nextflow-io/nextflow/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22) or [`good first issue`](https://github.com/nextflow-io/nextflow/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22+). Feel free to ask for help if you get stuck while trying to implement a fix!

## New Features

Before contributing a new feature, please submit a new feature proposal on the [issues page](https://github.com/nextflow-io/nextflow/issues) and discuss it with the community.

Submitting a proposal helps identify possible overlaps with other planned features and avoid potential misunderstandings, conflicts, and wasted effort.

## Code Changes

When submitting a contribution, you will be required to sign a [Developer Certificate of Origin (DCO)](https://developercertificate.org/) to certify that you are the author of the source code or otherwise have the right to submit it to the project.

Contributor signatures are provided by adding a `Signed-off-by` line to the commit message as shown below, or by using the `-s` option with [`git commit`](https://help.github.com/articles/signing-commits/). For example:

```
This is my commit message

Signed-off-by: Random J Developer <random@developer.example.org>
```

The process is automatically managed by the [Probot](https://probot.github.io/apps/dco/) app for GitHub.

For more information about working on the Nextflow source code, visit the [Nextflow docs](https://nextflow.io/docs/latest/developer/).
