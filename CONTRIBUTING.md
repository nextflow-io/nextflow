# CONTRIBUTING TO NEXTFLOW

This guide documents the best way to make various types of contributions to Nextflow,
including what is required before submitting a code change.

Contributing to Nextflow doesn't just mean writing code. Helping new users on the mailing list,
testing releases and bug fixes, and improving documentation are all essential and valuable contributions. In fact, proposing
significant code changes usually first requires gaining experience and credibility within the
community by helping in other ways. This is also a guide to becoming an effective contributor.


## Contributing by Helping Other Users

A great way to contribute to Nextflow is to help answer user questions on the [discussion forum](https://groups.google.com/forum/#!forum/nextflow)
or the [Gitter channel](https://gitter.im/nextflow-io/nextflow). There are always many new Nextflow users;
taking a few minutes to help answer a question is a very valuable community service.

Contributors should ideally subscribe to these channels and follow them in order to keep up to date
on what's happening in Nextflow. Answering questions is an excellent and visible way to help the
community and also demonstrates your expertise.


## Contributing Documentation Changes

To propose a change to release documentation (that is, the docs that appear under http://docs.nextflow.io),
edit the documentation source files in Nextflow's [docs/](https://github.com/nextflow-io/nextflow/tree/master/docs)
directory, whose README file shows how to build the documentation locally to test your changes.

Then open a pull request with the proposed changes.


## Contributing Bug Reports

Filling a bug report is likely the simplest and most useful way to contribute to the project.
It helps us to identify issues and provide patches and therefore to make Nextflow more stable
and useful.

Report a bug using the "New issue" button in the
[issues page](https://github.com/nextflow-io/nextflow/issues) of this project.

A good bug report should include a minimal executable test case able to replicate the reported bug.

Follow the instructions in the bug [report template](https://github.com/nextflow-io/nextflow/blob/master/.github/issue_template.md) that is shown when filling the bug report out.

## Contributing Bug Fixes

Contributing bug fixes is the best way to gain experience and credibility within the community
and also to become an effective project contributor.

If you are a novice with the Nextflow code base, start by looking at issues marked
with the [help wanted](https://github.com/nextflow-io/nextflow/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22)
label.

If you have doubts on how to fix an issue, ask for help from senior contributors commenting
in the issue page.

## Contributing New Features

Before contributing a new feature, submit a new feature proposal in the
[issues page](https://github.com/nextflow-io/nextflow/issues) of the project and discuss it
with the community.

This is important to identify possible overlaps with other planned features and avoid misunderstandings and conflicts.

## Contributing Code Changes

When submitting a contribution, you will be required to sign a [Developer Certificate of Origin (DCO)](https://developercertificate.org/) to certify that you are the author of the source code or otherwise you have the right to submit it to the project. 

Contributor signatures are provided by adding a `Signed-off-by` line to the commit message 
as shown below, or by using the `-s` option for the [git commit command](https://help.github.com/articles/signing-commits/).

```
This is my commit message

Signed-off-by: Random J Developer <random@developer.example.org>
```

The process is automatically managed by the [Probot](https://probot.github.io/apps/dco/) app for GitHub.


## IDE settings

The suggested development environment is [IntelliJ IDEA](https://www.jetbrains.com/idea/download/). See the [README](https://github.com/nextflow-io/nextflow/#intellij-idea) for a short primer on how to import
and configure Nextflow to work with it.

Nextflow does not impose a strict code formatting style, however the following setting should be applied:

* Use spaces for indentation
* Tab size: 4
* Indent: 4
* Use single class import
* Class count to use import with `*`: 99
* Names count to use static import with `*`: 99
* Imports layout:
    * \<blank line>
    * `import org.junit.*`
    * `import spock.lang.*`
    * \<blank line>
    * `import java.*`
    * `import javax.*`
    * \<blank line>
    * *all other imports*
    * *all other static imports*

New files must include the appropriate license header boilerplate and the author name(s) and contact email(s) ([see for example](https://github.com/nextflow-io/nextflow/blob/e8945e8b6fc355d3f2eec793d8f288515db2f409/modules/nextflow/src/main/groovy/nextflow/Const.groovy#L1-L15)).
