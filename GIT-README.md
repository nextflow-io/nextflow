# GIT README 

Shortcuts to recurrent Git commands 

## Subtree  

The `tests` directory is a Git subtree created with the 
following commands: 

    git remote add tests git@github.com:nextflow-io/tests.git
    git subtree add --squash --prefix=tests/ tests integration


To pull changes from the [tests repo](https://github.com/nextflow-io/tests) use this command: 

    git subtree pull --squash --prefix=tests/ tests integration

To push changes to to [tests repo](https://github.com/nextflow-io/tests) use this command: 

    git subtree push --prefix=tests/ tests integration


Read more [here](https://andrey.nering.com.br/2016/git-submodules-vs-subtrees/).