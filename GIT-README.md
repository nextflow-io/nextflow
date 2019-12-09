# GIT README 

Shortcuts to recurrent Git commands 

## Sync remote fork 

```bash
git fetch public
git checkout master
git merge public/master
```

List branch remote

    git branch -vv 

Checkout from a remote branch 

    git co -b <local name> upstream/master

Push to a remote upstream branch

    git push <remote> <local branch>:<remote branch>

eg:

    git push upstream foo:master

Read more [here](https://help.github.com/articles/syncing-a-fork/).

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

## Stash shortcuts

    git stash list
    git stash pop
    git stash pop stash@{1}
    git showtool stash@{0}
    git stash drop
    git stash drop stash@{1}
    git stash clear
    git diff stash
    git diff stash@{1} [other]

## Misc 

Find a commit in any branch introducing a change

    git log -S <whatever> --source --all

## GPG keys 

To sign Git commits with a GPG key on Mac use [GPG Suite](https://gpgtools.org/), import your key, then: 

    git config --global gpg.program /usr/local/MacGPG2/bin/gpg2
    git config --global user.signingkey <your key> 
    git config --global commit.gpgsign true 
    git config --global format.signoff true ## TO AVOID TO SPECIFY -S option each time

Read more: 
https://gist.github.com/danieleggert/b029d44d4a54b328c0bac65d46ba4c65

