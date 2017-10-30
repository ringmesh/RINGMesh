#!/bin/bash
# function to make a commit on a branch in a Travis CI build
# be sure to avoid creating a Travis CI fork bomb
# see https://github.com/travis-ci/travis-ci/issues/1701
function travis-branch-commit() {
    local head_ref branch_ref
    head_ref=$(git rev-parse HEAD)
    if [[ $? -ne 0 || ! $head_ref ]]; then
        err "failed to get HEAD reference"
        return 1
    fi
    branch_ref=$(git rev-parse "$TRAVIS_BRANCH")
    if [[ $? -ne 0 || ! $branch_ref ]]; then
        err "failed to get $TRAVIS_BRANCH reference"
        return 1
    fi
    if [[ $head_ref != $branch_ref ]]; then
        msg "HEAD ref ($head_ref) does not match $TRAVIS_BRANCH ref ($branch_ref)"
        msg "someone may have pushed new commits before this build cloned the repo"
        return 0
    fi
    if ! git checkout "$TRAVIS_BRANCH"; then
        err "failed to checkout $TRAVIS_BRANCH"
        return 1
    fi

    if ! git add --all "$TRAVIS_BUILD_DIR"; then
        err "failed to add modified files to git index"
        return 1
    fi
    # make Travis CI skip this build
    if ! git commit -m "Travis CI update"; then
        err "failed to commit updates"
        return 1
    fi
    if ! git push https://$GITHUB_TOKEN@github.com/$TRAVIS_REPO_SLUG "$TRAVIS_BRANCH" > /dev/null 2>&1; then
        err "failed to push git changes"
        err "$GITHUB_TOKEN"
        return 1
    fi
}

function msg() {
    echo "travis-commit: $*"
}

function err() {
    msg "$*" 1>&2
}

travis-branch-commit
