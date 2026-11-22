<!--
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
# Hipparchus Release Guide

This release guide lists the steps that are needed to release a new version of Hipparchus.  
When in doubt, ask a question on the ["Hipparchus development" section of the forum](https://forum.orekit.org/c/hipparchus-development/16).

Three types of versions can be released:
 - **Major version**: numbered xx.0, a version that can introduce both new features and bug corrections.  
 This is the only type of version where APIs are allowed to evolve, creating eventual incompatibilities with former versions.
 - **Minor version**: numbered xx.y, with y > 0, a version that can introduce both new features and bug corrections but where breaking APIs is forbidden.  
 A version xx.y must be perfectly compatible with its xx.0 major version counterpart.
 - **Patch version**: numbered xx.y.z, a version where only bug corrections are allowed.  
 Once again APIs incompatibilities with version xx.y or xx.0 are not allowed.  
 A patch version is the only type of version where a vote from the PMC isn't required to publish the release.

# 0. Prerequisites

The following tasks are one-time actions (or at least they should only be renewed infrequently).
They just prepare the credentials for performing the release, but these credentials will remain
valid for months and for several releases.

## Forum API Key

The release process shell script described below creates the vote topic automatically
on the Orekit forum https://forum.orekit.org. In order to do this, it needs an API key
to post on behalf of the release manager. There are no automatic processes for setting
up such a key, candidate release managers should just send direct messages to the forum
administrators, and they will proceed with generating the API key and sending it, using
any means they see fit.

## Signing key for non-maven artifacts

Release managers should set up a gpg key for signing non-maven artifacts. This can be done
using

    gpg --gen-key

This key must be provided to other Hipparchus server admins so they can add it to the
set of trusted keys file at https://www.hipparchus.org/KEYS. There are instructions at the
top of this file explaining how to export the key.

## Access to the Hipparchus server upload area

In order to upload the non-maven artifacts to the Hipparchus server, release managers
must set up a SSH key and provide its public part to the Hipparchus server admins so
they can authorize access. This can be done using a command similar to:

    ssh-keygen -t ed25519 \
               -C "Access to Hipparchus server upload area" \
               -f .ssh/id_ed25519_hipparchus_upload

## GitHub personal access token

The release process shell script described below creates releases in the Hipparchus GitHub
mirror repository. In order to do this, it needs a personal access token with repos and workflow
access. This access token can be created by users directly from their GitHub account (User Settings →
Access → Credentials → Personal Access Tokens (classic)).

Note that the access token is only visible once from the creation page, so it must be saved
before leaving the page. The GitHub personal access token have an expiration date, hence this step
must be renewed from time to time.

## Sonatype central portal

Binary artifacts are published in the maven central repository, which is managed by Sonatype.
The `org.hipparchus` namespace is used for these artifacts, and only people that have
been identified as allowed to publish releases can push binary artifacts for publication.

In order to publish, candidate release managers should therefore create an account on
Sonatype central repository using the [signup page](https://my.sonatype.com/). Once the
account has been created they should ask existing release managers to add them to the
allowed list for the `org.hipparchus` namespace.

After creating the account on Sonatype central portal, release managers should create
a User Token. This is done by login to [central portal login page](https://central.sonatype.com/),
selecting the `View User Tokens` item from the menu drop list at the upper right of the page (small
button that displays the start of the username) and selecting the `Generate User Token` button in
the `Tokens` page. This will generate a token with a small random Id and a long random password.

Note that Hipparchus does *not* use the `central-publishing-maven-plugin`. The reason is
that we need to provide public access to the binary artifacts for the voting process. As
of summer 2025, the central portal does not allow that, it restricts the visibility of
staged artifacts before they are officially published. So the Hipparchus process is to use
the Orekit Nexus repository at `https://packages.orekit.org/` during the voting period,
and to transfer the artifacts from the Orekit Nexus repository to central portal only
after the vote has been tallied and publication has been allowed.

As we don't use the `central-publishing-maven-plugin`, there is no need to copy the xml
snippet that is generated by central portal into the `$HOME/.m2/settings.xml` file. Users
could just store these credentials in any password manager they prefer, the script will
ask for them when needed.

The Sonatype [publish portal guide](https://central.sonatype.org/publish/publish-portal-guide/)
explains everything in detail, but this release guide is enough to perform a release.

## Note on maven plugins versions

Maven plugins should be updated from time to time, but it is probably
unwise to do it at release time, it is too late, so these updates
should happen well before the release. All maven plugin versions are
gathered at one place, in a set of properties in `hipparchus-parent/pom.xml`:

    <!-- Project specific plugin versions -->
    <hipparchus.spotbugs-maven-plugin.version>4.10.5.0</hipparchus.spotbugs-maven-plugin.version>
    <hipparchus.jacoco-maven-plugin.version>0.8.16</hipparchus.jacoco-maven-plugin.version>
    <hipparchus.maven-assembly-plugin.version>3.8.0</hipparchus.maven-assembly-plugin.version>
    ...

You can find the latest version of the plugins using the search feature at
[https://central.sonatype.com/](https://central.sonatype.com/). The
properties name all follow the pattern `hipparchus.some-plugin-name.version`, the
plugin name should be used in the web form to check for available versions.

Beware that in some cases, the latest version cannot be used due to
incompatibilities. For example, when a plugin was not recently updated and
conflicts appear with newer versions of its dependencies.

Beware also that some plugins use configuration files that may need update too.
This is typically the case with `maven-checkstyle-plugin` and
`spotbugs-maven-plugin`. The `src/conf/checkstyle.xml` and
`hipparchus-<some-module>/src/conf/spotbugs-exclude-filter.xml` files may need to be checked.

Before committing these changes, you have to check that everything works. So
run the following command, taking care to use the Java version that matches
the compatibility requirements (Java 21 at time of writing):

    mvn clean
    LANG=C JAVA_HOME=/path/to/correct/java/version PATH=${JAVA_HOME}/bin:$PATH mvn -Prelease site site:stage

If something goes wrong, either fix it by changing the plugin configuration or
roll back to an earlier version of the plugin.

Browse the generated site starting at page `target/staging/index.html` and check
that everything is rendered properly.

When everything runs fine and the generated site is OK, then you can commit the
changes:

    git add hipparchus-parent/pom.xml src/conf/checkstyle.xml hipparchus-*/src/conf/spotbugs-exclude-filter.xml
    git commit -m "Updated maven plugins versions."


# 1. Verify the status of start branch

Major and minor versions are released starting from the develop branch. Patch versions
are released from a dedicated branch (typically named patch-X.Y.Z where X.Y.Z is the
candidate version number).

Before anything, check on the
[continuous integration](https://sonar.orekit.org/dashboard?id=hipparchus%3Ahipparchus&branch=develop)
that everything is fine on develop branch for major or minor release; adapt the URL for patch
release which starts from dedicated branches and not from the develop branch.

Start from the
[milestone](https://github.com/Hipparchus-Math/hipparchus/milestones) of
the version, note the issues in the scope that are not closed, and
find the associated pull requests.

Then, for each pull request:
* Check that the pipeline was a success (i.e. all stages passed)
* If it exists, find the SonarQube report for the pull request, either in a short-lived branch of
  [Hipparchus main CI report](https://sonar.orekit.org/dashboard?id=hipparchus%3Ahipparchus) or on
  the developer's own [sub-project](https://sonar.orekit.org/projects?sort=-analysis_date)
  on SonarQube.  

On the start branch, check that:
* All tests pass;
* Code coverage is up to the requirements;
* There are no bugs, vulnerabilities or code smells.

If not, fix the warnings and errors first!

It is also necessary to check on the
[GitHub actions](https://github.com/Hipparchus-Math/hipparchus/actions)
that everything is fine on develop branch (i.e., all stages are passed); here again,
adapt the URL for patch release which starts from dedicated branches and not from
the develop branch.

### 1.1. Merging the remaining merge requests

For patch releases, it is possible to include merge requests if they were created from the
`release-X.Y` branch, but not if they were created from the `develop` branch as it would
include evolutions and not only bug fixes.

Note that if there aren't any conflicts in GitHub you can directly do
the merge from GitHub; just make sure that the target branch is
`release-X.Y` and not "develop".

Find the pull request on the repository, it should be in branch `origin/merge-requests/XXX`
with XXX the number of the pull request

    git merge --no-ff origin/merge-requests/XXX

The `--no-ff` option will force a merge commit to be made, instead of doing a fast-forward on the branch.  
You can keep the automatic merge message unless you want to add some content to it.

Eventually, resolve any conflict and commit the result.

### 1.2. Cherry-picking the commits

If the developers started from develop branch (or any other branch),
then the pull request branch may contain code that should not be added to the
release.  You will have to cherry-pick the appropriate commits and add
them to a dedicated branch.  It is advised to use an IDE to do the
cherry-picking, although the command lines below will help you.

Find the pull request on the repository and the commits you are interested in.

Create a dedicated branch for the issue in your local repository:

    git checkout -b issue-YYY

Where YYY is the issue number that the pull request XXX fixes.  
If the branch already exists, give it a different name like `integrate-issue-YYY` or `integrate-PR-XXX`

Make a list of the IDs of the commits you want to add, example `A B C D` (in the order they were committed).

Cherry-pick the commits in a chronological order:

    git cherry-pick A B C D

Eventually, resolve any conflict and commit the result.

Return to the start branch and merge the branch issue-YYY:

    git checkout <name of the start branch>
    git merge --no-ff issue-YYY

# 2. Merging the module-specific changes.xml files

Hipparchus is a multi-module product hence there are several module-specific
`hipparchus-<some-module>/src/changes/changes.xml` files, and one top-level
`src/changes/changes.xml` aggregator file. Copying the entries from the module-specific
files to the top-level aggregation file is not yet managed by the release scripts below
(contributions welcome). It has to be done manually beforehand. Note that the release date*
and release type (major, minor or patch) in the various `changes.xml` files will be
computed automatically by the script, so it should remain set to `TBD` when manually
editing the files.

# 3.Run the prepare-release.sh script

There is a `prepare-release.sh` shell script in the `scripts`
directory that automatically performs all the release steps by
itself. It asks questions and waits for confirmation at various
stages. It creates branches, checks them out, performs some automated
edition on several files, commits the changes, sets up tags that will
be recognized by Hipparchus Continuous Integration setup, and pushes
everything to Hipparchus GitHub repository.

When releasing a major or minor version, the script automatically
posts a vote topic on the Orekit forum so that community and PMC
members can vote on the release. When releasing a patch version, the
script skips the vote topic creation and continues up to promoting the
release from staging to published.

The shell script performs the release on a temporary branch that is
merged to the release branch after the last automatic file edition.

When the temporary branch is remotely merged to the release branch on
the Hipparchus GitHub repository, the Continuous Integration is automatically
triggered. It performs a full build, signs the artifacts using the
release key and puts them in a staging repository, where users will
later be able to review them during the voting period.

Several steps (creating the release branch, merging the temporary
branch into the release branch, creating the binary artifacts…) are
performed remotely, from the release manager computer to the Hipparchus
GitHub repository (this is the reason why release managers need a GitHub
personal access token).  Several of these steps take some time, for
example when the merge request is created, it is not executed
immediately. There are several waiting loops that monitor the progress
of these steps, polling the forge regularly. These loops display a
timestamped waiting message evey 5 seconds (every 30 seconds for the
long pipeline run loop) until the operation is performed, and stop the
release process if the operation exceeds a timeout. Branch creation,
merge request creation, and merge request completion have a 10 minutes
timeout but should typically be completed in one or two
minutes. Continuous Integration pipeline start has a 30 minutes
timeout as runners may be busy with already running compilations. Once
started, Continuous Integration pipeline run time has a 120 minutes
timeout. The final step (creating the vote thread) is triggered only
if the Continuous Integration pipeline succeeds.

Note that the binary artifacts are *not* generated locally on the
release manager computer, they are generated (and signed) by the
Continuous Integration server.

If something goes wrong, or if the release manager answers "no" to one
of the confirmation questions triggered by the script, or if the release
manager interrupts the script by a signal like hitting ctrl-C, then the
script removes all the branches and tags it created (both locally and on
the Hipparchus GitHub repository) and then stops. This ensures everything
remains clean. The last chance to stop the release process corresponds to
the question about posting the vote thread to the forum.

This script must be run from the command line on a computer with
several Linux utilities (git, sed, xsltproc, sort, tail, curl, stty,
base64, jq, md5sum and scp), with the git worktree already set to the start branch (i.e.,
develop or a patch branch):

    sh scripts/prepare-release.sh

Here are the steps the script performs on its own, asking
for user confirmation before any commit:

    - perform safety checks (files and directories present, utilities available, java version)
    - ask for GitHub personal access token
    - check if the release is a major, minor or patch release
      using the -SNAPSHOT version number from the current branch top-level `pom.xml`
    - set up a release candidate number, starting from 1 or incrementing existing RC tags
    - for first release candidate of major or minor release, create a release-X.Y branch from develop
      (reuse existing branch otherwise)
    - checkout the release-X.Y branch
    - create a RCn-X.Y branch
    - checkout the RCn-X.Y branch
    - merge the start branch into the RCn-X.Y branch
    - drop -SNAPSHOT version number from all modules `pom.xml` files and commit the change
    - compute candidate release date, allocating 3 days for the vote
    - update all `changes.xml` files with release date and release type and commit the change
    - update site structure menus and downloads page and commit the changes
    - generate the non-maven artifacts (binary, sources and site), this takes a lot of time
    - upload the non-maven artifacts to Hipparchus server
    - push the RCn-X.Y branch to origin
    - if the release-X.Y branch does not exist on the Hipparchus GitHub repository, create it remotely
    - trigger a merge request from RCn-X.Y branch to release-X.Y branch
      (this will trigger full build, signing and deployment to Orekit Nexus instance)
    - switch to release-X.Y branch
    - pull merged branch from origin
    - tag the RCn-X.Y branch with tag X.Y-RCn and push the tag to origin
    - monitor continuous integration run
    - if the release is a patch release, call immediately the `successful-vote.sh` script
    - if the release is a major or minor release, create the vote topic automatically on the forum
      (the API key for posting to the forum will be asked for at this point)

## Managing non-maven artifacts

The non-maven artifacts are the web site and the source and binary archives. They are
generated and uploaded to the Hipparchus server by the script above, but only to the
release manager home directory. Some manual operations need to be performed on this
server:
    - a new `staging/downloads` directory must be created in the web server
    - the non-maven artifacts should be moved there
    - the generated website www-X.Y.hipparchus.org.tar.bz2 must be unpacked in the `staging` directory
    - the Apache server configuration must be updated to allow this `staging` directory to be served
      (it is not served when no vote is going on)
    - the `apache2ctl restart` command must be run to relaunch the Apache server after the configuration change

These operations must be performed on the Hipparhus server itself, so they are *not* performed by the script.

## 4. Voting period

For major and minor releases, the script generates the vote topic on
the forum by itself. The voting period is 3 days, so release managers
should wait until this time has elapsed before concluding. Of course,
if one of the voters raises concerns that should be addressed, the
vote can be canceled before the 5-day duration.

### 4.1. Failed or canceled vote

If the vote fails or is canceled, the only thing to do is edit the title
of the vote thread on the forum to add a [CANCELED] or [FAILED] label.
This is done manually directly on the forum web pages. The failed
release candidate can be kept, it will just not be published.

Once the problem has been addressed and a new release candidate can be
created, the release process should just be started again using the
`prepare-release.sh` script.  Note that the script keeps track of the
release candidate number on its own as it looks for the tags in the
git repository and as release candidate tags are not deleted if the
release preparation script has completed.

### 4.2. Successful vote

When the vote for a release candidate succeeds (or automatically when
the release is a patch release), the maven artifacts must be published
to the central portal. A final tag (without reference to release
candidate number) must be created. The repository must also be prepared
for the next development cycle. This is done by running the
`successful-vote.sh` script, taking care to be in the release-X.Y
branch:

    sh scripts/successful-vote.sh

Here are the steps the script performs on its own, asking
for user confirmation before any commit:

    - perform safety checks (files and directories present, utilities available)
    - ask for GitHub personal access token
    - check the branch to see it is really a release branch
    - check the tag for the successful release candidate
    - retrieve maven artifacts from the Orekit Nexus repository
    - publish maven artifacts to the central portal
      (the API key for connecting to the central portal will be asked for at this point)
    - tag the repository
    - push the tag to origin
    - trigger a merge request from release branch to main branch
    - pull merged main branch from origin
    - trigger a merge request from main branch to develop branch
    - pull merged develop branch from origin
    - create a temporary preparation branch for next development cycle
    - update version number in all modules `pom.xml` files
    - update version number in all modules `changes.xml` files
    - push temporary preparation branch to origin
    - trigger a merge request from preparation branch to develop branch
    - pull merged develop branch from origin
    - delete temporary preparation branch
    - close milestone for released version on the GitHub repository
    - create milestone for next version on the GitHub repository if it doesn't exist yet
    - create release on the GitHub repository
    - create the annoucement topic automatically on the forum
      (the API key for posting to the forum will be asked for at this point)
