#!/bin/sh

cleanup_at_exit()
{
    # rewind git repositories (both local and in the github forge)
    if test ! -z "$start_branch" && test "$(git branch --show-current)" != "$start_branch" ; then
        # we need to clean up the branches and rewind everything
        git reset -q --hard
        git checkout $start_branch
        delete_local_and_remote_branch "$github_token" "$rc_branch"
        delete_local_and_remote_tag    "$github_token" "$rc_tag"
        if test "$delete_release_branch_on_cleanup" = "true" ; then
            delete_local_and_remote_branch "$github_token" "$release_branch"
        fi
        git fetch --prune ${github_origin}
        echo "everything has been cleaned, branch set back to $start_branch" 1>&2
    fi
}

# run everything from top project directory
cd $(dirname $0)/..

# load common functions
. scripts/functions.sh

# safety checks
safety_checks hipparchus-core/src/main/java/org/hipparchus/analysis

# set the GitHub remote repository name
github_origin=$(find_github_origin)

# get users credentials
github_token=$(enter_private_credentials "enter your github private token")

start_branch=$(find_current_branch)
start_sha=$(find_current_sha)

# extract version numbers
pom_version=$(find_pom_version hipparchus-parent/pom.xml)
changes_version=$(find_changes_version)
release_version=$(echo $pom_version | sed 's,-SNAPSHOT,,')
release_type=$(release_type $release_version)
test "$release_version"  = "$changes_version" || complain "wrong version in changes.xml ($changes_version instead of $release_version)"
echo "current version is $pom_version"
echo "release version will be $release_version, a $release_type release"
request_confirmation "do you agree with these version numbers?"

# compute release candidate number
next_rc=$(expr $(last_rc $release_version) + 1)
rc_tag=$(rc_tag $release_version $next_rc)

# reuse existing release branch for patch release or new release candidate, create it otherwise
release_branch=$(release_branch $release_version)
if test "$release_type" = patch -o $next_rc -gt 1 ; then
    test ! -z $(git rev-parse --quiet --verify $release_branch) || complain "branch $release_branch doesn't exist, stopping"
    delete_release_branch_on_cleanup="false"
else
    git branch $release_branch || complain "branch $release_branch already exist, stopping"
    delete_release_branch_on_cleanup="true"
fi
git checkout $release_branch

# create release candidate branch
rc_branch="RC${next_rc}-${release_version}"
create_local_branch $rc_branch $start_branch

# modify pom
echo
echo "dropping -SNAPSHOT from version number"
for module in $(modules); do
    cp -p $(module_dir ${module})/pom.xml $tmpdir/original-${module}-pom.xml
    xsltproc --stringparam release-version "$release_version" \
             scripts/update-pom-version.xsl \
             $tmpdir/original-${module}-pom.xml \
             > $(module_dir ${module})/pom.xml
done
echo
git diff
echo
request_confirmation "commit modules pom.xml?"
for module in $(modules); do
    git add $(module_dir $module)/pom.xml
    git commit -m "Dropped -SNAPSHOT in version number for official release."
done

# compute release date in the future
release_date=$(TZ=UTC date -d "+3 days" +"%Y-%m-%d")

# update changes.xml
for module in $(modules); do
    if test -f $(module_dir $module)/src/changes/changes.xml ; then
        cp -p $(module_dir $module)/src/changes/changes.xml $tmpdir/original-${module}-changes.xml
        xsltproc --stringparam release-version     "$release_version" \
                 --stringparam release-date        "$release_date" \
                 --stringparam release-description "$release_version is a $release_type release." \
                 scripts/update-changes.xsl \
                 $tmpdir/original-${module}-changes.xml \
                 > $(module_dir $module)/src/changes/changes.xml
    fi
    echo
done
git diff
echo
request_confirmation "commit modules changes.xml?"
for module in $(modules); do
    test -f $(module_dir $module)/src/changes/changes.xml && git add $(module_dir $module)/src/changes/changes.xml
done
git commit -m "Updated changes.xml for official release."

# update downloads page
first_table=$(sed -n -e '/^---$/=' -e '/^---$/q' src/site/markdown/downloads.md)
sed -i -e "${first_table} a\\
\\
|     package       | link                                                                   | hash                                              | signature                                         |\\
|-------------------|------------------------------------------------------------------------|---------------------------------------------------|---------------------------------------------------|\\
|  source (zip)     | [hipparchus-${release_version}-src.zip](downloads/hipparchus-${release_version}-src.zip)             | [md5](downloads/hipparchus-${release_version}-src.zip.md5)       | [sig](downloads/hipparchus-${release_version}-src.zip.asc)       |\\
|  source (tar.bz2) | [hipparchus-${release_version}-src.tar.bz2](downloads/hipparchus-${release_version}-src.tar.bz2)     | [md5](downloads/hipparchus-${release_version}-src.tar.bz2.md5)   | [sig](downloads/hipparchus-${release_version}-src.tar.bz2.asc)   |\\
|  binary (zip)     | [hipparchus-${release_version}-bin.zip](downloads/hipparchus-${release_version}-bin.zip)             | [md5](downloads/hipparchus-${release_version}-bin.zip.md5)       | [sig](downloads/hipparchus-${release_version}-bin.zip.asc)       |\\
|  binary (tar.bz2) | [hipparchus-${release_version}-bin.tar.bz2](downloads/hipparchus-${release_version}-bin.tar.bz2)     | [md5](downloads/hipparchus-${release_version}-bin.tar.bz2.md5)   | [sig](downloads/hipparchus-${release_version}-bin.tar.bz2.asc)   |\\
version ${release_version} downloads (release date: ${release_date})\\
\\
---" src/site/markdown/downloads.md

# update site structure menus
for module in $(modules); do
    last_api_doc=$(sed -n -e "/API docs.*href=/=" $(module_dir $module)/src/site/site.xml | tail -1)
    sed -i -e "${last_api_doc} a\\
      <item name=\"${release_version} API docs\"                  href=\"https://hipparchus.org/apidocs-${release_version}/index.html\"/>" \
      $(module_dir $module)/src/site/site.xml
done

# update main page
release_out_line=$(sed -n -e "/## Release is out!/=" src/site/markdown/index.md | tail -1)
sed -i -e "$release_out_line,\\$d" src/site/markdown/index.md
cat <<EOT >> src/site/markdown/index.md
\#\# ${release_version} Release is out!

Hipparchus ${release_version} is now available for download from the
[Hipparchus download page](downloads.html) or on Maven central and its
mirrors under the groupId **org.hipparchus**.

$(extract_release_description)
EOT

# update README
h_prop="myprojectname.hipparchus.version"
sed -i -e "s,<${h_prop}>.*</${h_prop}>,<${h_prop}>${release_version}</${h_prop}>," \
    README.md

echo
git diff
echo
request_confirmation "commit documentation updates?"
git add src/site/markdown/downloads.md \
        src/site/markdown/index.md \
        README.md
for module in $(modules); do
  git add $(module_dir $module)/src/site/site.xml
done
git commit -m "Updated documentation for official release."

# push to origin
echo
test -z "$(github_get "$github_token" "git/ref/heads/$rc_branch" | jp --raw-output .object.sha)" \
  || complain "branch $rc_branch already exists in ${github_origin}"
request_confirmation "push $rc_branch branch to ${github_origin}?"
git push $github_origin $rc_branch

# generate site (this takes *a lot* of time)
LANG=C mvn -Prelease site site:stage

# create MD5 hash
for f in target/hipparchus-*.zip target/hipparchus-*.tar.bz2 ; do
  md5sum < $f > $f.md5
done

# create site archive
(cd target/staging ; tar cjf ../www-${release_version}.hipparchus.org.tar.bz2 .)

# ask release manager to upload files
request_confirmation "upload non-maven artifacts to Hipparchus site?"
scp -p \
  target/hipparchus-${release_version}-bin.zip
  target/hipparchus-${release_version}-bin.zip.asc
  target/hipparchus-${release_version}-bin.zip.md5`
  target/hipparchus-${release_version}-bin.tar.bz2
  target/hipparchus-${release_version}-bin.tar.bz2.asc
  target/hipparchus-${release_version}-bin.tar.bz2.md5
  target/hipparchus-${release_version}-src.zip
  target/hipparchus-${release_version}-src.zip.asc
  target/hipparchus-${release_version}-src.zip.md5
  target/hipparchus-${release_version}-src.tar.bz2
  target/hipparchus-${release_version}-src.tar.bz2.asc
  target/hipparchus-${release_version}-src.tar.bz2.md5
  target/www-${release_version}.hipparchus.org.tar.bz2
  ${hipparchus_server_fqdn}:~
echo "please prepare staging site on Hipparchus server so people can review the uploaded artifacts"

# make sure we can merge in a release branch on the origin server
create_remote_branch "$github_token" $release_branch $start_sha $github_origin

# remotely merge release candidate branch into release branch
merge_remote_branch "$github_token" $rc_branch $release_branch "preparing release $release_version"

# switch to release branch and delete release candidate branch
git fetch --prune ${github_origin}
git checkout $release_branch
git branch --set-upstream-to ${github_origin}/$release_branch $release_branch
git pull ${github_origin}
git branch -d $rc_branch

echo ""
request_confirmation "create tag $rc_tag?"
git tag $rc_tag -m "Release Candidate $next_rc for version $release_version."
git push $github_origin $rc_tag
echo ""

# monitor continuous integration pipeline triggering (10 minutes max)
merge_sha=$(git rev-parse --verify HEAD)
monitor_pipeline "$github_token" $merge_sha $release_branch || complain "pipeline did not succeed"

if test "$release_type" = "patch" ; then
    # for patch release, there are no votes, we jump directly to the publish step
    sh scripts/successful-vote.sh
else
    # create vote topic
    vote_date=$(TZ=UTC date -d "+3 days" +"%Y-%m-%dT%H:%M:%SZ")
    topic_title="[VOTE] Releasing Hipparchus $release_version from release candidate $next_rc"
    release_description="$(extract_release_description)"
    topic_raw="This is a vote in order to release version $release_version of the Hipparchus library.
Version $release_version is a $release_type release.

$release_description

The release candidate $next_rc can be found on the GitHub repository
as tag $rc_tag in the $release_branch branch:
https://github.com/${github_owner}/${github_project}/releases/tag/${release_version}

The maven artifacts are available in the Orekit Nexus repository at:
https://packages.orekit.org/#browse/browse:maven-releases:org%2Fhipparchus

The vote will be tallied on $vote_date (UTC time)"

    echo ""
    echo "proposed vote topic for the forum:"
    echo ""
    echo "$topic_raw"
    request_confirmation "OK to post vote topic on forum?"

    read -p "enter your forum user name " forum_username
    forum_api_key=$(enter_private_credentials "enter your forum API key")

    development_category=$(find_forum_category $forum_username $forum_api_key "Hipparchus development")
    post_url=$(post_to_forum $forum_username $forum_api_key $development_category "$topic_title" "$topic_raw")
    echo ""
    echo "vote topic posted at URL: https://forum.orekit.org/$post_url"

    echo ""
    echo "please ping PMC members so they are aware of the vote."
    echo "Their vote is essential for a release as per project governance."
fi
