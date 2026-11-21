#!/bin/sh

cleanup_at_exit()
{
    if test ! -z "$start_branch" ; then
        git reset -q --hard
        git checkout $start_branch
        delete_local_and_remote_tag    "$github_token" "$release_tag"
        delete_local_and_remote_branch "$github_token" "$preparation_branch"
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

# set the github remote repository name
github_origin=$(find_github_origin)

# get users credentials
github_token=$(enter_private_credentials "enter your github private token")

# check branch
start_branch=$(find_current_branch)
release_version=$(find_pom_version)
release_type=$(release_type $release_version)
release_branch=$(release_branch $release_version)
test "$start_branch" = "$release_branch" || complain "$start_branch is not a release branch"

# check tag for successful release candidate
last_rc=$(last_rc $release_version)
test $last_rc -ge 1 || complain "there was no release candidate for $release_version"
rc_tag=$(rc_tag $release_version $last_rc)
release_tag="$release_version"
test -z "$(git tag -l \"${release_tag}\")" || complain "tag $release_tag already exists"

# retrieve artifacts from Orekit Nexus repository
bundle_dir=$tmpdir/bundle/org/hipparchus;
mkdir -p ${bundle_dir}/hipparchus/${release_version}
base_url="https://packages.orekit.org/repository/maven-releases/org/hipparchus"
for module in $(modules | sed 's,hipparchus-parent,,'); do
    mkdir -p ${bundle_dir}/$(module_dir ${module})/${release_version}
    for artifact in "-cyclonedx.json" "-cyclonedx.xml" ".pom" ; do
        for suffix in "" ".asc" ".md5" ".sha1" ; do
           artifact_with_suffix=$(module_dir ${module})/${release_version}/${module}-${release_version}${artifact}${suffix}
           curl --silent --output ${bundle_dir}/${artifact_with_suffix} $base_url/${artifact_with_suffix}
        done
    done
    if "${module}" != "hipparchus" ; then
        for artifact in "-javadoc.jar" "-sources.jar" ".jar" ; do
            for suffix in "" ".asc" ".md5" ".sha1" ; do
               artifact_with_suffix=$(module_dir ${module})/${release_version}/${module}-${release_version}${artifact}${suffix}
               curl --silent --output ${bundle_dir}/${artifact_with_suffix} $base_url/${artifact_with_suffix}
            done
        done
    fi
done
(cd $tmpdir/bundle ; tar czf hipparchus-${release_version}.tar.gz org)

request_confirmation "upload maven artifacts to central portal?"

# get user credentials
read -p "enter your central portal user token name" central_username
central_password=$(enter_private_credentials "enter your central portal user token password")
central_bearer=$(echo ${central_username}:$central_password | base64)

# upload maven artifacts to central portal
deployment_id=$(cd $tmpdir/bundle ;
                curl --silent \
                     --request POST \
                     --header "Authorization: Bearer $central_bearer" \
                     --form bundle=@hipparchus-${release_version}.tar.gz \
                      "https://central.sonatype.com/api/v1/publisher/upload?name=hipparchus-${release_version}&publishingType=USER_MANAGED")

echo "deployment_id: $deployment_id"

# wait for validation
status=""
timeout=0
while test "$status" != "\"VALIDATED\"" ; do
  current_date=$(date +"%Y-%m-%dT%H:%M:%S")
  echo "$current_date deployment status is \"$status\", waiting for validation…"
  sleep 5
  timeout=$(expr $timeout + 5)
  test $timeout -lt 600 || complain "deployment not validated after 10 minutes, exiting"
  status=$(curl --silent \
                --request POST \
                --header "Authorization: Bearer $central_bearer" \
                "https://central.sonatype.com/api/v1/publisher/status?id=${deployment_id}" \
           | jq .deploymentState)
done
echo "deployment has been validated"
echo ""

# publish maven artifacts to central portal
request_confirmation "publish maven artifacts to central portal?"
curl --silent \
     --request POST \
     --header "Authorization: Bearer $central_bearer" \
     "https://central.sonatype.com/api/v1/publisher/deployment/${deployment_id}"

request_confirmation "create tag $release_tag?"
git tag $release_tag -m "Version $release_version."

# push release tag to origin
request_confirmation "push $release_tag tag to ${github_origin}?"
git push $github_origin $release_tag

request_confirmation "merge $release_branch to main?"

# merge release branch to main branch
merge_remote_branch "$github_token" $release_branch main false "Release $release_version is now the main official release."
git checkout main
git pull $github_origin main
echo ""

# merge main branch to develop branch
request_confirmation "merge main to develop?"
merge_remote_branch "$github_token" main develop false "Merging main branch back to develop."
git checkout develop
git pull $github_origin develop
echo ""

release_description="$(extract_release_description)"

# prepare next development cycle
next_version=$(next_version $release_version)
preparation_branch="prepare-$next_version"
request_confirmation "prepare branch $preparation_branch for next version?"
create_local_branch $preparation_branch

for module in $(modules); do
    cp -p $(module_dir ${module})/pom.xml $tmpdir/original-${module}-pom.xml
    xsltproc --stringparam release-version "${next_version}-SNAPSHOT" \
             scripts/update-pom-version.xsl \
             $tmpdir/original-${module}-pom.xml \
             > $(module_dir ${module})/pom.xml
    if test -f $(module_dir ${module})/src/changes/changes.xml ; then
        cp -p $(module_dir ${module})/src/changes/changes.xml $tmpdir/original-${module}-changes.xml
        xsltproc --stringparam release-version "$release_version" \
                 --stringparam next-version    "$next_version" \
                 scripts/prepare-next-release-in-changes.xsl \
                 $tmpdir/original-${module}-changes.xml | \
                 sed 's,<release \([^/]*\)/>,<release \1>\n    </release>\n    ,' \
                 > $(module_dir ${module})/src/changes/changes.xml
    fi
done

echo
git diff
echo
request_confirmation "commit all modules pom.xml and changes.xml and push $preparation_branch to ${github_origin}?"
for module in $(modules); do
    git add $(module_dir ${module})/pom.xml
    test -f $(module_dir ${module})/src/changes/changes.xml && git add $(module_dir ${module})/src/changes/changes.xml
fi
git commit -m "Incremented version number for next-release."
git push $github_origin $preparation_branch

# merge to develop branch
request_confirmation "merge $preparation_branch branch to develop branch at ${github_origin}?"
merge_remote_branch "$github_token" $preparation_branch develop true "preparing development for version $next_version."
git checkout develop
git pull --prune $github_origin
git branch -d $preparation_branch
echo ""

# manage milestones
request_confirmation "close milestone $release_version?"
release_milestone_id=$(get_milestone_id "$github_token" $release_version)
close_milestone "$github_token" $release_milestone_id
if test -z "$(get_milestone_id "$github_token" $next_version)" ; then
    create_milestone "$github_token" $next_version
fi

# create GitHub release
request_confirmation "create release on GitHub repository?"
github_release_url=$(curl --silent \
                          --request POST \
                          --header "Accept: application/vnd.github+json" \
                          --header "Authorization: Bearer $github_token" \
                          --header "X-GitHub-Api-Version: 2022-11-28" \
                          --data   "{ \"tag_name\": \"${release_tag}\", \"name\": \"$release_tag\", \"body\":\"$release_description\"}" \
                         ${github_api}/releases \
                     | jq --raw-output .html_url)
echo "GitHub release URL: $github_release_url"
echo ""

# create announcement topic
topic_title="Hipparchus $release_version released"
topic_raw="The Hipparchus team is pleased to announce the release of Hipparchus version $release_version.
Version $release_version is a $release_type release.

$release_description

The maven artifacts are available in maven central

The sources and binaries can be retrieved from the forge releases page:
https://github.com/${github_owner}/${github_project}/releases"

echo ""
echo "proposed announcement topic for the forum:"
echo ""
echo "$topic_raw"
request_confirmation "OK to post announcement topic on forum?"

read -p "enter your forum user name " forum_username
forum_api_key=$(enter_private_credentials "enter your forum API key")

announcements_category=$(find_forum_category $forum_username $forum_api_key "Hipparchus announcements")
post_url=$(post_to_forum $forum_username $forum_api_key $announcements_category "$topic_title" "$topic_raw")
echo ""
echo "vote topic posted at URL: https://forum.orekit.org/$post_url"

echo "please navigate to https://github.com/${github_owner}/${github_project}/releases"
echo "and check release notes are OK (it may take some time for the releases to show up there)"
