#!/bin/sh

tmpdir=$(mktemp -d /tmp/hipparchus-script.XXXXXX)
trap "rm -fr $tmpdir" 0
trap "exit 1" 1 2 15

#TODO: update owner before merging to Hipparchus official repository
#github_owner=Hipparchus-Math
github_owner=maisonobe
github_project=hipparchus
github_api=https://api.github.com/repos/${github_owner}/${github_project}
github_api_version=2026-03-10
forum_fqdn=forum.orekit.org
hipparchus_server_fqdn=ganymede.hipparchus.org

# note that the following function refers to a cleanup_at_exit function
# that is *not* defined here but must be defined in the calling script
complain()
{
    echo "$1" 1>&2
    cleanup_at_exit
    # we don't use "exit" here because sometimes we call this from command substitution like in:
    # some_variable=$(some_function_which_could_call_complain)
    # in this case, we still want to exit the entire script, not just the function
    kill -HUP $$
}

request_confirmation()
{
    answer=""
    while test "$answer" != "yes" && test "$answer" != "no" ; do
        read -p "$1 (enter yes to continue, no to stop the release process) " answer
    done
    test $answer = "yes" || complain "release process stopped at user request"
}

safety_checks()
{
    for cmd in git sed xsltproc sort tail curl stty base64 jq sha256sum gpg scp ; do
        which -s $cmd || complain "$cmd command not found"
    done
    git -C . rev-parse 2>/dev/null || complain "$(pwd) does not contain a git repository"
    test -f pom.xml                || complain "$(pwd)/pom.xml not found"
    test -d $1                     || complain "$(pwd)/$1 not found"
}

find_github_origin()
{
  git remote -v | sed -n "s,\([^ \t]*\)\t*.*@github.com:${github_owner}/${github_project}.git.*(push).*,\1,p"
}

find_current_branch()
{
    local start=$(git branch --show-current)
    test -z "$(git status --porcelain)" || complain "there are uncommitted changes in the branch"
    echo "current branch is $start" 1>&2
    echo $start
}

find_current_sha()
{
    git rev-parse --verify HEAD
}

find_pom_version()
{
    xsltproc scripts/get-pom-version.xsl $1
}

find_changes_version()
{
    xsltproc scripts/get-changes-version.xsl src/changes/changes.xml
}

extract_release_description()
{
    xsltproc scripts/changes2release.xsl src/changes/changes.xml
}

release_type()
{
  echo $1 | sed -e "s,^[0-9]*\.0$,major," -e "s,^[0-9]*\.[1-9][0-9]*$,minor," -e "s,^[0-9]*\.[0-9]*\.[0-9]*$,patch,"
}

release_branch()
{
  echo $1 | sed 's,\([0-9]\+.[0-9]\+\)[.0-9]*,release-\1,'
}

last_rc()
{
    local last="$(git tag -l "${1}-RC*" | sed 's,.*-RC,,' | sort -n | tail -1)"
    if test -z "$last" ; then
        echo 0
    else
        echo $last
    fi
}

rc_tag()
{
    echo "${1}-RC$2"
}

next_version()
{
    local current_major current_minor next_minor
    current_major=$(echo $1 | sed 's,\..*,,')
    current_minor=$(echo $1 | sed 's,[^.]*\.\([^.]\).*,\1,')
    next_minor=$(expr $current_minor + 1)
    echo ${current_major}.${next_minor}
}

modules()
{
  echo "hipparchus" \
       "hipparchus-clustering" "hipparchus-core"      "hipparchus-coverage" \
       "hipparchus-fft"        "hipparchus-filtering" "hipparchus-fitting"  \
       "hipparchus-geometry"   "hipparchus-ode"       "hipparchus-optim"    \
       "hipparchus-parent"     "hipparchus-samples"   "hipparchus-stat"
}

module_dir()
{
    if test "$1" = "hipparchus" ; then
        echo "."
    else
        echo "$1"
    fi
}

create_local_branch()
{
    test -z "$(git branch --list $1)" || complain "branch $1 already exists, stopping"
    git branch $1
    git checkout $1
    git merge --no-ff --no-commit $2
    if test ! -z "$(git status --porcelain)" ; then
        git status
        request_confirmation "commit merge from $2?"
        git commit -m "merging $2 to $1"
    fi
}

create_remote_branch()
{

    local created_branch="" timeout=0 current_date

    created_branch=$(github_get "$1" "git/refs/heads/$2" | jq --raw-output .object.sha)
    if test "$created_branch" = "null" ; then
        # release branch does not exist on origin yet, create it
        echo
        echo "creating remote branch $4/$2"
        github_post $1 "{ \"ref\": \"refs/heads/$2\", \"sha\": \"$3\" }" git/refs

        # waiting for remote branch to be available
        while test "$created_branch" = null ; do
          current_date=$(date +"%Y-%m-%dT%H:%M:%S")
          echo "${current_date} branch $2 not yet available in $4, waiting…"
          sleep 5
          timeout=$(expr $timeout + 5)
          test $timeout -lt 600 || complain "branch $2 not created in $4 after 10 minutes, exiting"
          created_branch=$(github_get $1 git/refs/heads/$2 | jq --raw-output .object.sha)
        done
        echo "branch $2 has been created"
        echo ""

    fi

}

delete_remote_branch()
{

     echo
     echo "deleting remote branch $2"
     github_delete $1 git/refs/heads/$2
     echo "branch $2 has been deleted"
     echo ""

}

merge_remote_branch()
{
    echo "merging $2 into $3"
    github_post $1 "{  \"head\": \"$2\", \"base\": \"$3\", \"commit_message\": \"$4\"}" merges
}

get_milestone_number()
{
    github_get $1 milestones | jq ".[] | select(.title==\"$2\") | .number"
}

delete_milestone()
{
    github_delete $1 milestones/$2
}

create_milestone()
{
    github_post $1 "{ \"title\": \"$2\" }" milestones
}

monitor_run()
{
    local run_id="" timeout=0 run_conclusion="null" current_date

    while test -z "$run_id" ; do
        current_date=$(date +"%Y-%m-%dT%H:%M:%S")
        echo "${current_date} waiting for pipeline to be triggered…"
        sleep 5
        timeout=$(expr $timeout + 5)
        test $timeout -lt 1800 || complain "pipeline not started after 30 minutes, exiting"
        run_id=$(github_get $1 actions/runs | jq ".workflow_runs[] | select(.head_sha==\"$2\" and .head_branch==\"$3\") | .id")
    done
    echo "run $run_id triggered"

    # monitor continuous integration action run (2 hours max)
    # conclusion is "null" while status is "in_progress" (i.e. not "completed")
    while test "$run_conclusion" != "success" -a "$run_conclusion" != "failure"  -a "$run_conclusion" != "cancelled" ; do
      current_date=$(date +"%Y-%m-%dT%H:%M:%S")
      echo "${current_date} pipeline ${run_id} status: ${run_conclusion}, waiting…"
      sleep 30
      timeout=$(expr $timeout + 30)
      test $timeout -lt 7200 || complain "pipeline not completed after 2 hours, exiting"
      run_conclusion=$(github_get $1 actions/runs/run_id | jq --raw-output ".conclusion")
    done

    if test "$run_conclusion" = "success" ; then
        return 0
    else
        return 1
    fi

}

enter_private_credentials()
{
    local credentials=""
    while test -z "$credentials" ; do
        echo "$1" 1>&2
        stty_orig=$(stty -g)
        stty -echo
        read credentials
        stty $stty_orig
    done
    echo "$credentials"
}

delete_local_and_remote_branch()
{
    if test ! -z "$2" ; then
        if test ! -z "$(git branch --list $2)" ; then
            git branch -D $2
            test "$(github_get "$1" "git/refs/heads/$2" | jq --raw-output .object.sha)" = "null" \
                || github_delete "$1" "git/refs/heads/$2"
        fi
    fi
}

delete_local_and_remote_tag()
{
    if test ! -z "$2" ; then
        if test ! -z "$(git tag --list $2)" ; then
            git tag -d $2
            test "$(github_get "$1" "git/refs/tags/$2" | jq --raw-output .object.sha)" = "null" \
                || github_delete "$1" "git/refs/tags/$2"
        fi
    fi
}

github_get()
{
  curl \
    --silent \
    --request GET \
    --header "Accept: application/vnd.github+json" \
    --header "Authorization: Bearer $1" \
    --header "X-GitHub-Api-Version: ${github_api_version}" \
  ${github_api}/$2
}

github_delete()
{
  curl \
    --silent \
    --output /dev/null \
    --request DELETE \
    --header "Accept: application/vnd.github+json" \
    --header "Authorization: Bearer $1" \
    --header "X-GitHub-Api-Version: ${github_api_version}" \
  ${github_api}/$2
}

github_post()
{
  curl \
    --silent \
    --output /dev/null \
    --request POST \
    --header "Accept: application/vnd.github+json" \
    --header "Authorization: Bearer $1" \
    --header "X-GitHub-Api-Version: ${github_api_version}" \
    --data   "$2" \
  ${github_api}/$3
}

find_forum_category()
{
    curl --silent \
         --request GET \
         --header "Content-Type: application/json" \
         --header "Api-Username: $1" \
         --header "Api-Key: $2" \
         --header "Accept: application/json" \
         https://${forum_fqdn}/categories \
    | jq ".category_list.categories[] | select(.name==\"$3\") | .id"
}

post_to_forum()
{
    curl --silent \
         --request POST \
         --header "Content-Type: application/json" \
         --header "Accept: application/json" \
         --header "Api-Username: $1"\
         --header "Api-Key: $2" \
         --data   "{ \"title\": \"$3\", \"raw\": \"$5\", \"category\": \"$3\" }" \
         https://${forum_fqdn}/posts.json \
    | jq --raw-output .post_url
}
