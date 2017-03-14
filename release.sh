#!/bin/sh
#
# Licensed to the Hipparchus under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# Hipparchus licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------
#
# Script to create a full release candidate (all modules).
#
# Run from top-level directory in full hipparchus checkout.
# Signed and hashed .bz2 files are created in /target.
#

# the keyID is not configured in the script anymore; it should be in
# each release manager own settings.xml file as a property as follows:
# <profiles>
#   <profile>
#      <id>hipparchus-gpg-settings</id>
#      <activation>
#        <file>
#          <exists>${basedir}/hipparchus-core</exists>
#        </file>
#      </activation>
#      <properties>
#        <gpg.keyname>the 16 hexadecimal digits representing the manager keyID</gpg.keyname>
#      </properties>
#    </profile>
#  </profiles>


# Config
# Set to true to generate top-leve release notes draft
notes=false

# Set to true to drop -SNAPSHOT from versions in poms
kill_snap=false

# Drop -SNAPSHOT from version. Set kill_snap=true to execute this
# on a newly created release branch.
if $kill_snap
then
    find . -name pom.xml -exec sed -i '' 's/\-SNAPSHOT//g' {} \;
    # Pause to commit the change
    read -p "Commit pom changes to release branch, then press [Enter] to continue."
fi

# Top-level release notes. Set notes=true to execute this, typically just for
# the first RC.
if $notes
then
	mkdir -p hipparchus-parent/src/changes
	cp src/changes/* hipparchus-parent/src/changes/
	mvn -Dversion=$version changes:announcement-generate -pl hipparchus-parent
	cp hipparchus-parent/target/announcement/release-notes.vm RELEASE-NOTES.txt
	rm -rf hipparchus-parent/src
    # Maybe hack the release notes a little ....
    read -p "Edit, add module-specific changes and check RELEASE-NOTES.txt in. Then Press [Enter] key start building tarballs."
fi

# Stage the release artifacts
mvn clean
mvn deploy -DskipStagingRepositoryClose=true -Prelease

# Add hashes and remove pom files in target
cd target
for f in *.zip
do
	md5 $f > $f.md5
done
for f in *.bz2
do
	md5 $f > $f.md5
done
rm -rf archive-tmp*
rm *.pom
rm *.pom.asc
cd ..
