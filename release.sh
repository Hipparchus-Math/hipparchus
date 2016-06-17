# 
# Script to create a full release candidate (all modules).
#
# Run from top-level directory in full hipparchus checkout.
# Signed and hashed .bz2 files are created in /target.
#

# Config
keyID=0xCC6CB50E6A59DD12
version=1.0
rc=2
notes=false

# Temporarily drop -SNAPSHOT from version
find . -name pom.xml -exec sed -i '' 's/\-SNAPSHOT//g' {} \;

# Release notes - only top-level, refer to the online changelogs for modules
if $notes
	mkdir -p hipparchus-parent/src/changes
	cp src/changes/* hipparchus-parent/src/changes/
	mvn -Dversion=$version changes:announcement-generate -pl hipparchus-parent
	cp hipparchus-parent/target/announcement/release-notes.vm RELEASE-NOTES.txt
	rm -rf hipparchus-parent/src
    # Maybe hack the release notes a little ....
    read -p "Press [Enter] key start building tarballs"
fi

# Stage the release artifacts
mvn deploy -Dgpg.keyname=$keyID -Dwagon.skip=true -Ddescription=Hipparchus_${version}_RC${rc} -DskipStagingRepositoryClose=true -Prelease

# Cleanup target
# Add hashes and remove pom files 
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

# Undo pom changes
find . -name pom.xml -exec git checkout {} \;


