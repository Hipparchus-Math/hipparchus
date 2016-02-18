#!/usr/bin/env python

import argparse
import os
import tempfile
import filecmp
import shutil
import re
import atexit

def is_valid_directory(parser, arg):
    if not os.path.isdir(arg):
        parser.error('directory {} does not exist!'.format(arg))
    else:
        # File exists so return the directory
        return arg

def fnmatch(fname, extensions):
	for ext in extensions:
		if fname.endswith(ext):
			return True
	return False

def processDir(rootDir, args):
	for dirName, subdirList, fileList in os.walk(rootDir):
		for fname in fileList:
			if fnmatch(fname, args.ext):
				fullName = os.path.join(dirName, fname)
				processFile(fullName, args)
		
		# do not traverse into subdirs that match an ignore
		for ignore in args.ignore:
			index = subdirList.index(ignore)
			if index >= 0:
				del subdirList[index]


def processFile(file, args):
	g = os.path.basename(file)
	tempfile = os.path.join(args.tempdir, g)
	with open(file, "r") as sources:
		lines = sources.readlines()
	with open(tempfile, "w") as sources:
		for line in lines:
			sources.write(packageReplace(line, args))

	if filecmp.cmp(file, tempfile) == False:
		if args.dry_run == True:
			print file
			os.remove(tempfile)
		else:
			if args.save == True:
				print "saving file" + file
				shutil.move(file, file + ".orig")

			shutil.move(tempfile, file)

def packageReplace(line, args):
	for pattern in args.packages_subst_pattern:
		p = pattern.split(':')
		line = re.sub(p[1], p[2], line)
	return line

def parsePackageSubstitutionRules(args):
	with open(args.packages_subst, "r") as sources:
		lines = sources.readlines()

	for line in lines:
		pattern = re.sub(r'^\s*"([^"]*)"\s*"([^"]*)"\s*', r's:\1:\2', line)
		if pattern.startswith('s:'):
			args.packages_subst_pattern.append(pattern)


# start main
# TODO: handle class substitutions

parser = argparse.ArgumentParser(prog='migrate')
parser.add_argument('--dir', action='append', required=True, type=lambda x: is_valid_directory(parser, x))
parser.add_argument('--ignore', action='append', default=[])
parser.add_argument('--ext', action='append', required=True)
parser.add_argument('--nosave', action='store_false', dest='save', default=True)
parser.add_argument('--packages-subst', default=None)
parser.add_argument('--classes-subst', default=None, required=True)
parser.add_argument('--dry-run', action='store_true', default=False)

args = parser.parse_args()

tempdir = tempfile.mkdtemp()
atexit.register(lambda dir=tempdir: shutil.rmtree(dir))
args.tempdir = tempdir

args.packages_subst_pattern = []
if args.packages_subst == None:
	print 'WARNING: package substitutions disabled at user request'
else:
	parsePackageSubstitutionRules(args)

#print args

# process the files
for dir in args.dir:
	print 'Processing files in directory ' + dir
	processDir(dir, args)
