#!/usr/bin/env python2.7
from setuptools import setup, find_packages

def readme():
	with open('README.rst') as f:
		return f.read()

setup(	name = 'GRAViTy',
	version = '1.1.0',
	description = 'GRAViTy: Genome Relationships Applied to Virus Taxonomy',
	long_description = readme(),
	classifiers = [
		#'Development Status :: 3 - Alpha',
		#'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2.7',
		#'Topic :: Text Processing :: Linguistic',
	],
	keywords = [
		'Viruses',
		'Virus classification', 'Virus classifier',
		'ICTV', 'International Committee on Taxonomy of Viruses',
		'GRAViTy',
		'Metagenomics',
	],
	#url = 'https://github.com/PAiewsakun/GRAViTy',
	author = 'Pakorn Aiewsakun, Peter Simmonds',
	author_email = 'pakorn.aiewsakun@gmail.com, peter.simmonds@ndm.ox.ac.uk',
	#license = 'MIT',
	packages = ['GRAViTy','GRAViTy.Utilities'],
	#packages = find_packages(),
	scripts = [
		'bin/booster',
		'bin/ReadGenomeDescTable',
		'bin/PPHMMDBConstruction',
		'bin/RefVirusAnnotator',
		'bin/GRAViTyDendrogramAndHeatmapConstruction',
		'bin/MutualInformationCalculator',
		'bin/UcfVirusAnnotator',
		'bin/VirusClassificationAndEvaluation',
		'bin/GRAViTy_Pipeline_I',
		'bin/GRAViTy_Pipeline_II',
	],
	#dependency_links = ['http://github.com/user/repo/tarball/master#egg = package-1.0'],
	install_requires = [
		'Biopython<=1.76',
		'numpy',
		'ete3',
		'matplotlib',
		'scikit-learn',
		#'dcor',
		'scipy',
		'dendropy'
	],
	include_package_data = True,
	zip_safe = False)

