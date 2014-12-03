from sys import exit
import pysam
from argparse import ArgumentParser
from random import randint
from time import time
import os
from shutil import copy
from collections import defaultdict
import subprocess
from glob import glob
import numpy as np	
import logging
logging.basicConfig(format='[%(asctime)s][%(funcName)s][%(levelname)s] - %(message)s', level=logging.DEBUG)
logging.basicConfig(format='[%(asctime)s][%(funcName)s][%(levelname)s] - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def averageCoverage(BamFileList, chrom):

	for bamfile in BamFileList:
		samfile = pysam.Samfile(bamfile, 'rb')
		flag = checkChr(samfile, chrom)
		if flag == True:
			totalCoverage = 0

			for base in samfile.pileup(chrom):
				totalCoverage += base.n
				lastBasePos = base.pos

			writeSummary(totalCoverage/float(lastBasePos+1), bamfile, chrom)

def writeSummary(coverage, bamPath, chromosome):

	with open('coverage_summary'+str(chromosome)+'.csv', 'a') as fout:
		fout.write("%s,%s\,%s\n" % (bamPath, chromosome, coverage))

def checkChr(samfile, chrom):
	
	flag = False
	for x in samfile.header['SQ']:
		if x['SN'] == chrom:
			flag = True
	return flag

def chooseChrs(chrToAnalyze, BamFileList):

	chromosomesToAnalyze = []
	samfile = pysam.Samfile(BamFileList[0], 'rb')

	if chrToAnalyze.upper() == 'ALL':
		logger.info('Analyzing all chromosomes:')
		for x in samfile.header['SQ']:
			chromosomesToAnalyze.append(x['SN'])

	elif chrToAnalyze.upper() == 'LONG':
		logger.info('Analyzing the largest chromosome:')
		maxRefLen = (max([x['LN'] for x in samfile.header['SQ']])) * 1.0
		for x in samfile.header['SQ']:
			if x['LN'] == maxRefLen:
				chromosomesToAnalyze.append(x['SN'])

	else:
		flag = False
		for x in samfile.header['SQ']:
			if chrToAnalyze == x['SN']:
				flag = True
		if flag == True:
			logger.info('Analyzing chromosome:')
			chromosomesToAnalyze.append(chrToAnalyze)
		else:
			logger.info('Chosen chromosome is not an option in the given bam files, please select another.')
			exit()

	print ', '.join(chromosomesToAnalyze)
	return chromosomesToAnalyze


def main(chrToAnalyze, BamFileList):
	
	chromosomesToAnalyze = chooseChrs(chrToAnalyze, BamFileList)
		
	for chrom in chromosomesToAnalyze:
		averageCoverage(BamFileList, chrom)



if __name__ == "__main__":

	dsFactor = 0

	parser = ArgumentParser()
	
	parser.add_argument('--chr', dest='chrToAnalyze', type=str, default='LONG', help='Choose to analyze largest chromosome ("LONG"), choose to analyze all chromosomes ("ALL"), or specify the chromosome name to analyze')
	parser.add_argument('--bamfiles', dest='BamFileList', nargs='+', type=str, required=True, help='Bam file(s) to analyze/plot (if ONLY plotting HIST files, indicicate here which ones you want to plot.')

	args = parser.parse_args()

	main(args.chrToAnalyze, args.BamFileList)