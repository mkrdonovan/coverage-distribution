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

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


class downSampleBam(object):

	def __init__(self, BamFileList, ignoreSmallCoverages, chrom):

		self.BamFileList = BamFileList
		self.FinalBamList = []
		self.ignoreSmallCoverages = ignoreSmallCoverages
		
		self.referenceGenome = chrom
		self.maxRefLength = 0
		self.coverages = []
		self.depthCounter = defaultdict(int)

		
	def run(self, dsCoverage, picardPath, histogramFolder):
		
		histFileList = self.generateHistogram(histogramFolder)
		return histFileList


	def generateHistogram(self, outputFolderName):

		histogramFileList = []
		
		for bam in self.BamFileList:

			pysam.index(bam)
			base = os.path.basename(bam)
			histogramFile = os.path.join(outputFolderName, self.referenceGenome + "_" + base +"_hist.txt")
			
			with open(histogramFile, 'w') as fout:
			
				histogramFileList.append(histogramFile)

				samfile = pysam.Samfile(bam, 'rb')
				logger.info('Generating histogram for %s' %(base))
				
				for x in samfile.header['SQ']:
					if x['SN'] == self.referenceGenome:
						self.maxRefLength = (x['LN'])*1.0

				a = []
				pc = samfile.pileup(self.referenceGenome)
				for p in pc:
					a.append((p.pos, (p.n)*1.0))

				ad = [x[1] for x in a]
				
				for depth in ad:
					self.depthCounter[depth]+=1

				for depth, numBases in self.depthCounter.items():

					fout.write('%i\t%f\n' %(depth, numBases/self.maxRefLength))
				
				self.depthCounter = defaultdict(int)

		return histogramFileList


class makeCoveragePlot(object):

	def __init__(self, outputFolder, histogramFileList, refGen):
		
		self.outputFolder = str(outputFolder)
		self.histogramFileList = histogramFileList
		self.refgen = refGen
		# self.BamsToPlot = []

		self.linePatterns = ['#ed6161','#76baf5', '#ebb970', '#74b993', '#c19ccd', '#9d0000', '#003d71', '#cf7b00', '#00803a', '#75009b', '#f13232', '#42a2f5', '#f5a127', '#20cd6e', '#b96ad3']
		self.patternCounter = 0
		self.IQRlist = []

	def run(self, plotOut, referencegenome):
		
		if self.histogramFileList:
			logger.info('Beginning to generate plots....')
			self.plot(plotOut, referencegenome)
			self.outputIQRfile(self.IQRlist, referencegenome + '_IQR_out.txt')

	def plot(self, plotOut, referencegenome):

		fig = plt.figure('Project_Coverage_Plot_' + referencegenome)
		# plt.rcParams['font.family']='Sawasdee'
		ax = fig.add_subplot(1,1,1)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()
		plt.title('Chromosome: %s' %self.refgen)
		ax.spines["left"].axis.axes.tick_params(direction="outward")
		ax.spines["bottom"].axis.axes.tick_params(direction="outward")
		depthLimits=[]
		ylimit=[]

		for histFile in self.histogramFileList:

			IQRlisttemp = []
			
			coverageCoordinates = self.getCoverageCoordinates(histFile)
			depth = coverageCoordinates[0]
			percent = (coverageCoordinates[1])
			depth.insert(0,0)
			percent.insert(0,0)

			basehist = os.path.basename(histFile).split("_")
			ax.plot(depth, percent, self.linePatterns[self.patternCounter], linewidth = 2.5, label = "_".join(basehist[:-1]))
			leg = ax.legend(loc = 'upper right', bbox_to_anchor=(0,0,1,1), prop ={'size': 6}, frameon=False, shadow=False)
			
			self.patternCounter += 1

			if self.patternCounter == len(self.linePatterns):
				self.patternCounter = 0

			
			IQR = str(self.getCDF(percent))
			IQRlisttemp.append(os.path.basename(histFile))
			IQRlisttemp.append(IQR)
			self.IQRlist.append(IQRlisttemp)

			perDepthList= zip(depth, percent)
			for coordinate in reversed(perDepthList):
				if coordinate[1] >= 0.01:
					depthLimits.append(int(coordinate[0]))
				if coordinate[1] > 0:
					ylimit.append(float(coordinate[1]))
		
		plt.xlim([0, max(depthLimits) + (max(depthLimits)*.20)])
		plt.ylim([0, max(ylimit) + (max(ylimit)*.15)])
		plt.xlabel('Depth of Coverage')
		plt.ylabel('% Bases')	
		fig.savefig(os.path.join(self.outputFolder,referencegenome + "_" + plotOut), format = 'PDF')

	def getCoverageCoordinates(self, histFILE):
		
		hisList=[]
		logger.info('Plotting histogram file')
		
		
		with open(histFILE, 'rb') as fobj:
			for line in fobj:

				hisList.append(((line.strip()).split("\t")))

		depth=[]
		percent=[]
		for x in hisList:
			depth.append(x[0])
			percent.append((float(x[1]))*100.)
		
		return depth, percent

	def getCDF(self, percent):
	
		percent = np.array(percent, dtype=float)
		cdf = np.cumsum(percent)
		x = np.interp([0.25, 0.75], cdf, range(len(cdf)))
		IQR = x[1] - x[0]
		
		return IQR

	def outputIQRfile(self, IQRlist, outputName):


		fn = open(os.path.join(self.outputFolder,outputName), 'w')

		fn.write("#IQR Summary\tIQR\n")

		for iqrItem in IQRlist:
			fn.write("\t".join(iqrItem)+"\n")

		fn.close()

def checkPaths(BamFileList, picardPath):

	if not os.path.exists(picardPath):
		logger.info('Picard path specified does not exist. Exiting program!')
		exit()

	for bamFile in BamFileList:
		if not os.path.exists(bamFile):
			logger.info('%s does not exist. Removing this file from the list of bam files to be analyzed....' %(bamFile))
			BamFileList.remove(bamFile)

	for bamFile in BamFileList:
		if not (glob(bamFile+'*bai'))[0]:
			logger.info('%s is not indexed. Indexing now....' %(bamFile))
			pysam.index(bamFile)

	return BamFileList

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


# def removeDuplicates(BamFileList, tempDir):
	
# 	logger.info('Removing optical and PCR duplicates....')
# 	noDupBamFileList = []
# 	for bamfile in BamFileList:

# 		fName = os.path.join(tempDir, os.path.basename(os.path.splitext(bamfile)[0])+"_noDup")
# 		noDupBamFileList.append(fName)
# 		pysam.rmdup(bamfile, fName)
# 		pysam.index(fName)
	
# 	return noDupBamFileList


def downsampleAllBams(BamFileList, picardPath, dsCoverage, ignoreSmallCoverages, outputFolderName, chrToAnalyze, histogramFolder, plotFile):

	BamFileList = checkPaths(BamFileList, picardPath)

	if not os.path.exists(histogramFolder):
		os.makedirs(histogramFolder)

	

	chromosomesToAnalyze = chooseChrs(chrToAnalyze, BamFileList)
	# noDupBamFileList = removeDuplicates(BamFileList, tempDir)
	for chrom in chromosomesToAnalyze:
		ds = downSampleBam(BamFileList, ignoreSmallCoverages, chrom)
		histogramFileList = ds.run(dsCoverage, picardPath, histogramFolder)
		
		plot = makeCoveragePlot(outputFolderName, histogramFileList, chrom)
		plot.run(plotFile, chrom)

	# except KeyboardInterrupt:

	# 	logger.info('KeyboardInterrupt, removing temp file....')
	# 	rmDIR = subprocess.Popen(['rm', '-r', tempDir], stderr=subprocess.PIPE)
	# 	output, error = rmDIR.communicate()

	# if os.path.exists(tempDir):
	# 	logger.info('Removing temp file....')
	# 	rmDIR = subprocess.Popen(['rm', '-r', tempDir], stderr=subprocess.PIPE)
	# 	output, error = rmDIR.communicate()
	# 	if error:
	# 		logger.error('Not removing temp file...')


def plotHistsOnly(bamFileList, outputFolder, plotFile, coverage, chrToAnalyze, histogramFolder):

	# histogramFileList = []
	# chromosomesToAnalyze = chooseChrs(chrToAnalyze, bamFileList)

	# for bamf in bamFileList:
	# 	base = os.path.basename(os.path.splitext(bamf)[0])

	# 	for chrm in chromosomesToAnalyze:

	# 		histogramFile = glob(os.path.join(histogramFolder, chrm + "*" + base + "_hist.txt"))
	# 		histogramFileList.append(histogramFile[0])

	# samfile = pysam.Samfile(bamFileList[0], 'rb')	
	plotAgain = makeCoveragePlot(outputFolder, bamFileList, ['chrY'])
	plotAgain.plot(plotFile, chrToAnalyze)



if __name__ == "__main__":

	dsFactor = 0

	parser = ArgumentParser()
	
	parser.add_argument('--generate-histograms', dest='downsample', action='store_true')
	parser.add_argument('--replot-histograms', dest='plotHist', action='store_true', help='if you only want to plot histrogram files you have already generated, input the output_folder where the histogram files are and indicate which bamfiles you want to plot')

	parser.add_argument('--ignore-small-coverages', dest='ignoreSmallCoverages', action='store_true', help = 'Ignore samples (do not downsample, but still plot) that have coverages that are smaller than your choosen ds coverage')
	parser.add_argument('--coverage', dest='dsCoverage', type=int, required=False, default=0, help='Coverage after downsampling')

	parser.add_argument('--out_folder', dest='outputFolderName', type=str, required=True, help='folder where your plot, IQR file, and histrogram files will be stored.')
	parser.add_argument('--chr', dest='chrToAnalyze', type=str, default='LONG', help='Choose to analyze largest chromosome ("LONG"), choose to analyze all chromosomes ("ALL"), or specify the chromosome name to analyze')
	parser.add_argument('--picard', dest='picardPath', type=str, required=False, default='/illumina/thirdparty/picard-tools/picard-tools-1.85/', help='Path to picard tools software')
	
	parser.add_argument('--plot-file', dest='plotOut', type=str, required=False, default='Coverage_Plot.pdf', help='file you want to output you plot to')
	parser.add_argument('--bamfiles', dest='BamFileList', nargs='+', type=str, required=True, help='Bam file(s) to analyze/plot (if ONLY plotting HIST files, indicicate here which ones you want to plot.')

	args = parser.parse_args()

	histogramFolder = os.path.join(args.outputFolderName, 'histograms')

	if args.downsample == False and args.plotHist == False:
		logger.info('Did not choose an action. Please select --generate-histograms or --plot-histograms.')

	if args.downsample == True:
		downsampleAllBams(args.BamFileList, args.picardPath, args.dsCoverage, args.ignoreSmallCoverages, args.outputFolderName, args.chrToAnalyze, histogramFolder, args.plotOut)
	
	if args.plotHist == True:
		if not os.path.exists(histogramFolder):
			logger.info('Histogram files not yet generated. First select --downsample to generate histograms for a desired coverage.')
		else:
			plotHistsOnly(args.BamFileList, args.outputFolderName, args.plotOut, args.dsCoverage, args.chrToAnalyze, histogramFolder)

