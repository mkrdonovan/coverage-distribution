##UPDATED:7/18/14
##Histogram is generated using pysam pileup, rather than genomeCoverageBed

import sys
import pysam
import argparse
import random
import time
import os
import shutil
import collections
import subprocess
import glob
import numpy as np
import math	
import logging
logging.basicConfig(format='[%(asctime)s][%(funcName)s][%(levelname)s] - %(message)s', level=logging.DEBUG)
logging.basicConfig(format='[%(asctime)s][%(funcName)s][%(levelname)s] - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


class downSampleBam(object):

	def __init__(self, BamFileList, tempDir, ignoreSmallCoverages):

		self.BamFileList = BamFileList
		self.tempDir = tempDir
		self.ignoreSmallCoverages = ignoreSmallCoverages
		
		self.referenceGenome = ''
		self.maxRefLength = 0
		self.coverages = []
		self.depthCounter = collections.defaultdict(int)

		
	def run(self, dsCoverage, picardPath, outputFolderName):
		
		self.coverageRun()
		dsCoverageToUse = self.checkDsCoverage(dsCoverage)
		self.downsample(dsCoverageToUse, picardPath)
		histFileList = self.generateHistogram(outputFolderName)
		return self.referenceGenome, histFileList


	def coverageRun(self):

		for bamfile in self.BamFileList:
			
			samfile = pysam.Samfile(bamfile, 'rb')
			self.findReferenceGenome(samfile)

			logger.info('Calculating coverage of: %s' %(os.path.basename(bamfile)))

			avgCoverage = self.averageCoverage(samfile)

			logger.info('Coverage: %s' %(avgCoverage))

			self.coverages.append(avgCoverage)


	def findReferenceGenome(self, samfile):
		
		maxRefLen = max([x['LN'] for x in samfile.header['SQ']])
		self.maxRefLength = maxRefLen * 1.

		for x in samfile.header['SQ']:
			if x['LN'] == maxRefLen:
				self.referenceGenome = x['SN']

		return self.referenceGenome


	def averageCoverage(self, samfile):

		totalCoverage = 0

		for base in samfile.pileup(self.referenceGenome):
			totalCoverage += base.n
			lastBasePos = base.pos

		return totalCoverage/float(lastBasePos+1)


	def checkDsCoverage(self, dsCoverage):
		
		if dsCoverage == 0:

			logger.info('Using minimum coverage: %s to downsample....' %(min(self.coverages)))
			return min(self.coverages)

		elif dsCoverage > min(self.coverages):
			
			if self.ignoreSmallCoverages == True:
				logger.info('Using given coverage:%sX to downsample... all samples with a coverage smaller than this will not be downsampled!' %(dsCoverage))
				return dsCoverage
			else:
				logger.info('Ignoring given coverage and using the minimum coverage: %s to downsample....' %(min(self.coverages)))
				return min(self.coverages)

		else:

			logger.info('Using given coverage:%sX to downsample....' %(dsCoverage))
			return dsCoverage
		

	def downsample(self, dsCoverageToUse, picardPath):
		
		x = 0 
		for bamfile in self.BamFileList:
			
			samfile = pysam.Samfile(bamfile, 'rb')
			probability = (dsCoverageToUse*1.0) / self.coverages[x]
			
			if probability == 1:
				logger.info('%s does not require downsampling... already at minimum coverage of %s' %(os.path.basename(bamfile), min(self.coverages)))
				shutil.copy(bamfile, os.path.join(self.tempDir, os.path.basename(os.path.splitext(bamfile)[0])+"_MinCoverage_"+str(int(round(dsCoverageToUse)))+"X"))
			
			elif probability > 1:
				logger.info('%s coverage is lower than the specified (%sX) coverage to downsample to... leaving at lower coverage (%sX)' %(os.path.basename(bamfile), dsCoverageToUse, self.coverages[x]))
				shutil.copy(bamfile, os.path.join(self.tempDir, os.path.basename(os.path.splitext(bamfile)[0])+"_underSpecCoverage_"+str(int(round(self.coverages[x])))+"X"))
			
			else:
				
			
				try:
					logger.info('Calling Picard for downsampling %s to %sX' %(os.path.basename(bamfile), dsCoverageToUse))
					subprocess.check_call(['java', '-Xmx1g', '-jar', os.path.join(picardPath, 'DownsampleSam.jar'), 'INPUT='+bamfile, 'OUTPUT='+os.path.join(self.tempDir, os.path.basename(os.path.splitext(bamfile)[0])+"_ds_"+str(int(round(dsCoverageToUse)))+"X"), 'PROBABILITY='+str(probability)], stdout=subprocess.PIPE)
					# logger.info('Picard called ... waiting for the downsampling to finish ...')
				
				except:

					print "Error reading BAM file:", sys.exc_info()[0]
					logger.info('%s will not be downsampled or plotted....' %(os.path.basename(bamfile)))
					os.remove(os.path.join(self.tempDir, os.path.basename(os.path.splitext(bamfile)[0])+"_ds_"+str(int(round(dsCoverageToUse)))+"X"))


				# logger.info('Calling Picard for downsampling %s to %sX' %(os.path.basename(bamfile), dsCoverageToUse))
				# dwnsmpl = subprocess.Popen(['java', '-Xmx1g', '-jar', os.path.join(picardPath, 'DownsampleSam.jar'), 'INPUT='+bamfile, 'OUTPUT='+os.path.join(self.tempDir, os.path.basename(os.path.splitext(bamfile)[0])+"_ds_"+str(int(round(dsCoverageToUse)))+"X"), 'PROBABILITY='+str(probability)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				# logger.info('Picard called ... waiting for the downsampling to finish ...')
				# output, error = dwnsmpl.communicate()
				
			x+=1

	def generateHistogram(self, outputFolderName):

		histogramFileList = []
		
		for bamfDS in glob.glob(os.path.join(str(self.tempDir), '*X')):

			pysam.index(bamfDS)
			base = os.path.basename(bamfDS)
			
			histogramFile = os.path.join(outputFolderName, (base)+"_hist.txt")
			
			with open(histogramFile, 'w') as fout:
			
				histogramFileList.append(histogramFile)

				samfile = pysam.Samfile(bamfDS, 'rb')
				logger.info('Generating histogram for %s' %(base))
				
				a = []
				pc = samfile.pileup(self.referenceGenome)
				for p in pc:
					a.append((p.pos, (p.n)*1.0))

				ad = [x[1] for x in a]
				
				for depth in ad:
					self.depthCounter[depth]+=1

				
				for depth, numBases in self.depthCounter.items():

					fout.write('%i\t%f\n' %(depth, numBases/self.maxRefLength))
				
				self.depthCounter = collections.defaultdict(int)

		return histogramFileList

class makeCoveragePlot(object):

	def __init__(self, tempDir, outputFolder, histogramFileList):
		
		self.tempDir = str(tempDir)
		self.outputFolder = str(outputFolder)
		self.histogramFileList = histogramFileList
		# self.BamsToPlot = []

		self.linePatterns = ['#ed6161','#76baf5', '#ebb970', '#74b993', '#c19ccd', '#9d0000', '#003d71', '#cf7b00', '#00803a', '#75009b', '#f13232', '#42a2f5', '#f5a127', '#20cd6e', '#b96ad3']
		self.patternCounter = 0
		self.IQRlist = []

	def run(self, plotOut, referencegenome):
		
		if self.histogramFileList:
			logger.info('Beginning to generate histograms....')
			self.plot(plotOut, referencegenome)
			self.outputIQRfile(self.IQRlist, 'IQR_out.txt')

	def plot(self, plotOut, referencegenome):

		fig = plt.figure('Project_Coverage_Plot')
		ax = fig.add_subplot(1,1,1)
		
		depthLimits=[]
		ylimit=[]

		for histFile in self.histogramFileList:

			IQRlisttemp = []
			
			coverageCoordinates = self.getCoverageCoordinates(histFile, referencegenome)
			depth = coverageCoordinates[0]
			percent = (coverageCoordinates[1])

			basehist = os.path.basename(histFile)
			
			ax.plot(depth, percent, self.linePatterns[self.patternCounter], linewidth = 2.5, label = os.path.basename(os.path.splitext(basehist)[0]))
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
		fig.savefig(os.path.join(self.outputFolder,plotOut), format = 'PDF')

	def getCoverageCoordinates(self, histFILE, referencegenome):
		
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


		fn = open(os.path.join(self.outputFolder,"IQR.txt"), 'w')

		fn.write("#IQR Summary\tIQR\n")

		for iqrItem in IQRlist:
			fn.write("\t".join(iqrItem)+"\n")

		fn.close()


def checkPaths(BamFileList, picardPath):

	if not os.path.exists(picardPath):
		logger.info('Picard path specified does not exist. Exiting program!')
		sys.exit()

	for bamFile in BamFileList:
		if not os.path.exists(bamFile):
			logger.info('%s does not exist. Removing this file from the list of bam files to be analyzed....' %(bamFile))
			BamFileList.remove(bamFile)

	for bamFile in BamFileList:
		if not (glob.glob(bamFile+'*bai'))[0]:
			logger.info('%s is not indexed. Indexing now....' %(bamFile))
			pysam.index(bamFile)

	return BamFileList

def main(BamFileList, picardPath, dsCoverage, plotOut, ignoreSmallCoverages, outputFolderName):
	
	BamFileList = checkPaths(BamFileList, picardPath)

	tempDir = ('tmp_' + str(random.randint(0,int(time.time()))))
	os.makedirs(tempDir)

	
	if not os.path.exists(outputFolderName):
		os.makedirs(outputFolderName)

	try:
		ds = downSampleBam(BamFileList, tempDir, ignoreSmallCoverages)
		out = ds.run(dsCoverage, picardPath, outputFolderName)

		referencegenome = out[0]
		histogramFileList = out[1]

		plot = makeCoveragePlot(tempDir, outputFolderName, histogramFileList)
		plot.run(plotOut, referencegenome)

	except KeyboardInterrupt:

		logger.info('KeyboardInterrupt, removing temp file....')
		rmDIR = subprocess.Popen(['rm', '-r', tempDir], stderr=subprocess.PIPE)
		output, error = rmDIR.communicate()


	if os.path.exists(tempDir):
		logger.info('Removing temp file....')
		rmDIR = subprocess.Popen(['rm', '-r', tempDir], stderr=subprocess.PIPE)
		output, error = rmDIR.communicate()
		if error:
			logger.error('Not removing temp file...')


def plotHistsOnly(bamFileList, outputFolder, plotFile, coverage):

	histogramFileList = []

	for bamf in bamFileList:
		base = os.path.basename(os.path.splitext(bamf)[0])
		histogramFile = glob.glob(os.path.join(outputFolder, (base)+'*'+str(coverage)+'*_hist.txt'))
		histogramFileList.append(histogramFile[0])

	samfile = pysam.Samfile(bamFileList[0], 'rb')
	getRefGenome = downSampleBam(bamFileList, '', False)
	referenceGenome = getRefGenome.findReferenceGenome(samfile)

	plotAgain = makeCoveragePlot(plotFile, outputFolder, histogramFileList)
	plotAgain.plot(plotFile, referenceGenome)




if __name__ == "__main__":

	dsFactor = 0

	parser = argparse.ArgumentParser()
	
	parser.add_argument('--only_plot', action='store_true', help='if you only want to plot histrogram files you have already generated, input the output_folder where the histogram files are and indicate which bamfiles you want to plot')
	parser.add_argument('--output_folder', dest='outputFolderName', type=str, required=True, help='folder where your plot, IQR file, and histrogram files will be stored.')
	parser.add_argument('--picard', dest='picardPath', type=str, required=False, default='/illumina/thirdparty/picard-tools/picard-tools-1.85/', help='Path to picard tools software')
	parser.add_argument('--ignore_small_coverages', action='store_true', help = 'Ignore samples (do not downsample, but still plot) that have coverages that are smaller than your choosen ds coverage')
	parser.add_argument('--coverage', dest='dsCoverage', type=int, required=False, default=0, help='Coverage after downsampling')
	parser.add_argument('--plot_file', dest='plotOut', type=str, required=False, default='Coverage_Plot.pdf', help='file you want to output you plot to')
	parser.add_argument('--bamfiles', dest='BamFileList', nargs='+', type=str, required=True, help='Bam file(s) to analyze/plot (if ONLY plotting HIST files, indicicate here which ones you want to plot.')

	args = parser.parse_args()


	if args.only_plot == True:
		plotHistsOnly(args.BamFileList, args.outputFolderName, args.plotOut, args.dsCoverage)
	else:
		main(args.BamFileList, args.picardPath, args.dsCoverage, args.plotOut, args.ignore_small_coverages, args.outputFolderName)