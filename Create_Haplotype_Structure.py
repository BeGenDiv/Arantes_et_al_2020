import argparse
import gzip
import os


parser = argparse.ArgumentParser(description="")
parser.add_argument("popMap", help="population map used for stacks",type=str)
parser.add_argument("mappedCounts", help="tab separated file containing the readcounts per locus(for each sample e.g. sampleID<tab>locID<tab>readCounts)",type=str)
parser.add_argument("popOutDir", help="output directory of the stacks population module",type=str)
parser.add_argument("outDir", help="output directory",type=str)
parser.add_argument("--readCOV", help="minimum read coverage needed for a locus [default=6]", default=6, type=int)
parser.add_argument("--popCOV", help="total number of populations that needs to cover a locus to include it in the structure output [default=4]", default=4, type=int)
parser.add_argument("--intraCOV", help="percentage of each population that needs to cover a locus to include it in the structure output [default=0.7]", default=0.7, type=float)
parser.add_argument("--intraCOVtotal", help="enable this option to have a strict int value as input for --intraCOV --> e.g. all need to by covered by at least 2 samples [default=disabled]", action="store_true")
args = parser.parse_args()


############################
####### Parse popmap #######
############################
sampleDICT={}
popDICT_popSize={}
popDICT_ID={}
popDICT={}
totalSamples=0

print("1) Parsing popmap")
myReader=open(args.popMap)
# get sample popIds and pop sizes
for line in myReader:
	sID, sPop = line.strip().split("\t")

	if popDICT_popSize.get(sPop, "nope") == "nope":
		popDICT_popSize[sPop] = 1
	else:
		popDICT_popSize[sPop] += 1

	if sampleDICT.get(sID, "nope")=="nope":
		sampleDICT[sID]=sPop
	#else:
	#	print("Weird sample %s is already in dict!!!"%(sID))
	if popDICT.get(sPop, "nope") == "nope":
		popDICT[sPop]=[sID]
	else:
		popDICT[sPop].append(sID)

#print(popDICT_popSize)
popCounts=1
for pop in popDICT_popSize:
	if popDICT_ID.get(pop, "nope")== "nope":
		popDICT_ID[pop] = popCounts
		popCounts+=1
	else:
		print("Weird pop %s is already in popDICT_ID!!!"%(pop))

myReader.close()


print("\t#Populations:%i"%(len(popDICT_popSize)))
for pop in popDICT_popSize:
	print("\t\t%s: %i"%(pop, popDICT_popSize[pop]))
	totalSamples+=popDICT_popSize[pop]


###################################################################################
####### PARSE: populations.sumstats.tsv just to get the contig to locus IDs #######
###################################################################################
print("2) Parsing populations.sumstats.tsv")
# to store original locusID to stacksLocusID
contigID_locusID_DICT={}

for file in os.listdir(args.popOutDir):

	if file.endswith("populations.sumstats.tsv"):
		
		myReader=open(args.popOutDir+"/"+file)
		for line in myReader:
			if not line.startswith("#"):
				splitted=line.strip().split()
				curLoc = int(splitted[0])
				contigID=splitted[1]

				#save contig and stacks locus ID
				if contigID_locusID_DICT.get(curLoc, "nope") == "nope":
					contigID_locusID_DICT[curLoc]=contigID


		myReader.close()


###################################################
###### GET MAPPED COUNTS PER REF PER SAMPLE #######
###################################################
print("3) Parsing readCounts")
mappedReader=open(args.mappedCounts)

readCountDICT={}
#1. key locusID --> 2. key popID --> contains list with all mapped counts for each sample of a pop per locus --> for easy checking of popCOVintra
refPop_countDICT={}

for line in mappedReader:
	sID, contigID, rCounts = line.strip().split("\t")

	
	if readCountDICT.get(sID, "nope") == "nope":
		readCountDICT[sID]={}
		readCountDICT[sID][contigID]=int(rCounts)/2
	else:
		if readCountDICT[sID].get(contigID, "nope") == "nope":
			#divided by 2 because of paired end data
			readCountDICT[sID][contigID]=int(rCounts)/2
		else:
			print("Weird already in count dict...")




	#get population for a species
	#if just to check if a pop is in sampleDICT --> its not for the 3 bad samples so we dont have to create a new readcounts file
	if sampleDICT.get(sID, "nope") != "nope":
		curPop = sampleDICT[sID]
		if refPop_countDICT.get(contigID, "nope")=="nope":
			refPop_countDICT[contigID]={}
			#divided by 2 because paired end reads R1 and R2 only overlap with other R1 or R2
			refPop_countDICT[contigID][curPop]=[int(rCounts)/2]
		else:
			if refPop_countDICT[contigID].get(curPop, "nope")=="nope":
				refPop_countDICT[contigID][curPop]=[int(rCounts)/2]
			else:
				refPop_countDICT[contigID][curPop].append(int(rCounts)/2)



mappedReader.close()

#print(readCountDICT)


##################################################
####### PARSE:  populations.haplotypes.tsv #######
##################################################
print("4) Parsing populations.haplotypes.tsv")
# parse haplotypes for all samples across all loci
for file in os.listdir(args.popOutDir):

	if file.endswith("populations.haplotypes.tsv"):
		
		myReader=open(args.popOutDir+"/"+file)
		myWriter=open(args.outDir+"/populations.cov_filtered.haplotypes.tsv", "w")
		myWriter2=open(args.outDir+"/populations.cov_filtered_noNs.haplotypes.tsv", "w")
		myRemovedW1=open(args.outDir+"/populations.cov_filtered.haplotypess_removed_COV.tsv", "w")
		myRemovedW2=open(args.outDir+"/populations.cov_filtered.haplotypess_removed_Ns.tsv", "w")
		myHeader=[]
		locusDICT={}
		header=""


		myRemovedW1.write("#sample\tlocID\tstacksLocID\thaplo\treadCOV\n")
		myRemovedW2.write("#sample\tlocID\tstacksLocID\thaplo\treadCOV\n")

		for line in myReader:
			if line.startswith("#"):
				myWriter.write(line)
				myWriter2.write(line)
				header=line.strip().split("\t")
			else:
				splitted=line.strip().split("\t")
				curLoc=int(splitted[0])

				
				if contigID_locusID_DICT.get(curLoc, "nope") != "nope":

					curCont=contigID_locusID_DICT[curLoc]
					myWriter.write("%s\t%s"%(splitted[0], splitted[1]))
					myWriter2.write("%s\t%s"%(splitted[0], splitted[1]))
					for i in range(2, len(splitted)):
						curHaplo=splitted[i]
						curSample=header[i]

						if readCountDICT[curSample].get(curCont, "nope") == "nope":
							myWriter.write("\t-")
							myWriter2.write("\t-")

						else:
							#print(curSample, curLoc, curCont, readCountDICT[curSample][curCont], args.readCOV)
							if readCountDICT[curSample][curCont] >= args.readCOV:
								myWriter.write("\t%s"%(splitted[i]))
								if not "N" in splitted[i].upper():
									myWriter2.write("\t%s"%(splitted[i]))
								else:
									myWriter2.write("\t-")
									#write removed haplotype by N
									myRemovedW2.write("%s\t%s\t%s\t%s\t%i\n"%(curSample, curCont, curLoc, splitted[i], readCountDICT[curSample][curCont]))
							else:
								myWriter.write("\t-")
								myWriter2.write("\t-")
								#write removed haplotype by COV
								myRemovedW1.write("%s\t%s\t%s\t%s\t%i\n"%(curSample, curCont, curLoc, splitted[i], readCountDICT[curSample][curCont]))


					myWriter.write("\n")
					myWriter2.write("\n")
		

		myWriter.close()
		myWriter2.close()
		myRemovedW1.close()
		myRemovedW2.close()


##################################################
####### FILTER:  populations.haplotypes.tsv ######
##################################################
print("5) Filtering new populations.cov_filtered_noNs.haplotypes.tsv")

myReader=open(args.outDir+"/populations.cov_filtered_noNs.haplotypes.tsv")
#stores the frequency of all haplotypes at a specific locus
locusHaploCOUNTS={}
locusDICT={}

for line in myReader:
	#print(line)
	#consensusLocus=False

	if line.startswith("#"):
		myHeader=line.strip().split("\t")
	else:
		splitted=line.strip().split("\t")

		curLoc=int(splitted[0])



		haploDICT={}

		#other lines can be ignored --> no variation at those loci if loci was not in populations.sumstats.tsv
		if contigID_locusID_DICT.get(curLoc, "nope") != "nope":
			curCont=contigID_locusID_DICT[curLoc]

			# counts the number of samples for each pop which cover this locus
			popCovDICT={}
			#initialize with count 0 for each population
			for pop in popDICT:
				popCovDICT[pop]=0

			for i in range(2, len(splitted)):

				curSample=myHeader[i]
				curHaplo=splitted[i]
				curPop=sampleDICT[curSample]

				#count for each pop how many samples are covering this locus
				if not curHaplo=="-" and not curHaplo=="consensus":
					popCovDICT[curPop]+=1


				if haploDICT.get(curSample, "nope") == "nope":
					haploDICT[curSample] = curHaplo
				else:
					print("Weird sample %s already in haploDICT for this locus %s!!!"%(curSample, curLoc))

			#print(popCovDICT)
			#check if covered at least by 1 indiv of args.popCOV populations
			popCoverage=0
			#gets set to false if intraCOV is violated for a species
			intraCOVfilter=True
			for pop in popCovDICT:
				# count popCOV
				if popCovDICT[pop] > 0:
					popCoverage+=1
				#filter intraCOV
				if not args.intraCOVtotal:
					if popCovDICT[pop] < args.intraCOV*popDICT_popSize[pop]:
						intraCOVfilter=False
						#print("Locus %s removed violated intraCOV for pop %s covered by %i/%f total: %i"%(curLoc, pop, popCovDICT[pop], args.intraCOV*popDICT_popSize[pop], popDICT_popSize[pop]))
				else:
					if popCovDICT[pop] < args.intraCOV:
						intraCOVfilter=False
						#print("Locus %s removed violated intraCOV for pop %s covered by %i/%f total: %i"%(curLoc, pop, popCovDICT[pop], args.intraCOV, popDICT_popSize[pop]))


			#if popCOV fulfilled
			if popCoverage>=args.popCOV and intraCOVfilter==True:
				print(popCoverage, args.popCOV)

				if locusDICT.get(curLoc, "nope") == "nope":
					print(curLoc)
					locusDICT[curLoc]=haploDICT
				else:
					print("Weird locus %s is already in the locusDICT!!!"%(curLoc))
			else:
				print("Locus %s removed"%(curLoc))

myReader.close()


######################################################
####### WRITE: populations.haplotypes.structure ######
######################################################
print("6) Writing new populations.filtered.haplotypes.structure")
#1. sampleID --> 2. key locus --> structre value at this term
#print(locusDICT)
#create int IDs for all the different haplotypes at a locus
locusDICT_Haplo_str={}
for loc in locusDICT:
	haploCounter=1
	locHaplos={}

	for sample in locusDICT[loc]:
		if locusDICT[loc][sample] != "-" and locusDICT[loc][sample] != "consensus":
			haplo1, haplo2 = locusDICT[loc][sample].split("/")
			if locHaplos.get(haplo1, "nope") == "nope":
				locHaplos[haplo1] = haploCounter
				haploCounter+=1
			if locHaplos.get(haplo2, "nope") == "nope":
				locHaplos[haplo2] = haploCounter
				haploCounter+=1
	locusDICT_Haplo_str[loc]=locHaplos



structureWriter=open(args.outDir+"/populations.filtered.haplotypes_new.structure", "w")


structureWriter.write("\t")
for key in sorted(locusDICT.keys()):
	curCONT=contigID_locusID_DICT[key]
	structureWriter.write("\t%s"%(key))
structureWriter.write("\n")

mafViolations1=0
mafViolations2=0
#loop over samples



#only those that cover a locus properly
finalDICT={}
for sID in sampleDICT:

	#write haplo for 1. allele
	#loop over loci
	tmp=0
	for key in sorted(locusDICT.keys()):
		curCONT=contigID_locusID_DICT[key]
		if tmp==0:
			structureWriter.write("%s\t%s"%(sID, popDICT_ID[sampleDICT[sID]]))
			tmp=1

		if locusDICT[key][sID] == "-" or locusDICT[key][sID] == "consensus":
			structureWriter.write("\t%i"%(-9))

		else:
			slHaplo1 = locusDICT[key][sID].split("/")[0]
			slHaplo1_ID = locusDICT_Haplo_str[key][slHaplo1]
			structureWriter.write("\t%i"%(slHaplo1_ID))
			if finalDICT.get(key, "nope")=="nope":
				finalDICT[key]={}
				finalDICT[key][sID]=1
			else:
				finalDICT[key][sID]=1


	if tmp==1:
		structureWriter.write("\n")


	
	#write haplo for 2. allele
	#loop over loci
	tmp=0
	for key in sorted(locusDICT.keys()):
		curCONT=contigID_locusID_DICT[key]
		if tmp==0:
			structureWriter.write("%s\t%s"%(sID, popDICT_ID[sampleDICT[sID]]))
			tmp=1
			
		if locusDICT[key][sID] == "-" or locusDICT[key][sID] == "consensus":

			structureWriter.write("\t%i"%(-9))

		else:		
			slHaplo2 = locusDICT[key][sID].split("/")[1]			
			slHaplo2_ID = locusDICT_Haplo_str[key][slHaplo2]

			structureWriter.write("\t%i"%(slHaplo2_ID))
			if finalDICT.get(key, "nope")=="nope":
				finalDICT[key]={}
				finalDICT[key][sID]=1
			else:
				finalDICT[key][sID]=1

	if tmp==1:
		structureWriter.write("\n")


structureWriter.close()


######################################################
####### FILTER: populations.samples.fa ###############
######################################################
print("6) Filtering the fasta file")
for file in os.listdir(args.popOutDir):

	if file.endswith("populations.samples.fa"):
		myReader=open(args.popOutDir+"/"+file)
		myWriter=open(args.outDir+"/"+"populations.samples_filtered.fa", "w")

		write=False
		for line in myReader:
			if line.startswith(">"):
				sID=line.strip().split("[")[1].split(";")[0]
				sLoc=int(line.strip().split("_")[1])
				if finalDICT.get(sLoc, "nope") != "nope":
					if finalDICT[sLoc].get(sID, "nope")==1:
						write=True
						myWriter.write(line)
			else:
				if write==True:
					myWriter.write(line)
					write=False

		myReader.close()
		myWriter.close()
