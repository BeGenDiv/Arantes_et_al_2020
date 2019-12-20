import argparse
import gzip
import os


class myFQ:
	def __init__(self, header, seq, comm, qual):
		self.header=header
		self.seq=seq
		self.comm=comm
		self.qual=qual


parser = argparse.ArgumentParser(description="script to remove PCR duplicates from paired end fastq files (barcodes must be present in the fastq ID line for each sequence)")
parser.add_argument("-r1", "--read1", help="forward reads in fastq format(can be zipped)",type=str, required=True)
parser.add_argument("-r2", "--read2", help="forward reads in fastq format(can be zipped)",type=str, required=True)
parser.add_argument("-o", "--out", help="output directory",type=str, required=True)
args=parser.parse_args()


if not os.path.isdir(args.out):
	os.mkdir(args.out)


#open files
if args.read1.endswith(".gz"):
	r1Reader=gzip.open(args.read1)	
else:
	r1Reader=open(args.read1)

if args.read2.endswith(".gz"):
	r2Reader=gzip.open(args.read2)	
else:
	r2Reader=open(args.read2)


r1out=args.read1.split("/")[-1].split(".")[0]+"_PCRfiltered.fastq"
r2out=args.read2.split("/")[-1].split(".")[0]+"_PCRfiltered.fastq"

oWriter1=open(args.out+"/"+r1out, "w")
oWriter2=open(args.out+"/"+r2out, "w")

oligoDICT={}
oligoDICT_R1={}
oligoDICT_R2={}

readToOligo={}


lineCount=0
seqCount=0
dupCount=0
write=False
for line1 in r1Reader:
	line2=r2Reader.readline()

	lineCount+=1

	if lineCount==1:
		rID=line1.strip().split(" ")[0]

		cOli=line1.strip().split("+")[-1]

		myR1=myFQ(line1.strip(), "", "", "")
		myR2=myFQ(line2.strip(), "", "", "")
	
	elif lineCount==2:
		seq1=line1.strip()
		seq2=line2.strip()
		
		myR1.seq=seq1
		myR2.seq=seq2

		identSeq=cOli+seq1[:100]+seq2[:100]

		if oligoDICT.get(identSeq, "nope")=="nope":
			oligoDICT[identSeq]=1
			oligoDICT_R1[identSeq]=[]
			oligoDICT_R2[identSeq]=[]
			write=True


		else:
			oligoDICT[identSeq]+=1
			dupCount+=1


	elif lineCount==3:
		myR1.comm=line1.strip()
		myR2.comm=line2.strip()


	elif lineCount==4:
		myR1.qual=line1.strip()
		myR2.qual=line2.strip()

		oligoDICT_R1[identSeq].append(myR1)
		oligoDICT_R2[identSeq].append(myR2)

		if write==True:
			oWriter1.write(myR1.header+"\n")
			oWriter1.write(myR1.seq+"\n")
			oWriter1.write(myR1.comm+"\n")
			oWriter1.write(myR1.qual+"\n")

			oWriter2.write(myR2.header+"\n")
			oWriter2.write(myR2.seq+"\n")
			oWriter2.write(myR2.comm+"\n")
			oWriter2.write(myR2.qual+"\n")

		lineCount=0
		seqCount+=1
		write=False



r1Reader.close()
r2Reader.close()


statWriter=open(args.out+"/"+args.read1.split("/")[-1].split(".")[0].strip("_1")+"_filterPCRdups_Stats.tsv", "w")


r1ClustWriter=open(args.out+"/"+args.read1.split("/")[-1].split(".")[0]+"_PCRfiltered_CLUSTERS.fastq", "w")
r2ClustWriter=open(args.out+"/"+args.read2.split("/")[-1].split(".")[0]+"_PCRfiltered_CLUSTERS.fastq", "w")




for oligo in oligoDICT_R1:
	#print(oligo, len(oligoDICT_R1[oligo]))
	if len(oligoDICT_R1[oligo]) > 1:
		for read1 in oligoDICT_R1[oligo]:
			r1ClustWriter.write(read1.header+"\n")
			r1ClustWriter.write(read1.seq+"\n")
			r1ClustWriter.write(read1.comm+"\n")
			r1ClustWriter.write(read1.qual+"\n")
		
		for read2 in oligoDICT_R2[oligo]:
			r2ClustWriter.write(read2.header+"\n")
			r2ClustWriter.write(read2.seq+"\n")
			r2ClustWriter.write(read2.comm+"\n")
			r2ClustWriter.write(read2.qual+"\n")


		r1ClustWriter.write("######################################################\n")
		r2ClustWriter.write("######################################################\n")











print("Total number paired reads: %i"%(seqCount))
statWriter.write("Total number paired reads: %i\n"%(seqCount))
print("Total number of PCR duplicates: %i"%(dupCount))
statWriter.write("Total number of PCR duplicates: %i\n"%(dupCount))

countDICT={}

for ele in oligoDICT:
	oCounts=oligoDICT[ele]
	if countDICT.get(oCounts, "nope") == "nope":
		countDICT[oCounts]=1
	else:
		countDICT[oCounts]+=1


statWriter.write("######################################################\n")
for key in sorted(countDICT.iterkeys()):
   # print "%s: %s" % (key, countDICT[key])
    statWriter.write("%s\t%s\n"%(key, countDICT[key]))
statWriter.write("######################################################\n")
for key, value in sorted(oligoDICT.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    #print "%s: %s" % (key, value)
    statWriter.write("%s\t%s\n"%(key, value))


