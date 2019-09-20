"""
Author: Muhammad Zohaib Anwar
License: GPL v3.0\n\n

"""

import subprocess
import sys
import argparse
import os
import os.path as path
from Bio import SeqIO
import os
import numpy

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-ref", "--reffile", help= "Reference Genome or Assembled contig file", required=True)
parser.add_argument("-sam", "--inputfile", help= "Input Alignment in SAM format", required=True)
parser.add_argument("-o", "--outputfile", help= "Output matrix .tsv format", required=True)
parser.add_argument("-c", "--CIGAR", help= "75M, 101M, 151M, 251M (etc)", required=True)
parser.add_argument("-r", "--readtype", help= "single-read (se) or paired end reads (pe)", required=True)
parser.add_argument("-l", "--length", help= "Maximum length of reads", required=True, type = int, default = 75)
parser.add_argument("-covDNA", "--covfileDNA", help= "Coverage of hotspot regions in DNA file", required= True)
parser.add_argument("-covcDNA", "--covfilecDNA", help= "Coverage of hotspot regions in cDNA file", required= True)
parser.add_argument("-cov", "--cov", help= "Minimum coverage after CPM Normalization", required=True, type = float, default = 1.0)
parser.add_argument("-maxcov", "--maxcov", help = "Maximum coverage after CPM normalization", required= True, type = float, default = 100)

args = parser.parse_args()
if not args.reffile: print("No reference file provided")
if not args.inputfile: print("No input file provided")
if not args.outputfile: print("No output file provided")
if not args.CIGAR: print("No CIGAR provided")
if not args.readtype: print("Please provide read type - se or pe")
if not args.covfileDNA: print("Please provie a coverage file")
if not args.covfilecDNA: print("Please provie a coverage file")

if not (args.inputfile or args.outputfile or CIGAR or length or readtype or covfileDNA or covfilecDNA): sys.exit(1)

def fasta2dict(filename):
	dic = {}
	cur_scaf = ''
	cur_seq = []
	infile=open(str(filename),'r')
	for line in open (str(filename), mode='r'):
		if line.startswith(">") and cur_scaf == '':
			cur_scaf = line.split(' ')[0]
			cur_scaf = cur_scaf.strip()
			cur_scaf = cur_scaf.replace(">","")
		elif line.startswith(">") and cur_scaf != '':
			dic[cur_scaf] = ''.join(cur_seq)
			cur_scaf = line.split(' ')[0]
			cur_scaf = cur_scaf.strip()
			cur_scaf = cur_scaf.replace(">","")
			cur_seq = []
		else:
			cur_seq.append(line.rstrip())
		dic[cur_scaf] = ''.join(cur_seq)
	return dic

def calgc(seq):
	a = seq.count('A')
	c = seq.count('C')
	g = seq.count('G')
	t = seq.count('T')
	gc= float(c + g) / (a+c+g+t)
	return gc

def check_error(query,template, orientation, error_dic, gc_error):
	gc=calgc(seq=template)
	gc*=100
	if 0.0 <= gc <= 5.0:
		gc_row=0
	if 5.01 <= gc <= 10.0:
		gc_row=1
	if 10.01 <= gc <= 15.0:
		gc_row=2
	if 15.01 <= gc <= 20.0:
		gc_row=3
	if 20.01 <= gc <= 25.0:
		gc_row=4
	if 25.01 <= gc <= 30.0:
		gc_row=5
	if 30.0 <= gc <= 35.0:
		gc_row=6
	if 35.01 <= gc <= 40.0:
		gc_row=7
	if 40.0 <= gc <= 45.0:
		gc_row=8
	if 45.01 <= gc <= 50.0:
		gc_row=9
	if 50.0 <= gc <= 55.0:
		gc_row=10	
	if 55.01 <= gc <= 60.0:
		gc_row=11
	if 60.0 <= gc <= 65.0:
		gc_row=12
	if 65.01 <= gc <= 70.0:
		gc_row=13
	if 70.0 <= gc <= 75.0:
		gc_row=14
	if 75.01 <= gc <= 80.0:
		gc_row=15
	if 80.0 <= gc <= 85.0:
		gc_row=16
	if 85.01 <= gc <= 90.0:
		gc_row=17
	if 90.01 <= gc <= 95.0:
		gc_row=18
	if 95.01 <= gc <= 100.0:
		gc_row=19


	for i in range(0,len(query)):
		if query[i] != template[i]:
			if orientation == 0:
				index = int(i)
			else:
				index = len(query)-(i+1)
			target = 0
			ref = template[i]
			alternate = query[i]

			if ref == 'A' and alternate == 'T':
				target = 1
				error_dic["A-T"] +=1
			elif ref == 'A' and alternate == 'G':
				target = 2
				error_dic["A-G"] +=1
			elif ref == 'A' and alternate == 'C':
				target = 3
				error_dic["A-C"] +=1
			elif ref == 'A' and alternate == 'N':
				target = 4
				error_dic["A-N"] +=1
			elif ref == 'T' and alternate == 'A':
				target = 5
				error_dic["T-A"] +=1
			elif ref == 'T' and alternate == 'G':
				target = 6
				error_dic["T-G"] +=1
			elif ref == 'T' and alternate == 'C':
				target = 7
				error_dic["T-C"] +=1
			elif ref == 'T' and alternate == 'N':
				target = 8
				error_dic["T-N"] +=1
			elif ref == 'G' and alternate == 'A':
				target = 9
				error_dic["G-A"] +=1
			elif ref == 'G' and alternate == 'T':
				target = 10
				error_dic["G-T"] +=1
			elif ref == 'G' and alternate == 'C':
				target = 11
				error_dic["G-C"] +=1
			elif ref == 'G' and alternate == 'N':
				target = 12
				error_dic["G-N"] +=1
			elif ref == 'C' and alternate == 'A':
				target = 13
				error_dic["C-A"] +=1
			elif ref == 'C' and alternate == 'T':
				target = 14
				error_dic["C-T"] +=1
			elif ref == 'C' and alternate == 'G':
				target = 15
				error_dic["C-G"] +=1
			elif ref == 'C' and alternate == 'N':
				target = 16
				error_dic["C-N"] +=1
			elif ref == 'N' and alternate == 'A':
				target = 17
				error_dic["N-A"] +=1
			elif ref == 'N' and alternate == 'T':
				target = 18
				error_dic["N-T"] +=1
			elif ref == 'N' and alternate == 'G':
				target = 19
				error_dic["N-G"] +=1
			elif ref == 'N' and alternate == 'C':
				target = 20
				error_dic["N-C"] +=1

			else:
				print(index, ref, alternate)
			error[index,target] = error[index,target] +1
			gc_error[gc_row, target] =gc_error[gc_row, target] +1

#def coverage(DNA_file, cDNA_file, coverage):
def coverage(DNA_file, cDNA_file, coverage, max_coverage):

	f = open(DNA_file, 'r')
	lines=f.readlines()
	f.close()

	f = open(cDNA_file, 'r')
	lines1=f.readlines()
	f.close()

	DNA_nodes={}
	cDNA_nodes={}

	for line in lines:
		line=line.strip()
		if line.split("\t")[0] not in DNA_nodes.keys() and float(line.split("\t")[3]) >= coverage and float(line.split("\t")[3]) <= max_coverage:
		#if line.split("\t")[0] not in DNA_nodes.keys():
			DNA_nodes[line.split("\t")[0]]=line.split("\t")[1]+","+line.split("\t")[2]
		elif line.split("\t")[0] in DNA_nodes.keys() and float(line.split("\t")[3]) >= coverage and float(line.split("\t")[3]) <= max_coverage:
		#else:
			DNA_nodes[line.split("\t")[0]]=DNA_nodes[line.split("\t")[0]]+";"+line.split("\t")[1]+","+line.split("\t")[2]

	for line in lines1:
		line=line.strip()
		if line.split("\t")[0] not in cDNA_nodes.keys() and float(line.split("\t")[3]) >= coverage:
			cDNA_nodes[line.split("\t")[0]]=line.split("\t")[1]+","+line.split("\t")[2]
		elif line.split("\t")[0] in cDNA_nodes.keys() and float(line.split("\t")[3]) >= coverage:
			cDNA_nodes[line.split("\t")[0]]=cDNA_nodes[line.split("\t")[0]]+";"+line.split("\t")[1]+","+line.split("\t")[2]
	
	combined = {}
	for x in DNA_nodes:
		if x in cDNA_nodes.keys():
			i=DNA_nodes[x].split(";")
			j=cDNA_nodes[x].split(";")
			#combined[x]=sorted(list(set(i).intersection(j)))
			combined[x]=sorted(list(set(i).union(j)))		
	
	return combined

def checkcov(read,pos,bins):
	result = ""
	if read in bins:
		#print(read)
		#print("yes read in bin",str(bins[read]))
		for x in bins[read]:
			#print(bins[read])
			if float(x.split(",")[0]) <= pos <= float(x.split(",")[1]):
				result = "yes"
				break
			else:
				result="no"   
	return result

def normalize(filename, total):
	f = open(filename, 'r')
	lines = f.readlines()
	f.close()

	if "gc" in filename:
		f = open(filename.replace(".tsv","_normalized.tsv"),'w')
		f.write(lines[0])
		for i in range(1,len(lines)):
			lines[i] = lines[i].strip()
			cells = lines[i].split("\t")
			if int(cells[0]) < 10:
				cells[0] = "Bin00"+str(cells[0])
			else:
				cells[0] = "Bin0"+str(cells[0])
			s = cells[0]
			for i in range(1,len(cells)):
				x = 0.0
				x = (float(cells[i])/(total*75)) * 1000000
				s = s +"\t"+ str(x)
			f.write(s+"\n")
		f.close()
	else:
		f = open(filename.replace(".tsv","_normalized.tsv"),'w')
		f.write(lines[0])
		for i in range(1,len(lines)):
			lines[i] = lines[i].strip()
			cells = lines[i].split("\t")
			if int(cells[0]) <10:
				cells[0] = "P00"+str(cells[0])
			else:
				cells[0] = "P0"+str(cells[0])
			s = cells[0]
			for i in range(1,len(cells)):
				x = 0.0
				x = (float(cells[i])/(total*75)) * 1000000
				s = s +"\t"+ str(x)
			f.write(s+"\n")
		f.close()


Scriptdir = os.path.realpath(__file__)
refdir, ref_file = os.path.split(args.reffile)
outputdir, outputfile = os.path.split(args.outputfile)
inputdir, inputfile = os.path.split(args.inputfile)
CIG = args.CIGAR
read_length = args.length
sequence_type = args.readtype
DNAcovdir, DNAcovfile = os.path.split(args.covfileDNA)
cDNAcovdir, cDNAcovfile = os.path.split(args.covfilecDNA)
error = numpy.zeros((int(read_length),21,), dtype=int)
gc_error = numpy.zeros((20,21,), dtype=int)
c = args.cov
max_c = args.maxcov

if __name__ == "__main__":

	import csv
	coveragedic=coverage(DNA_file=DNAcovdir +"/"+ DNAcovfile, cDNA_file= cDNAcovdir +"/" + cDNAcovfile, coverage=c, max_coverage=max_c)
	#print(coveragedic)
	Contig_dic = fasta2dict(filename = refdir + "/"+ ref_file)
	print("Fasta file read and Contig Dictionary made\n\n")
	error_dic={
		"A-T" : 0, "A-G" : 0, "A-C" : 0, "A-N" : 0,
		"T-A" : 0, "T-G" : 0, "T-C" : 0, "T-N" : 0,
		"G-A" : 0, "G-T" : 0, "G-C" : 0, "G-N" : 0,
		"C-A" : 0, "C-T" : 0, "C-G" : 0, "C-N" : 0,
		"N-A" : 0, "N-T" : 0, "N-G" : 0, "N-C" : 0
		}	

	print("Reading the Sam file now\n\n")

	f = open(inputdir + "/" + inputfile ,"r")
	lines = f.readlines()
	f.close()

	print("Sam file Read\n\n")
	
	match = 0
	mismatch = 0
	FLAG=[]
	CIGAR= ""
	for i in range(0,read_length):
		error[i,0]=i
	for i in range(0,20):
		gc_error[i,0]=i	
	for line in lines:
		line = line.strip()
		if not line.startswith("@"):
			CIGAR = line.split("\t")[5].strip()
			#print(CIGAR)
			if CIGAR == CIG:
				attributes= line.split("\t")
				qreadID = attributes[0]
				Flag = attributes[1]
				temp_read = Contig_dic[attributes[2]]
				select=checkcov(read = str(attributes[2]), pos=int(attributes[3]), bins=coveragedic)
				if select=='yes':
					#print("yes")
					temp_read = str(temp_read[int(attributes[3])-1:int(attributes[3])+(read_length-1)])
					query_read = str(attributes[9])
					if sequence_type == "se":
						if Flag in ['0','1024'] and not query_read == temp_read:
							check_error(query = query_read, template = temp_read, orientation = 0,error_dic=error_dic, gc_error=gc_error)
							mismatch = mismatch + 1
						if Flag in ['16','1040'] and not query_read == temp_read:
							check_error(query = query_read, template = temp_read, orientation = 1, error_dic=error_dic, gc_error=gc_error)
							mismatch = mismatch + 1
						if Flag in ['0', '16', '1024', '1040'] and query_read == temp_read:
							match = match + 1

					elif sequence_type == "pe":
						if Flag in ['99'] and not query_read == temp_read:
							check_error(query = query_read, template = temp_read, orientation = 0, error_dic=error_dic, gc_error=gc_error)
							mismatch = mismatch +1
						if Flag in ['147'] and not query_read == temp_read:
							check_error(query = query_read, template = temp_read, orientation = 1, error_dic=error_dic, gc_error=gc_error)
							mismatch = mismatch + 1
						else:
							match = match +1


	print("Matched Reads = \t"+str(match))
	print("error Reads = \t"+str(mismatch))
	Total = match+mismatch
	print(Total)

	f= open(outputdir+"/"+outputfile , 'w')
	f.write("Position\tA_T\tA_G\tA_C\tA_N\tT_A\tT_G\tT_C\tT_N\tG_A\tG_T\tG_C\tG_N\tC_A\tC_T\tC_G\tC_N\tN_A\tN_T\tN_G\tN_C\n")
	f.close()

	with open(outputdir+"/"+outputfile , 'ab') as abc:
		numpy.savetxt(abc, error,fmt="%s", delimiter="\t")

	normalize(filename = outputdir+"/"+outputfile, total = int(Total))

	f= open(outputdir+"/"+outputfile.replace(".tsv", "_gc.tsv") , 'w')
	f.write("GC_bin\tA_T\tA_G\tA_C\tA_N\tT_A\tT_G\tT_C\tT_N\tG_A\tG_T\tG_C\tG_N\tC_A\tC_T\tC_G\tC_N\tN_A\tN_T\tN_G\tN_C\n")
	f.close()

	with open(outputdir+"/"+outputfile.replace(".tsv", "_gc.tsv") , 'ab') as abc:
		numpy.savetxt(abc, gc_error,fmt="%s", delimiter="\t")

	normalize(filename = outputdir+"/"+outputfile.replace(".tsv", "_gc.tsv"), total = int(Total))
		
	f = open(outputdir+"/"+outputfile.replace(".tsv", "_errortype.tsv") ,'w')
	f.write("ErrorType\tEPM\n")
	for key in error_dic.keys():
		f.write(str(key)+"\t"+str(float(error_dic[key]/Total)*1000000)+"\n")
	f.close()

	#keys = ['Error_Type','value']
	#csvfilename = os.path.join(outputdir+"/"+outputfile.replace(".tsv", "_errortype.tsv"))
	#with open(csvfilename, 'w') as output_file:
	#	dict_writer = csv.DictWriter(output_file, keys)
	#	dict_writer.writeheader()
	#	dict_writer.writerows(error_dic)
