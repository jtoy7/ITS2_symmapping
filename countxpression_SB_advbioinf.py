#!/usr/bin/env python
#----------------USAGE-----------------------
#by Dan Barshis
#edited by Sandrine Boissel 6/17 to deal with multiply matching reads, collapse isoform counts to genes and use input flags
#This script parses a .sam alignment file created using bowtie2 and returns a tab-delimited text file for every input file with the 
#number of good quality multiply and singly aligned reads (above a certain mapping quality threshold and read length or alignment score) 
#as well as a single file with summary stats for all input files. This script can also take a tab-delimited list of all isogroup and genes
#in the format isogroup\tgene to collapse counts to the gene level. 
#The command line requires 1 flag: 	-i anynumberofinputfiles
#and has a number of optional flags:  	-o outputstatsfilename 	(defaults to match_counts.txt)
#					-t threshold		(defaults to 20)
#					-l lengththreshold 	(defaults to 20)
#					-a ASthreshold		(defaults to 40)
#					-g isogroup_table	(defaults to False)
#as such with optional flags in paranthesis:
#countexpression_SB.py -i *.sam (-t N -l N -a N -o output.txt -g table.txt)  

def countxpression(infilename, threshold, lengththresh, asthreshold, zeroforheader1forno, outfilename, countsoutputname,isogroup_table):
        #Opens an infile specified by the user. Should be a .sam file 
        IN = open(infilename, 'r')

        threshold=float(threshold) #inputs a mapping quality threshold for counting reads as sucessfully mapped
        lenthresh=float(lengththresh)
        ASthreshold=float(asthreshold)

        #Opens an output text file as specified by user 
        OUT = open(outfilename, 'w')

        countsout = open(countsoutputname, 'a')

        linecount=0
        totreads=0
        notaligned=0
        multi=0
        single=0
        aligned=0
        goodmaps=0
        contigswzeros=0
        #Starts a for loop to read through the infile line by line
        multicontigs={}
        contigs={}
        for line in IN:
                line=line.rstrip()
                cols=line.split('\t')

                if cols[0]=="@SQ": #Generate dictionary of contig names
                        contigname=cols[1].split(':')[1]
                        contigs[contigname]=[0, 0, 0]
                if cols[0][0] != '@' : #to skip past header
                        columntocount=13
			if cols[2] == '*': #this is to count reads that are flagged as unaligned
				notaligned+=1
				totreads+=1
			else:
				if cols[12].startswith('XS'): #this is to find reads that had >1 alignment
					aligned+=1
					totreads+=1
					multi+=1
					AS=int(cols[11].split(':')[2])
					if AS >= ASthreshold:
						goodmaps+=1
						contigs[cols[2]][1]+=1
				else:	#find reads aligned only once
					aligned+=1
					totreads+=1
					single+=1
					if float(cols[4]) >= threshold and len(cols[9]) >= lenthresh:
						goodmaps+=1
						contigs[cols[2]][0]+=1
        if isogroup_table!=False:
		iso_to_gene={}
		i2g_file=open(isogroup_table,'r')
		for line in i2g_file:
			cols=line.strip().split()
			iso_to_gene[cols[0]]=cols[1]
		i2g_file.close()
		genes={}
	for item in contigs:
                contigs[item][2]=contigs[item][0]+contigs[item][1]
                if contigs[item][0] == 0:
			contigswzeros+=1
                if isogroup_table!=False:
			contig_gene=iso_to_gene[item]
			if contig_gene in genes:
				genes[contig_gene][0]+=contigs[item][0]
				genes[contig_gene][1]+=contigs[item][1]
				genes[contig_gene][2]+=contigs[item][2]
			else:
				genes[contig_gene]=contigs[item]                
        
        contigsmatched=len(contigs.keys())-contigswzeros
        propqualaligned=float(goodmaps)/float(totreads)

        if int(zeroforheader1forno) == 0:
                countsout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Filename', 'Total#Reads', 'NumContigsMatched', 'NumUnaligned', 'NumAligned', 'NumMultiAligned', 'NumSingleAligned', 'NumQualSingles', 'PropQualAligned'))
        countsout.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n' % (infilename, totreads, contigsmatched, notaligned, aligned, multi, single, goodmaps, propqualaligned))


        OUT.write('ContigName'+'\t'+'UniqueTotReads'+'\t'+'MultiTotReads'+'\t'+'totalreadsgoodmapped\n')
        l=[]

#this is for general sorting based on text sorting rules
        if isogroup_table==False:
		for key,value in contigs.items():
                	l.append((key,value))
	else:
		for key,value in genes.items():
                        l.append((key,value))
        l.sort()
        for item in l:
        #       print item
                OUT.write('%s\t%d\t%d\t%d\n' % (item[0], item[1][0], item[1][1], item[1][2])) #writes each line of the tuple as separate tab delimited text

        OUT.close()

if __name__=="__main__":
	import argparse
	parser = argparse.ArgumentParser(description='This script parses a .sam alignment file created using bowtie2 and returns a tab-delimited text file for every input file with the number of good quality multiply and singly aligned reads (above a certain mapping quality threshold and read length or alignment score) as well as a single file with summary stats for all input files. This script can also take a tab-delimited list o all isogroup and genes in the format isogroup\\tgene to collapse counts to the gene level.')
	parser.add_argument('-t', action="store", dest='threshold', default=20, type=int, help="Give a MAPQ cutoff, defaults to 20")
	parser.add_argument('-l', action="store", dest='lengththreshold', default=20, type=int, help="Give a match length cutoff, defaults to 20")
	parser.add_argument('-a', action="store", dest='ASthreshold', default=40, type=int, help="Give an alignment score cutoff, defaults to 40")
	parser.add_argument('infiles', nargs='+', help="Give all .sam input files")
	parser.add_argument('-o', action="store", dest='outfile', default='match_counts.txt', help="Give an output file name, defaults to match_counts.txt")
	parser.add_argument('-g', action="store", dest="isogroup_table", default=False, help="Give a file with an isogroup-gene name table")
	args=parser.parse_args()
	numfiles=0
	for name in args.infiles:
		numfiles+=1
                if numfiles==1:
        		countxpression(name, args.threshold, args.lengththreshold, args.ASthreshold, 0, name[:-4]+'_counts.txt',args.outfile,args.isogroup_table)
		if numfiles>1:
			countxpression(name, args.threshold, args.lengththreshold, args.ASthreshold, 1, name[:-4]+'_counts.txt',args.outfile,args.isogroup_table)

