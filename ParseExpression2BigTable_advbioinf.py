#!/usr/bin/env python

import sys
#sys.argv[1] Input file name #list of your "gene" names
#sys.argv[2] Output file name
#sys.argv[3] Text to add when no match
#sys.argv[4:] Any number of files to add columns from

nohittext=sys.argv[3]
columntoextract = 4

#this is for setting how the script sorts your contigs into order
#change the word to 'text' for a text-based sorting or 'coral' for a
#palumbi-lab coral-specific numerical sorting

def make_dict1(file):
	fin = open(file, 'r')
	dict={}
	headers=[]
	count=0
	for line in fin:
		count+=1
		line=line.rstrip()
		cols=line.split('\t') #for tab-delimited text files
		if count==1:
			headers=cols[0:]
			#count+=1
		if count > 1:
			dict[cols[0]]=cols[1:]
	
	return dict, headers
		
dictbase, dictheaders=make_dict1(sys.argv[1]) #Input data table
xpressionfiles=sys.argv[4:] #Any number of expression files
	
xpressfiles=[]
xtrahits=[]
for file in xpressionfiles:
	Xvalues, Xheaders=make_dict1(file)
	xpressfiles.append(file+'_'+str(Xheaders[columntoextract-1])) # used to denote the column that you're extracting
	xtracount=0
	matched = list(set(dictbase) & set(Xvalues))
	nomatch = list(set(dictbase) ^ set(Xvalues))
	in1not2 = list(set(nomatch) & set(dictbase))
	in2not1 = list(set(nomatch) & set(Xvalues))
#	For adding info when Xvalues is missing values that are in your dictbase
	if len(in1not2) >=1: 	
		for item in dictbase.keys():
			if Xvalues.has_key(item):
#				dictbase[item]+=Xvalues[item][1:3] #appends a list of strings
				dictbase[item].append(Xvalues[item][columntoextract-2])
			else:
				xtracount+=1
				dictbase[item].append(nohittext)
# 	For adding info when Xvalues has all the values in your dictbase or more
	if len(in1not2)==0 and in2not1>=1:
		for item in Xvalues.keys():
			if dictbase.has_key(item):
#				dictbase[item]+=Xvalues[item] #appends a list of strings
				dictbase[item].append(Xvalues[item][columntoextract-2]) #appends a single column of data as a string appends quality, single counts
	 		else:
				xtracount+=1
 	xtrahits.append(file+'='+str(xtracount))

o=open(str(sys.argv[2]), 'w') # New data table file name

o.write('\t'.join(dictheaders)+'\t'+'\t'.join(xpressfiles)+'\n') #used if you want a specific header and filename for each file

print 'Hits not matched' + '\t'.join(xtrahits)

############this if for sorting your contigs based on text order#############
l=[]
for key,value in dictbase.items():
	l.append((key,value))
l.sort()
for item in l:
	o.write('%s\t%s\n' % (item[0], '\t'.join(item[1]))) #writes each line of the tuple as separate tab delimited text
		
o.close()
