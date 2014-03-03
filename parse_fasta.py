#!/usr/bin/python

def load_fasta_to_dict(fh):
	dic=dict()
	temp=str()

	for line in fh:
		if '>' == line[0]:
			temp=line[1:].split()[0] #to work with mbk.goldstandard.fasta wich contains the ORF position in the header
		else:
			dic[temp]=line

	# for i in dic.keys():
	# 	print i+'\n'+dic[i],;
	return dic


####################################################################

####################################################################

if __name__ == "__main__":
	from sys import argv
	fh=open(argv[1])
	dic=load_fasta_to_dict(fh)
	intersect=open(argv[2])

	for line in intersect:
		a=line.split()
		if int(a[1])-20>0:
			b=int(a[1])-20
		else:
			b=0
		print '>'+a[0]+'\n'+dic[a[0]][b:int(a[2])+20]