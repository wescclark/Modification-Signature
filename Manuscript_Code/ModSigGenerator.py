import re
import sys
import array
import os.path
import io
import xlsxwriter
from decimal import *

filelist = open(sys.argv[1],'r')
splitfile = sys.argv[2]
countfile = sys.argv[3]
#flagfile = sys.argv[4]
outfile = open(sys.argv[4],'w')

seed=''
flaggedpos=[]
midict = {}

def get_seed(file):
    with file as infile:
        header = infile.readline()
        ladder1 = infile.readline()
        ladder2 = infile.readline()
        seed = infile.readline()
        seed = seed.rstrip()
        seed = seed.lstrip()
    return seed

def get_flagged_pos(file):
    posarray=[]
    with file as infile:
        for line in infile:
            hit = re.findall("(\S+)", line)
            #print(str(hit[1]))
            posarray.append(int(hit[1]))
    return posarray

def get_MI(file):
    mismatchdict = {}
    with file as infile:
        header = infile.readline()
        for line in infile:
            if len(line)==0:
                break
            position = 0
            nt = ''
            mismatch = 0
            mutation = 0
            stop = 0
            density = 0

            mihits = re.findall("(\S+)", line)
            position = str(mihits[0])
            nt = str(mihits[1])
            mismatch = str(mihits[7])
            mutation = str(mihits[8])
            stop = str(mihits[9])
            density = str(mihits[10])
            key = position

            Acount = float(mihits[2])
            Ccount = float(mihits[3])
            Gcount = float(mihits[4])
            Tcount = float(mihits[5])
            countsum = Acount+Ccount+Gcount+Tcount

            if countsum > 0:
                Aratio = Acount/countsum
                Cratio = Ccount/countsum
                Gratio = Gcount/countsum
                Tratio = Tcount/countsum
            else:
                Aratio = 0
                Cratio = 0
                Gratio = 0
                Tratio = 0
            print(position)
            mismatchdict[position]=[nt]
            #mismatch[position].append(nt)
            mismatchdict[position].append(mismatch)
            mismatchdict[position].append(mutation)
            mismatchdict[position].append(stop)
            mismatchdict[position].append(density)

            mismatchdict[position].append(str(Aratio))
            mismatchdict[position].append(str(Cratio))
            mismatchdict[position].append(str(Gratio))
            mismatchdict[position].append(str(Tratio))
            
    
        return mismatchdict
    
outfile.write('Source'+'\t'+'Position'+'\t'+'Sequence'+'\t'+'Nucleotide'+'\t'+'MI'+'\t'+'Mutation Rate'+'\t'+'Stop Rate'+'\t'+'Density'+'\t'+'A ratio'+'\t'+'C ratio'+'\t'+'G ratio'+'\t'+'T ratio'+'\n')
for line in filelist:
    filename = line
    filename = filename.rstrip()

    splitfilename = splitfile+filename
    #flagfilename = flagfile+filename
    countfilename = countfile+filename

    splitfilenew = open(splitfilename,'r')
    #flagfilenew = open(flagfilename,'r')
    countfilenew = open(countfilename,'r')

    #print(splitfilename)
    #print(flagfilename)
    print(countfilename)

    
    seed = get_seed(splitfilenew)
    #flaggedpos = get_flagged_pos(flagfilenew)
    midict = get_MI(countfilenew)

    print(seed)
    
    for n in range(len(midict)):
        #print(n)
        #print(seed[n-3:n+2])
        #print(midict[str(n)])

        filename = filename.rstrip()
        #print(n)
        if float(midict[str(n+1)][1]) > .4 and float(midict[str(n+1)][4]) > 100: 
            outfile.write(filename + '\t' + str(n+1)+'\t'+seed[n-2:n+3]+'\t')
            for l in midict[str(n+1)]:
                outfile.write(str(l)+'\t')
            print(*midict[str(n+1)],sep=', ')
            outfile.write('\n')

    splitfilenew.close()
    countfilenew.close()
    #flagfilenew.close()

filelist.close()
