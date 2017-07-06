import os
import sys
import re
import io
import array

infile = open(sys.argv[1],'r')
ofiledir = sys.argv[3]

aggregatefile = open(sys.argv[4],'a')
fn = sys.argv[1]

pos1 = fn[-9]
pos2 = fn[-8]
pos3 = fn[-6]
pos4 = fn[-5]

arrayA = [0 for x in range(20000)]
arrayG = [0 for x in range(20000)]
arrayC = [0 for x in range(20000)]
arrayT = [0 for x in range(20000)]
position = [0 for x in range(20000)]
stoprate = 0 
total = 0
mismatch = 0
stoperror = 0 
density = [0 for x in range(20000)]

header = ''

with infile as f:
    header  = infile.readline()
    header = header.rstrip()
    if '???' in header:
        header = header[:-3]
    outfile = open(sys.argv[2]+ header + '.txt','w')
    ofile = open(ofiledir + header + '.txt','w')
    ladder1 = infile.readline()
    ladder2 = infile.readline()
    seed = infile.readline()
    seed = seed.rstrip()
    seed = seed.lstrip()
    #print(header)
    seqs = f.read()
    seqout = io.StringIO(seqs)
    
   
    for line in seqout:
        space = ''
        if 'total' in line:
            break
        whitespace = len(line)-len(line.lstrip())
        hit = re.findall("(\S+)", line)
        count = hit[1] 
        count = count[1:]
        count = int(count)
        seq = hit[0]
        seq = seq.upper()
        seq = seq.lstrip()

        for i in range(whitespace):
            space = space + ' '

        newline = space+seq
      

        if newline[10] == pos1 and newline[11] == pos2 and newline[13] == pos3 and newline[14] == pos4:
            for n in range(len(seq)):
                if seq[n] == 'A':
                    arrayA[n+whitespace] = arrayA[n+whitespace] + count
                elif seq[n] == 'C':
                    arrayC[n+whitespace] = arrayC[n+whitespace] + count
                elif seq[n] == 'G':
                    arrayG[n+whitespace] = arrayG[n+whitespace] + count
                elif seq[n] == 'T':
                    arrayT[n+whitespace]=arrayT[n+whitespace] + count
                else:
                    continue
    seqout.seek(0)
    with seqout as infile:
        #header = infile.readline()
        #header = header.rstrip()
        #if '???' in header:
                #header = header[:-3]
        #outfile = open(sys.argv[2]+ header + '.txt','w')
        #ladder1 = infile.readline()
        #ladder2 = infile.readline()
        #seed = infile.readline()
        trigger = 'total:'
        trigger2 = ' x'
        for line in infile:
            space = ''
            if trigger in line:
                totalnum = re.findall(' x: (\d+)', line, flags=0)
                #print(totalnum)
                #totalnum = int(totalnum[0])
                #total = total + totalnum
                #print('here')
                #print(total)
            elif trigger2 in line:
                whitespace = len(line) - len(line.lstrip())
                hit = re.findall("(\S+)", line)
                #for n in hit:
                #    print(n)
                #readlen = len(line)-len(line.lstrip(' '))-len(line.rstrip())
                readlen = len(hit[0])
                count = hit[1] 
                count = count[1:]
                count = int(count)
                #print(readlen)
                #position = whitespace + 1

                seq = hit[0]
                seq = seq.upper()
                seq = seq.lstrip()

                

                for i in range(whitespace):
                    space = space + ' '
					
				newline = space+seq	

                #This increments the array at a given position for the count of a read. This forms an array with values where reads have occupied along a
                #plot of tRNA position.
                if newline[10] == pos1 and newline[11] == pos2 and newline[13] == pos3 and newline[14] == pos4:    
                    for x in range(count):    
                        position[whitespace] = position[whitespace]+1
                        for n in range(whitespace, readlen + whitespace):
                            density[n]=density[n]+1

sum = float(arrayA[0]+arrayC[0]+arrayG[0]+arrayT[0]+position[1])
densitysum = float(arrayA[0]+arrayC[0]+arrayG[0]+arrayT[0])


if sum > 0:
    floatA = float(arrayA[0]/sum)
    floatC = float(arrayC[0]/sum)
    floatG = float(arrayG[0]/sum)
    floatT = float(arrayT[0]/sum)
    stoprate = float(position[1]/sum)
    #print("boop")
    if (arrayA[0]+arrayC[0]+arrayG[0]+arrayT[0]) == 0:
                  mismatch = 0
                  mutation = 0
                  stoprate = 0
    elif seed[0] is 'A':
        mismatch = float(1 - floatA)
        mutation = float(1 - floatA - stoprate)
    elif seed[0] is 'C':
        mismatch = float(1 - floatC)
        mutation = float(1 - floatC - stoprate)
    elif seed[0] is 'G':
        mismatch = float(1 - floatG)
        mutation = float(1 - floatG - stoprate)
    else:
        mismatch = float(1 - floatT)
        mutation = float(1 - floatT - stoprate)
         
    #print(str(1) + '\t' + seed[0] + '\t'+ str(arrayA[0]) + '\t' + str(arrayC[0]) + '\t' + str(arrayG[0]) + '\t' + str(arrayT[0]) + '\t' + '0' + '\t' + str(mismatch) + '\t')
    #print(str(mismatch))
    outfile.write(str(1) + '\t' + seed[0] + '\t'+ str(arrayA[0]) + '\t' + str(arrayC[0]) + '\t' + str(arrayG[0]) + '\t' + str(arrayT[0]) + '\t' + str(position[1]) + '\t' + str(mismatch) + '\t' + str(mutation) + '\t' + str(stoprate) + '\t' + str(densitysum)+ '\n')

else:
    outfile.write(str(1) + '\t' + seed[0] + '\t'+ str(arrayA[0]) + '\t' + str(arrayC[0]) + '\t' + str(arrayG[0]) + '\t' + str(arrayT[0]) + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t'+ '0' + '\t' + '0' + '\n')

for i in range(1,len(seed),1):
    
    sum = float(arrayA[i]+arrayC[i]+arrayG[i]+arrayT[i]+position[i+1])
	densum = int(arrayA[i]+arrayC[i]+arrayG[i]+arrayT[i])
    
    if sum > 0:
        floatA = float(arrayA[i]/sum)
        floatC = float(arrayC[i]/sum)
        floatG = float(arrayG[i]/sum)
        floatT = float(arrayT[i]/sum)
        stoprate = float(position[i+1]/sum)
        #stoperror = mismatch - stoprate
        if (arrayA[i]+arrayC[i]+arrayG[i]+arrayT[i]) == 0:
            mismatch = 0
        elif seed[i] is 'A':
            mismatch = float(1 - floatA)
            mutation = float(1 - floatA - stoprate)
        elif seed[i] is 'C':
            mismatch = float(1 - floatC)
            mutation = float(1 - floatC - stoprate)
        elif seed[i] is 'G':
            mismatch = float(1 - floatG)
            mutation = float(1 - floatG - stoprate)
        else:
            mismatch = float(1 - floatT)
            mutation = float(1 - floatT - stoprate)
        #print(str(mismatch))
        #print(str(i+1) + '\t' + seed[i] + '\t'+ str(arrayA[i]) + '\t' + str(arrayC[i]) + '\t' + str(arrayG[i]) + '\t' + str(arrayT[i]) + '\t' + str(position[i-1]) + '\t' + str(mismatch) + '\t')
        stoperror = mismatch - stoprate

        outfile.write(str(i+1) + '\t' + seed[i] + '\t'+ str(arrayA[i]) + '\t' + str(arrayC[i]) + '\t' + str(arrayG[i]) + '\t' + str(arrayT[i]) + '\t' + str(position[i+1]) + '\t' + str(mismatch) + '\t' + str(mutation) + '\t' + str(stoprate) + '\t'+ str(densum)+'\n')

        if i == 12:

            aggregatefile.write(str(header) + '\t' + str(floatA) + '\t' + str(floatC) + '\t' + str(floatG) + '\t' + str(floatT) + '\t' + str(mutation) + '\n')             
        
        #if mismatch > .20 and mismatch < 1.0 and sum > 25 and stoperror > .15 and stoprate > .10:
        #print(header + " " + "at position " + str(i+1)+ " with identity " +seed[i] + " has a mismatch rate of " + str(mismatch) + " which includes a stop rate of " + str(stoprate) + " and has been counted " + str(sum)+" times") 
        #print(stoperror)
        #ofile.write(header + " " + "at position " + str(i+1)+ " with identity " +seed[i] + " has a mismatch rate of " + str(mismatch) + " which includes a stop rate of " + str(stoprate) + " and has been counted " + str(sum)+" times" + '\n')
    else:
        outfile.write(str(i+1) + '\t' + seed[i] + '\t'+ str(arrayA[i]) + '\t' + str(arrayC[i]) + '\t' + str(arrayG[i]) + '\t' + str(arrayT[i]) + '\t' + str(position[i+1]) + '\t' + '0' + '\t' + '0' + '\t'+ '0' + '\t'+ str(density[i]) + '\n')

outfile.close()
