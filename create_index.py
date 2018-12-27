#!/usr/bin/python
import argparse
import sys
import os.path
import subprocess
import itertools

#def hamming2(str1, str2):
#    assert len(str1) == len(str2)
#    return sum(d1 != d2 for d1, d2 in zip(str1, str2))

def hamming2(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
           diffs += 1
    return diffs
def genome2spacerfasta(nameofgenometofasta):
    files2genome={line.strip().split('\t')[1]:line.strip().split('\t')[0] for line in nameofgenometofasta}
#    print files2genome
    return files2genome

def comparefasta(fasta2r,mismatch):
    comparisonlist=list(itertools.combinations(fasta2r,2))
    print str(len(comparisonlist))+' pairwise comparisons'
    for keys in comparisonlist:
        if hamming2(keys[0], keys[1]) <=int(mismatch) and keys[1] in fasta2r and keys[0] in fasta2r:
            if len(fasta2r[keys[0]])> len(fasta2r[keys[1]]):
		fasta2r[keys[0]]=fasta2r[keys[0]]+fasta2r[keys[1]]
		fasta2r.pop(keys[1], None)
	    else:
		fasta2r[keys[1]]=fasta2r[keys[1]]+fasta2r[keys[0]]
		fasta2r.pop(keys[0], None)
#    print fasta2r
    print '0 mismatch filt:', len(fasta2r)
    return fasta2r

def spacerfastacluster(files2genome):
    files=files2genome.keys()
    fasta2genome={}
    for f in files:
#        print f
        with open(f, 'r') as spacerfasta:
            for line in spacerfasta:
#                print line
                if '>' not in line:
                    fasta2genome.setdefault(line.strip(),[]).append(files2genome[f])
#    print fasta2genome
    print '0 mismatch:',len(fasta2genome)
    return fasta2genome

def buildindex(fasta2r, files2genome,output):
    fasta2number={}
    f1=open('consolidatedspacers.fa', 'w')
    f2=open('spacer2genomes.txt','w')
    i=1
    for keys in fasta2r:
        f1.write('>'+str(i)+'\n')
        f1.write(keys+'\n')
        fasta2number[keys]=str(i)
        f2.write(str(i)+'\t'+str(fasta2r[keys])+'\n') 
        i=i+1
    f1.close()
    f2.close()

    f=open(output, 'w')
    genomes=files2genome.values()
    for genome in genomes:
       string=''
       for fasta in fasta2r.keys():
           if genome in fasta2r[fasta]:
               string=string+' '+fasta2number[fasta]
       f.write(genome+'\t'+string+'\n')
    f.close()

def main():
    parser = argparse.ArgumentParser(description="Do something.")
    parser.add_argument('-g', '--genomefile', type=argparse.FileType('r'), default=sys.stdin,required=True)
    parser.add_argument('-o', '--outputindex', required=True)
    parser.add_argument('-m', '--max_mismatch', required=True)
#    parser.add_argument('values', type=float, nargs='*')
    args = parser.parse_args()
    files2genome=genome2spacerfasta(args.genomefile)
    fasta2genome=spacerfastacluster(files2genome)
    fasta2r=comparefasta(fasta2genome,args.max_mismatch)
    buildindex(fasta2r, files2genome, args.outputindex)
         

if __name__== '__main__':
    main()
