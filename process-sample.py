# process-sample.py

# extract reads and make mitochondrial calls


import callmito_single
import os
import sys
import argparse
import time

# SETUP

parser = argparse.ArgumentParser(description='process-sample.py')

parser.add_argument('--ref', type=str,help='genome fasta with dictionary and .fai',required=True)
parser.add_argument('--finaldir', type=str,help='final dir for output',required=True)
parser.add_argument('--name', type=str,help='name of sample to process',required=True)
parser.add_argument('--cram', type=str,help='aligned cram file',required=True)
parser.add_argument('--coords',type=str,help='coordinates to extract, numts + chrM',required=True)
parser.add_argument('--mitoFa',type=str,help='mito fasta with index',required=True)
parser.add_argument('--mitoFaRotated',type=str,help='rotated mito fasta with index',required=True)
parser.add_argument('--chainfile',type=str,help='liftover chain fail to convert rotated to original',required=True)
parser.add_argument('--diagnosticTable',type=str,help='table of diagnostic SNPs',required=True)



args = parser.parse_args()

#####################################################################

myData = {} # dictionary for keeping and passing information

myData['finalDir'] = args.finaldir
myData['ref'] = args.ref
myData['sampleName'] = args.name

myData['cramFileName'] = args.cram
myData['coordsFileName'] = args.coords

myData['mitoFa'] = args.mitoFa
myData['mitoFaRotated'] = args.mitoFaRotated


myData['diagnosticTable'] = args.diagnosticTable



myData['chainFile'] = args.chainfile

# get sequence len
inFile = open(myData['mitoFa'] + '.fai','r')
line = inFile.readline()
line = line.rstrip()
line = line.split()
l = int(line[1])
myData['mitoLen'] = l
inFile.close()

myData['roteTake'] = 4000 # take 4000 first and last from the rotated
myData['minAlleleFreq'] = 0.5 # require >= 50% read support

# max coverage for downsampling
myData['maxCoverage'] = 5000


# check that have interval list file
myData['mitoFaIntervalList'] = myData['mitoFa'].replace('.fa','.interval_list')
myData['mitoFaRotatedIntervalList'] = myData['mitoFaRotated'].replace('.fa','.interval_list')

if os.path.isfile(myData['mitoFaIntervalList']) is False:
    print('ERROR! %s not found, please make interval list' % myData['mitoFaIntervalList'])
    sys.exit()

if os.path.isfile(myData['mitoFaRotatedIntervalList']) is False:
    print('ERROR! %s not found, please make interval list' % myData['mitoFaIntervalList'])
    sys.exit()



if myData['finalDir'][-1] != '/':
   myData['finalDir'] += '/'


if os.path.isdir(myData['finalDir']) is False:
    print('Error! output dir %s not does not exist' % myData['finalDir'])
    sys.exit()
    
# setup the output dir
myData['finalDirSample'] = myData['finalDir']  + myData['sampleName']

if os.path.isdir(myData['finalDirSample']) is False:
    print('making ',myData['finalDirSample'])
    cmd = 'mkdir ' + myData['finalDirSample']
    print(cmd)
    callmito_single.runCMD(cmd)
myData['finalDirSample'] += '/'    

myData['logFileName'] = myData['finalDirSample'] + myData['sampleName'] + '.mito.log'
myData['logFile'] = open(myData['logFileName'],'w')



# add initial infoto log
callmito_single.init_log(myData)
callmito_single.check_prog_paths(myData)


# get reads to extract
callmito_single.extract_reads(myData)



# align to each mito
callmito_single.align_to_mitos(myData)




callmito_single.run_coverage(myData)




callmito_single.down_sample(myData)



# run vcf
callmito_single.call_vars(myData)


# filter vcf
callmito_single.filter_germline(myData)



# make fasta and mask
callmito_single.make_fasta_germline(myData)

callmito_single.assign_haplogroup(myData)

myData['logFile'].close()