# callmito.py
# functions for going from aligned CRAM to mitochondrial calls

import sys
import subprocess
import os
import argparse
import time
import socket
import shutil
import gzip
import numpy as np

###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed')
        print(cmd)
        sys.exit(1)
###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD_output(cmd):
    val = subprocess.Popen(cmd, universal_newlines=True, shell=True, stdout = subprocess.PIPE)
    resLines = []
    for i in val.stdout:
       i = i.rstrip()
       resLines.append(i)
    return resLines
#############################################################################        
# setup paths to default programs to use and checks for required programs
def check_prog_paths(myData):        
    myData['logFile'].write('\nChecking for required programs...\n')
    
    for p in ['bwa','gatk','samtools','liftOver','bgzip','tabix','bcftools']:
        if shutil.which(p) is None:
            s = p + ' not found in path! please fix (module load?)'
            print(s, flush=True)
            myData['logFile'].write(s + '\n')
            myData['logFile'].close()        
            sys.exit()
        else:
            myData['logFile'].write('%s\t%s\n' % (p,shutil.which(p)))
            

    check_gatk_version(myData)

    myData['logFile'].flush()              


############################################################################# 
def check_gatk_version(myData):
    myData['tmpVersionName'] = myData['finalDirSample'] + 'tmp.gatk.version'

    gatk_v = runCMD_output('gatk --version > %s' % myData['tmpVersionName'])
    inFile = open(myData['tmpVersionName'],'r')    
    line = inFile.readline()
    inFile.close()
    line = line.rstrip()
    line = line.split()
    v = line[-1]
        
    s = 'GATK version is: %s' % v
    print(s,flush=True)
    myData['logFile'].write(s + '\n')

    if v != 'v4.2.5.0':
        s = 'ERROR! GATK v4.2.5.0 is required!.  Please fix'
        print(s, flush=True)
        myData['logFile'].write(s + '\n')
        myData['logFile'].close()        
        sys.exit()
############################################################################# 
def init_log(myData):
    k = list(myData.keys())
    k.sort()
    myData['startTime'] = time.localtime()
    myData['tStart'] = time.time()
    t = time.strftime("%a, %d %b %Y %H:%M:%S", myData['startTime'])        
    myData['logFile'].write(t + '\n')
    
    hn = socket.gethostname()
    myData['logFile'].write('Host name: %s\n' % hn)
    print('Host name: %s\n' % hn,flush=True)
    
    myData['logFile'].write('\nInput options:\n')
    for i in k:
        if i in ['logFile']:
            continue        
        myData['logFile'].write('%s\t%s\n' % (i,myData[i]))                
    myData['logFile'].flush()  
############################################################################# 
def parse_sam_line(myLine):
    res = {}
    res['seqName'] = myLine[0]
    res['flag'] = int(myLine[1])
    res['chrom'] = myLine[2]
    res['chromPos'] = int(myLine[3])
    res['mapQ'] = int(myLine[4])
    res['cigar'] = myLine[5]
    res['seq'] = myLine[9]
    
    #
    res['seqLen'] = len(myLine[9])
    
    
    res['cigarExpand'] = expand_cigar(res['cigar'])
    res['qual'] = myLine[10]
    res['mateChrom'] = myLine[6]
    res['matePos'] = myLine[7]
    
    res['fragLen'] = int(myLine[8])
    
    res['cigarCounts']={}
    res['cigarCounts']['M'] = 0
    res['cigarCounts']['D'] = 0
    res['cigarCounts']['I'] = 0
    res['cigarCounts']['S'] = 0
    res['cigarCounts']['H'] = 0

    
    if res['flag'] & 0x10 != 0:
        res['reverseStrand'] = True
    else:
        res['reverseStrand'] = False

    if res['flag'] & 0x4 != 0:
        res['unMapped'] = True
    else:
        res['unMapped'] = False

    if res['flag'] & 0x400 != 0:
        res['isDuplicate'] = True
    else:
        res['isDuplicate'] = False

    if res['flag'] & 0x100 != 0:
        res['notPrimaryAlignment'] = True
    else:
        res['notPrimaryAlignment'] = False

    if res['flag'] & 0x800 != 0:
        res['isSupplementaryAlignment'] = True
    else:
        res['isSupplementaryAlignment'] = False


    if res['flag'] & 0x1 != 0:
        res['isPaired'] = True
    else:
        res['isPaired'] = False


    if res['flag'] & 0x8 != 0:
        res['mateUnmapped'] = True
    else:
        res['mateUnmapped'] = False

    if res['flag'] & 0x40 != 0:
        res['isFirst'] = True
    else:
        res['isFirst'] = False

        
    for i in res['cigarExpand']:
        res['cigarCounts'][i[1]] += i[0]
        
    # check for proper seqlen to update, 2015-05-05
    if myLine[9] == '*':  #not actually sequence present in SAM line
        res['seqLen'] = res['cigarCounts']['M'] + res['cigarCounts']['I']  + res['cigarCounts']['S'] + res['cigarCounts']['H']
     
     
    res['otherTags'] = myLine[11:]
             
    return res
#####################################################################
#returns lists of [int,flag]
def expand_cigar(cigar):
    res = []
    if cigar == '*':
        return res
    digits = ['0','1','2','3','4','5','6','7','8','9']
    accumulate = ''
    i = 0
    while True:
        if i == len(cigar):
            break
        if cigar[i] in digits:
            accumulate += cigar[i]
            i += 1
        else:
            d = int(accumulate)
            res.append([d,cigar[i]])
            i += 1
            accumulate = ''
    return res
#####################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c   
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''    
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
##############################################################################
def get_seq_from_sam(samRec):
    name = samRec['seqName']
    if samRec['isFirst'] is True:
        readNum = 1
    else:
        readNum = 2
    
    # get if the seq...
    if samRec['reverseStrand'] is True:
        seq = samRec['seq']
        seq = revcomp(seq)
        qual =  samRec['qual']
        qual = qual[::-1] #reverse
    else:
        seq = samRec['seq']
        qual =  samRec['qual']
    return([name,readNum,seq,qual])         
############################################################################# 
# define critera to extract reads for remapping
def to_extract(samRec):
    if samRec['unMapped'] is True:
        return False

    if samRec['isDuplicate'] is True:
        return False

    if samRec['notPrimaryAlignment'] is True:
        return False
     
     # this is callmito_single -- we will only take single end reads!

    if samRec['isPaired'] is True:
        return False

    return True
###############################################################################        


###############################################################################        
def extract_reads(myData):
    searchDelta = 100 # space to search across
    
    myData['logFile'].write('\nstarting extraction of fastq\n')
    
    # due this one through a temp file name, as it may be very large and can break the pipe
    
    myData['tmpSamFileName'] = myData['finalDirSample'] + 'tmp.extract.sam'
    
    cmd = 'samtools view -T %s -M -L %s -o %s -O SAM %s  ' % (myData['ref'], myData['coordsFileName'],myData['tmpSamFileName'],myData['cramFileName'])
    print(cmd)
    runCMD(cmd)    
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].flush()
    myData['logFile'].write('DONE' + '\n')
    myData['logFile'].flush()
    print('DONE initial extraction',flush=True)
    
    # these dictionaries can be large, requires some more memory when there is a lot of mito reads
    myData['readData'] = {}  # dictionary to store all of read 1 and read 2 info
    myData['readsToExtract'] = {} # dictionary of reads to extract
    
    myData['logFile'].write('starting to read through extracted file' + '\n')
    myData['logFile'].flush()
  
    tmpIn = open(myData['tmpSamFileName'],'r')
    for samLine in tmpIn:
        samLine = samLine.rstrip()
        samLine = samLine.split()
        samRec = parse_sam_line(samLine)
        
        
        readNum = 1 # has to be -- only a single end read...
 
        if samRec['seqName'] not in myData['readData']:
            myData['readData'][samRec['seqName']] = ['Empty','Empty']
        # just save the samLine, which is the split list of the samfile line, this saves space    
        myData['readData'][samRec['seqName']][readNum-1] = samLine        
        if to_extract(samRec) is True:
            myData['readsToExtract'][samRec['seqName']] = 1   
    tmpIn.close()
    
    s = 'Have total of %i reads pass extraction criteria' % len(myData['readsToExtract'])    
    print(s,flush=True)
    myData['logFile'].write(s + '\n')

    # get how many need extraction
        
        
    myData['fastqOutName'] = myData['finalDirSample'] + 'read.fq.gz'
     
    out1 = gzip.open(myData['fastqOutName'],'wt')
    for rn in myData['readsToExtract']:
        samRec = parse_sam_line(myData['readData'][rn][0])
        seqInfo = get_seq_from_sam(samRec)
        s = seqInfo[2]
        q = seqInfo[3]
        
        out1.write('@%s\n%s\n+\n%s\n' % (rn,s,q))

    out1.close()
     
    s = 'reads written to output fastq files!'
    print(s,flush=True)
    myData['logFile'].write(s + '\n')          
    myData['logFile'].flush() 
    
    # free up memory
    myData['readsToExtract'].clear()
    myData['readData'].clear()
    
    
    cmd = 'rm ' + myData['tmpSamFileName']
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')          
    myData['logFile'].flush() 
    runCMD(cmd)


    
###############################################################################        
def align_to_mitos(myData):
    s = 'align to the two mitos'
    print(s,flush=True)
    myData['logFile'].write(s + '\n')          
    
    myData['mitoBam']= myData['finalDirSample'] + 'mito.bam'
    myData['mitoRotatedBam']= myData['finalDirSample'] + 'mitoRotated.bam'    

    saiTMP = myData['finalDirSample'] + 'TMP.sai'

    
    # align to mito
    rg = '@RG\\tID:norm\\tSM:%s\\tPL:Illumina' % (myData['sampleName'])
    rg = '\'' + rg + '\''
    print('rg is',rg)
    
    # do bwa aln/sampe thing...
    cmd = f"bwa aln -t 1 -l 1024 -n 0.01 -o 2 {myData['mitoFa']}  {myData['fastqOutName']} > {saiTMP}"
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].flush()                      
    runCMD(cmd)
    
    cmd = f"bwa samse -r {rg} {myData['mitoFa']} {saiTMP} {myData['fastqOutName']} | samtools view -F 4 -h -b - | samtools sort - > {myData['mitoBam']} "
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].flush()                      
    runCMD(cmd)
    
    
    
    
    # align to rotated mito
    rg = '@RG\\tID:rotate\\tSM:%s\\tPL:Illumina' % (myData['sampleName'])
    rg = '\'' + rg + '\''
    print('rg is',rg)

    cmd = f"bwa aln -t 1 -l 1024 -n 0.01 -o 2 {myData['mitoFaRotated']}  {myData['fastqOutName']} > {saiTMP}"
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].flush()                      
    runCMD(cmd)
    
    cmd = f"bwa samse -r {rg} {myData['mitoFaRotated']} {saiTMP} {myData['fastqOutName']} | samtools view -F 4 -h -b - | samtools sort - > {myData['mitoRotatedBam']} "
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].flush()                      
    runCMD(cmd)
    
    
    cmd = 'rm ' + saiTMP
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')
    myData['logFile'].flush()                      
    runCMD(cmd)
    
    
    


    # sort and markdups
    myData['mitoBamSort'] = myData['finalDirSample'] + 'mito.sort.bam'
    myData['mitoRotatedBamSort'] = myData['finalDirSample'] + 'mitoRotated.sort.bam'    
    
    cmd = 'gatk SortSam -SO coordinate -I %s -O %s ' % (myData['mitoBam'],myData['mitoBamSort'])
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)

    cmd = 'gatk SortSam -SO coordinate -I %s -O %s ' % (myData['mitoRotatedBam'],myData['mitoRotatedBamSort'])
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    myData['logFile'].flush()    

    runCMD(cmd)
    
    # mark duplicates
    myData['mitoBamSortMD'] = myData['finalDirSample'] + 'mito.sort.markdup.bam'
    myData['mitoRotatedBamSortMD'] = myData['finalDirSample'] + 'mitoRotated.sort.markdup.bam'    
    
    myData['mitoBamDupMet'] = myData['finalDirSample'] + 'mito.dup_metrics.txt'
    myData['mitoRotatedDupMet'] = myData['finalDirSample'] + 'mitoRotated.dup_metrics.txt'    
    
    cmd = 'gatk MarkDuplicates -I %s -O %s -M %s ' % (myData['mitoBamSort'],myData['mitoBamSortMD'],myData['mitoBamDupMet'])
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)

    cmd = 'gatk MarkDuplicates -I %s -O %s -M %s ' % (myData['mitoRotatedBamSort'],myData['mitoRotatedBamSortMD'],myData['mitoRotatedDupMet'])
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    
    # index    
    cmd = 'samtools index %s' % myData['mitoBamSortMD']
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)

    cmd = 'samtools index %s' % myData['mitoRotatedBamSortMD']
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    myData['logFile'].flush()    
    
###############################################################################        
def run_coverage(myData):
# get covergage
    
    myData['mitoHSmets'] =  myData['finalDirSample'] + 'mito.hsmets.txt'
    myData['mitoPerBp'] = myData['finalDirSample'] + 'mito.per-base.txt'

    myData['mitoRotatedHSmets'] =  myData['finalDirSample'] + 'mitoRotated.hsmets.txt'
    myData['mitoRotatedPerBp'] = myData['finalDirSample'] + 'mitoRotated.per-base.txt'


    
    cmd = 'gatk CollectHsMetrics -I %s -O %s -R %s -PER_BASE_COVERAGE %s --COVERAGE_CAP 50000 --SAMPLE_SIZE 1 -BI %s -TI %s ' % ( myData['mitoBamSortMD'],
              myData['mitoHSmets'],myData['mitoFa'],myData['mitoPerBp'],myData['mitoFaIntervalList'],myData['mitoFaIntervalList']  )
              
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)

    cmd = 'gatk CollectHsMetrics -I %s -O %s -R %s -PER_BASE_COVERAGE %s --COVERAGE_CAP 50000 --SAMPLE_SIZE 1 -BI %s -TI %s ' % ( myData['mitoRotatedBamSortMD'],
              myData['mitoRotatedHSmets'],myData['mitoFaRotated'],myData['mitoRotatedPerBp'],myData['mitoFaRotatedIntervalList'],myData['mitoFaRotatedIntervalList']  )
              
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    
    
    # now, do the conversion
    
    
    mitoDepth = {}
    inFile = open(myData['mitoPerBp'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        if line[0] == 'chrom':
            continue
        line[1] = int(line[1])
        line[3] = int(line[3])
        mitoDepth[line[1]] = line[3]
    inFile.close()

    mitoRotDepth = []
    inFile = open(myData['mitoRotatedPerBp'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        if line[0] == 'chrom':
            continue
        line[1] = int(line[1])
        line[3] = int(line[3])
        mitoRotDepth.append([line[0],line[1],line[3]])
    inFile.close()
    
    # make tmp rotated
    myData['mitoRotatedPerBpBED'] = myData['mitoRotatedPerBp'] + '.bed'
    myData['mitoRotatedPerBpBEDlift'] =  myData['mitoRotatedPerBpBED'] + '.lift'
    myData['mitoRotatedPerBpBEDliftFail'] =  myData['mitoRotatedPerBpBED'] + '.liftFail'    
    outFile = open(myData['mitoRotatedPerBpBED'],'w')
    for r in mitoRotDepth:
        k = str(r[1]) + ':' + str(r[2])
        nl = '%s\t%i\t%i\t%s\n' % (r[0],r[1]-1,r[1],k)
        outFile.write(nl)    
    outFile.close()
    
    # run liftover
    
    cmd = 'liftOver %s %s %s %s' % (myData['mitoRotatedPerBpBED'],myData['chainFile'],myData['mitoRotatedPerBpBEDlift'],myData['mitoRotatedPerBpBEDliftFail'])
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
        
    
    # read in the new one...
    mitoRotDepth = {}
    inFile = open(myData['mitoRotatedPerBpBEDlift'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        if line[0] == 'chrom':
            continue
        p = int(line[2])
        k = line[3]
        d = k.split(':')[1]
        d = int(d)
        mitoRotDepth[p] = d
    inFile.close()
    # now have to do the merge.....
    
    # now can output the depth
    myData['mitoMergePerBp'] =  myData['finalDirSample'] + 'mitoMerge.per-bp.txt'
    outFile = open(myData['mitoMergePerBp'],'w')
    
    tr = 0
    allDepth = []
    for i in range(1,myData['mitoLen']+1):
        if i <= myData['roteTake'] or i >= (myData['mitoLen']-myData['roteTake'] +1 ):
            tr += 1
            d = mitoRotDepth[i]
        else:
            d= mitoDepth[i]
        allDepth.append(d)
        outFile.write('%i\t%i\n' % (i,d))
    outFile.close()
    print('took %i from rotates' % tr)
    
    myData['meanDepth'] = np.mean(allDepth)
    myData['minDepth'] = min(allDepth)
    myData['maxDepth'] = max(allDepth)
    myData['medDepth'] = np.median(allDepth)
    
    myData['mitoMergePerBpStats'] = myData['mitoMergePerBp'].replace('.txt','.stats')
    outFile = open(myData['mitoMergePerBpStats'],'w')
    for i in ['meanDepth','medDepth','minDepth','maxDepth']:
        outFile.write('%s\t%i\n' % (i,myData[i]))
        print('%s\t%i' % (i,myData[i]))
    outFile.close()
    
    myData['logFile'].flush()    
    
    
###################################################################################################
def down_sample(myData):
# make new downsampled bams if coverage is too high
    seed = 1983 # so that both bams are downsamples the same
    s = 'mean depth is %f ' % myData['meanDepth']
    if myData['meanDepth'] <= myData['maxCoverage']:
        s += ' less than max of %f, OK!' % myData['maxCoverage']
        print(s,flush=True)
        myData['logFile'].write(s + '\n')  
        return; # no need to downsample
    # run downsample
    f = myData['maxCoverage'] / myData['meanDepth'] 
    s += ' more than max of %f, run downsample %f!' % (myData['maxCoverage'],f)
    print(s,flush=True)
    myData['logFile'].write(s + '\n')  

    myData['mitoBamOrig'] = myData['finalDirSample'] + 'ORIGINIAL.mito.sort.markdup.bam'
    myData['mitoRotatedBamOrig'] = myData['finalDirSample'] + 'ORIGINIAL.mitoRotated.sort.markdup.bam'    

    cmd = 'mv %s %s' % (myData['mitoBamSortMD'],myData['mitoBamOrig'] )
    print(cmd)
    myData['logFile'].write(cmd + '\n')  
    runCMD(cmd)   

    cmd = 'samtools index %s' % (myData['mitoBamOrig'] )
    print(cmd)
    myData['logFile'].write(cmd + '\n')  
    runCMD(cmd)   

    cmd = 'mv %s %s' % (myData['mitoRotatedBamSortMD'],myData['mitoRotatedBamOrig'] )
    print(cmd)
    myData['logFile'].write(cmd + '\n')  
    runCMD(cmd)   

    cmd = 'samtools index %s' % (myData['mitoRotatedBamOrig'] )
    print(cmd)
    myData['logFile'].write(cmd + '\n')  
    runCMD(cmd)   

    myData['logFile'].flush()

    # run subsample    
    cmd = 'gatk DownsampleSam -I %s -O %s -P %f -R %i '  % (myData['mitoBamOrig'],myData['mitoBamSortMD'],f,seed)
    print(cmd)
    myData['logFile'].write(cmd + '\n')  
    runCMD(cmd)   
    cmd = 'samtools index %s' % myData['mitoBamSortMD']
    print(cmd)
    myData['logFile'].write(cmd + '\n')  
    runCMD(cmd)   

    cmd = 'gatk DownsampleSam -I %s -O %s -P %f -R %i '  % (myData['mitoRotatedBamOrig'],myData['mitoRotatedBamSortMD'],f,seed)
    print(cmd)
    myData['logFile'].write(cmd + '\n')  
    runCMD(cmd)   
    cmd = 'samtools index %s' % myData['mitoRotatedBamSortMD']
    myData['logFile'].write(cmd + '\n')      
    runCMD(cmd)   
    print(cmd)
    
    myData['logFile'].flush()
    
###################################################################################################
def call_vars(myData):
# call the mitochondrial variants
    myData['mitoVCF'] = myData['finalDirSample'] + 'mito.vcf.gz'
    myData['mitoRotatedVCF'] = myData['finalDirSample'] + 'mitoRotated.vcf.gz'
    myData['mitoRotatedVCFLift'] = myData['finalDirSample'] + 'mitoRotated.LIFT.vcf.gz'
    myData['mitoRotatedVCFLiftFail'] = myData['finalDirSample'] + 'mitoRotated.LIFT-FAIL.vcf.gz'       
    
    myData['mitoMergeVCF'] =  myData['finalDirSample'] + 'mitoMerged.vcf'
    
    myData['mitoVCFFilter'] = myData['mitoVCF'] + '.filter.gz'
    myData['mitoRotatedVCFFilter'] = myData['mitoRotatedVCF'] + '.filter.gz'    
    

    cmd = 'gatk Mutect2 --max-reads-per-alignment-start 75 --max-mnp-distance 0 -R %s --mitochondria-mode -I %s --annotation StrandBiasBySample -O %s' % (myData['mitoFa'],myData['mitoBamSortMD'],myData['mitoVCF']) 
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)

    cmd = 'gatk Mutect2 --max-reads-per-alignment-start 75 --max-mnp-distance 0 -R %s --mitochondria-mode -I %s --annotation StrandBiasBySample -O %s' % (myData['mitoFaRotated'],myData['mitoRotatedBamSortMD'],myData['mitoRotatedVCF']) 
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)


    # filter..
    cmd = 'gatk FilterMutectCalls --mitochondria-mode -R %s -V %s -O %s ' % (myData['mitoFa'],myData['mitoVCF'],myData['mitoVCFFilter'] )
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)

    cmd = 'gatk FilterMutectCalls --mitochondria-mode -R %s -V %s -O %s ' % (myData['mitoFaRotated'],myData['mitoRotatedVCF'],myData['mitoRotatedVCFFilter'] )
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)


    # run liftover vcf    
    cmd = 'gatk LiftoverVcf -I %s -O %s -CHAIN %s -REJECT %s -R %s ' % (myData['mitoRotatedVCFFilter'],myData['mitoRotatedVCFLift'],myData['chainFile'],myData['mitoRotatedVCFLiftFail'],myData['mitoFa'] )    
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    
    # do the read in and merge
    # get rid of PASS annotation...
    
    
    # read in 
    liftedVCF = []
    inFile = gzip.open(myData['mitoRotatedVCFLift'],'rt')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        line[1] = int(line[1])
        #line[6] = '.' # keep filter results
        liftedVCF.append(line)
    inFile.close()
    print('read in %i from %s' % (len(liftedVCF),myData['mitoRotatedVCFLift']))

    mitoVCF = []
    inFile = gzip.open(myData['mitoVCFFilter'],'rt')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        line[1] = int(line[1])
        #line[6] = '.' # keep filter results
        mitoVCF.append(line)
    inFile.close()
    print('read in %i from %s' % (len(mitoVCF),myData['mitoVCF']))

    outFile = open(myData['mitoMergeVCF'],'w')
    
    # header
    inFile = gzip.open(myData['mitoVCFFilter'],'rt')
    for line in inFile:
        if line[0] == '#':
            outFile.write(line)
    inFile.close()
    
    # part 1
    for row in liftedVCF:
        if row[1] <= myData['roteTake']:
            nl = row
            nl = [str(i) for i in nl]
            nl = '\t'.join(nl) + '\n'
            outFile.write(nl)
    # middle part
    for row in mitoVCF:
        if row[1] > myData['roteTake'] and row[1] < (myData['mitoLen']-myData['roteTake'] +1 ):
            nl = row
            nl = [str(i) for i in nl]
            nl = '\t'.join(nl) + '\n'
            outFile.write(nl)
    
    # end part
    for row in liftedVCF:
        if row[1] >= (myData['mitoLen']-myData['roteTake'] +1 ):
            nl = row
            nl = [str(i) for i in nl]
            nl = '\t'.join(nl) + '\n'
            outFile.write(nl)
    outFile.close()
    
    # compress and tabix
    cmd = 'bgzip %s' % myData['mitoMergeVCF']
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    
    myData['mitoMergeVCF']+= '.gz'
    cmd = 'tabix -p vcf %s' % myData['mitoMergeVCF']
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    myData['logFile'].flush()    
    
    
###################################################################################################
def filter_germline(myData):
# filter out for germline calls
# parses output -- have already run  gatk FilterMutectCalls
    myData['mitoMergeVCFFilter'] =  myData['finalDirSample'] + myData['sampleName'] + '.mitoMerged.germline.filter.vcf'
    myData['mitoMergeNonRefFraction'] = myData['finalDirSample'] + myData['sampleName'] + '.nonRefFraction.txt'

    outStats = open(myData['mitoMergeNonRefFraction'],'w')

    inFile = gzip.open(myData['mitoMergeVCF'],'rt')
    outFile = open(myData['mitoMergeVCFFilter'],'w')
    for line in inFile:
        if line[0] == '#':
            outFile.write(line)
            continue
        line = line.rstrip()
        ol = line
        line = line.split()
        
        infoDict = parse_vcf_info(line[7])
        genoDict = parse_genotype(line[8],line[9])

        
                
        # check num alt alleles
        alts = line[4]
        alts = alts.split(',')

        # get dp
        dp = genoDict['AD']
        altIndexmaxAltAlleleFeq = 0 # is index of alleles, not of alt allels
        maxAltAlleleFeq = 0.0
        tot = 0
        for i in range(len(dp)):
            tot += int(dp[i])
        for i in range(1,len(dp)):
            d = int(dp[i])
            if tot == 0:
                f = 0.0
            else:
                f = d/tot
            if f > maxAltAlleleFeq:
                maxAltAlleleFeq = f
                altIndexmaxAltAlleleFeq = i
            

        # check filters, see if max alt allele passess filters
        # only look for strand_bias as an artefact
        AS_Filters = infoDict['AS_FilterStatus']        
        if 'strand_bias' in AS_Filters[altIndexmaxAltAlleleFeq-1]:
            s = 'fails strand bias,allele index is %i' % altIndexmaxAltAlleleFeq
            s += '\n' + ol
            print(s,flush=True)
            myData['logFile'].write(s + '\n')                                                
            continue

        outStats.write('%f\n' % maxAltAlleleFeq)        # print out the max alt alle freq
        if maxAltAlleleFeq < myData['minAlleleFreq']:  # max alt allele freq is no good, so skip it...
            continue
        


        
        line[6] = 'PASS'        
        # edit the gen
        gen = line[9]
        gen = gen.split(':')
        gen[0] = str(altIndexmaxAltAlleleFeq) + '/' + str(altIndexmaxAltAlleleFeq) # make it homozygous
        line[9] = ':'.join(gen)
        
        nl = '\t'.join(line) + '\n'
        outFile.write(nl)
            
    inFile.close()
    outFile.close()  
    # convert to gz
    outStats.close()

    cmd = 'bgzip %s' % myData['mitoMergeVCFFilter']
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    
    myData['mitoMergeVCFFilter']+= '.gz'
    cmd = 'tabix -p vcf %s' % myData['mitoMergeVCFFilter']
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)
    
    myData['logFile'].flush()    
    
    
###################################################################################################
def make_fasta_germline(myData):
    minMitoDepth = 3
    myData['mitoMergeMasked'] =  myData['finalDirSample'] + 'mask-regions.bed'
    myData['mitoMergeFasta'] = myData['finalDirSample'] + myData['sampleName'] + '.fa'
    
    tmpFa = myData['mitoMergeMasked'] + '.tmp.fa'
    
    
    # first setup regions to mask, includes hard coded regions
    # and any region with depth < 100
    outFile = open(myData['mitoMergeMasked'],'w')
    outFile.write('NC_002008.4\t15989\t16600\n')
    outFile.write('NC_002008.4\t15511\t15535\n')    

    
    alreadyMasked = {}
    for i in range(15990,15990+1):
        alreadyMasked[i] = 1
    for i in range(15512,15535+1):
        alreadyMasked[i] = 1
        
    failDepthMask = 0
    inFile = open(myData['mitoMergePerBp'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        pos = int(line[0])
        d = int(line[1])
        if d < minMitoDepth and pos not in alreadyMasked:
            failDepthMask +=1
            outFile.write('NC_002008.4\t%i\t%i\n' % (pos-1,pos))
    inFile.close()
    
    s = 'found %i positions that failed depth check %i' % (failDepthMask,minMitoDepth)
    print(s,flush=True)
    myData['logFile'].write(s + '\n')              

    s = 'checking for overlapping vcf intervals'
    print(s,flush=True)
    myData['logFile'].write(s + '\n')          
    
    intsToMask = []
    
    prevStart = 0
    prevEnd = 0
    inFile = gzip.open(myData['mitoMergeVCFFilter'],'rt')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        pos = int(line[1])
        ref = line[3]
        refLen = len(ref)
        posEnd = pos+refLen - 1
        
        if pos <= prevEnd and pos >= prevStart:
            intsToMask.append([prevStart,prevEnd])
            intsToMask.append([pos,posEnd])
        prevStart = pos
        prevEnd = posEnd
    inFile.close()
    
    s = 'found %i intervals to mask that overlap' % len(intsToMask)
    print(s,flush=True)
    myData['logFile'].write(s + '\n')   
    myData['logFile'].flush()    
           
    
    for r in intsToMask:
        posStart = r[0]-1
        posEnd = r[1]
        outFile.write('NC_002008.4\t%i\t%i\n' % (posStart,posEnd))
    outFile.close()
    
    # make fasta

    cmd = 'bcftools consensus --sample %s -f %s -m %s %s > %s ' % (myData['sampleName'], myData['mitoFa'],myData['mitoMergeMasked'],myData['mitoMergeVCFFilter'], tmpFa )
    print(cmd,flush=True)
    myData['logFile'].write(cmd + '\n')              
    runCMD(cmd)

    inFile = open(tmpFa,'r')
    outFile = open(myData['mitoMergeFasta'],'w')
    for line in inFile:
        if line[0] == '>':
            outFile.write('>%s\n' % myData['sampleName'])
        else:
            outFile.write(line)
    inFile.close()
    outFile.close()      
    myData['logFile'].flush()    
      
#############################################################################    
def assign_haplogroup(myData):
    # read in diagnostic table
    myData['mitoMergeHaploGroup'] = myData['finalDirSample'] + myData['sampleName'] + '.haplogroup.txt'
    
    inFile = open(myData['diagnosticTable'],'r')
    haplos = []
    for line in inFile:
        line = line.rstrip()
        line = line.split('\t')
        
        if line[0] == '#Haplogroup':
            header = line
        elif line[0] == '#haplogroup':
            break
        else:
            haplos.append(line)        
    definedHaps = []
    haploGroupOrder = line[1:]
    for line in inFile:
        line = line.rstrip()
        line = line.split('\t')
        hapName = line[0]
        hapSeq = line[1:]
        definedHaps.append([hapName,hapSeq])

    print('read in the hap seqs!',len(definedHaps))
    inFile.close()
       
    # read in all the SNPs
    inFile = gzip.open(myData['mitoMergeVCFFilter'],'rt')
    snps = {}
    for line in inFile:
        if line[0] == '#':
             continue
        line = line.rstrip()
        line = line.split()
        
        pos = line[1]
        ref = line[3]
        alt = line[4]
        i = ref + '-' + pos + '-' + alt
        snps[i] = 1
    inFile.close()
    
    print('read in snps to compare',len(snps))
    outFile = open(myData['mitoMergeHaploGroup'],'w')
    outFile.write('#Haplogroup\tNumber of SNPs\tSNPs\tNumber Present\n')
    
    pToSNPrecord = {}
    outputRows = []
    for hap in haplos:
        numFound = 0
        hapSnp = hap[2].split(';')
        for h in hapSnp:
            p = h.split('-')[1]
            pToSNPrecord[p] = h
            if h in snps:
                numFound += 1
        # print them all out I suppose..
        if numFound > 0:
            nl = hap
            nl.append(str(numFound))
            nl = '\t'.join(nl)
            outFile.write(nl + '\n')
            print(nl)
    
    # now need to make the hap seq to compare with....
    
    sampleSeq = []
    for p in haploGroupOrder:
        print(p,pToSNPrecord[p])
        if pToSNPrecord[p] in snps:
            allele = pToSNPrecord[p].split('-')[2]
        else:
            allele = pToSNPrecord[p].split('-')[0]
        sampleSeq.append(allele)
    print(sampleSeq)
    
    nl = ''.join(sampleSeq)
    outFile.write('\n#Haplotype Assignment\nsample haplotype:\t%s\n' % nl)
    
    # now need to go through and figure out the ones with the closest match
    dists = []
    for defHaps in definedHaps:
        hapName =defHaps[0]
        hapSeq = defHaps[1]
        
        numDiff = 0
        for i in range(len(hapSeq)):
            if hapSeq[i] != sampleSeq[i]:
                numDiff += 1
        dists.append([numDiff,hapName])
            
    # now need to sort
    dists.sort(key=lambda i: i[0])
    print(dists)
    
    if dists[0][0] < dists[1][0]:
        outFile.write('Best match:\t%s\tNum differences:\t%i\n' % (dists[0][1],dists[0][0]))
    else:
        for d in dists:
            if d[0] > dists[0][0]:
                continue
            outFile.write('Tied match:\t%s\tNum differences:\t%i\n' % (d[1],d[0]))    
    outFile.close()
#############################################################################    
# Makes a dictionary of the info field in a vcf file
# returns the dictionary
#example: DP=5;AF1=1;CI95=0.5,1;DP4=0,0,4,1;MQ=51
def parse_vcf_info(infoField):
    info = {}
    infoList = infoField.split(';')
    for field in infoList:
       if field.count('=') == 1:
           (name,vals) = field.split('=')[0:2]
           if name == 'AS_FilterStatus':
               vals = vals.split('|')
           else:
               vals = vals.split(',')
           if vals == 'true' or vals == 'True':
               vals = True
           if vals == 'false' or vals == 'False':
               vals = False
           info[name] = vals
       else:
           info[field] = 'PRESENT'
    return info
#############################################################################                 
def parse_genotype(formatField,genoTypeField):
# parse the genotype into a dictionary, based on format field
    genotypeInfo = {}
    formats = formatField.split(':')
    genFields = genoTypeField.split(':')
    for i in range(len(formats)):
        f= formats[i]
        g = genFields[i]
        g = g.split(',')
        genotypeInfo[f] = g
    return(genotypeInfo)
#############################################################################                 





