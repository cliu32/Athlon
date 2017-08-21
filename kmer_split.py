import pdb
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import numpy as np
import pandas as pd
from read_obj import read_obj
from read_obj import subread_obj
import argparse


def globalPaths():
    global alignFolder
    global refFolder
    global readsFolder
    global SNPFolder
    global trueSNPFolder
    global errorDir
    global consensusFolder
    global primerFolder
    global kallistoPseudoFolder
    global kallistoQuantFolder

    alignFolder='./alignment/'
    refFolder='./reference/'
    readsFolder='./data/'
    SNPFolder='./SNPCall/'
    trueSNPFolder='./trueSNP/'
    errorDir='./errorProfile/'
    consensusFolder='./consensus/'
    primerFolder='./primer/'
    kallistoPseudoFolder='./kallisto_pseudo/'
    kallistoQuantFolder='./kallisto_quant/'

def kmerDistr(seq,k):
    # seq is string or list of char
    kgrams=[seq_obj(''.join(x)) for x in zip(*[seq[i:] for i in xrange(k)])]
    h_rc={}
    [h_rc.update({aa.h:aa.h_rc}) for aa in kgrams]
    return counter_dict(kgrams),h_rc
    
def counter_dict(kgrams):
    cd={}
    for kmer in kgrams:
        try: 
            cd[kmer.h]+=1
        except KeyError:
            try:
                cd[kmer.h_rc]+=1
            except KeyError:
                cd[kmer.h]=1
    return cd

def matchReadPrimer(readLD,readRD,primers,h_rc):
    count=0
    primer_true=''
    current_score=0
    current_score_2nd=0
    for idx,row in primers.iterrows():
        fwD=row['fw_kmer']
        rvD=row['rv_kmer']
        args_list=[(readLD,fwD,row['fw_count'],h_rc),(readLD,rvD,row['rv_count'],h_rc),(readRD,fwD,row['fw_count'],h_rc),(readRD,fwD,row['rv_count'],h_rc)]
        score=max([scoreReadPrimer(*ar) for ar in args_list])
        if score>current_score:
            current_score_2nd=current_score
            current_score=score
            primer_true=idx
    return primer_true,current_score,current_score-current_score_2nd
        
def scoreReadPrimer(readDistr,primerDistr,count_primer,h_rc):
    count=0.
    for key,val in primerDistr.iteritems():
#        inRead=[key==x or h_rc[key]==x for x in readDistr.keys()]
#        if any(inRead):
#            assert sum(inRead)==1
#            count+=min(readDistr[key],val)
        try: 
            count+=min(readDistr[key],val)
        except KeyError:
            try:
                count+=min(readDistr[h_rc[key]],val)
            except KeyError:
                continue
    if count>count_primer: pdb.set_trace()
    return count/count_primer



def loadPrimer_temp(primerFN):
    primerFPath=primerFolder+primerFN
    hlaPrimers=pd.DataFrame.from_csv(primerFPath+'.csv',header=1)
    columns=['fw', 'rv']
    primers=pd.DataFrame(data=hlaPrimers,columns=columns)
    return primers

def splitReadsByPrimer_kmer(readsFN,primerFN,cutoff_list,k=7,endN=100,cutoff_option='score_diff'):
    primers=loadPrimer_temp(primerFN)
    cutoff=np.array(cutoff_list)
# test whether some of the files already exist
    fN_list=[readsFN+'_'+str(key)+'_k'+str(k)+'_endN'+str(endN)+'_cutoff'+str(cutoff[idx]) for key in primers.index for idx in xrange(len(cutoff))]
    doesExist=lambda i: os.path.isfile(readsFolder+fN_list[i]+'.fastq') and os.path.isfile(readsFolder+fN_list[i]+'.p')
    cutoff_notExist=[cutoff_list[i] for i in xrange(len(cutoff_list)) if not doesExist(i)]
    if len(cutoff_notExist)==0:
        print('NOTE_SPLIT: Reads split for '+readsFN+' already exist!!!!!!!!!!!')
        dics=[pickle.load(open(readsFolder+fN+'.p','rU')) for fN in fN_list]
        scores=[]
        scores_diff=[]
        mean_list=[]
        std_list=[]
        for dic in dics:
            scores.append(dic['score'])
            scores_diff.append(dic['score_diff'])
            mean_list.append(dic['mean'])
            std_list.append(dic['std'])
    else:
        print('NOTE_SPLIT: Creating new reads split for '+readsFN+'!!!!!!!')
        h_rc={} #h_rc is {h:rc} dict to preserve mapping between seq and its rc.
        results=[kmerDistr(x,k) for x in primers['fw']]
        primers['fw_kmer']=[x[0] for x in results]
        [h_rc.update(x[1]) for x in results]
        primers['fw_count']=[sum(x.values()) for x in primers['fw_kmer']]
        results=[kmerDistr(x,k) for x in primers['rv']]
        primers['rv_kmer']=[x[0] for x in results]
        [h_rc.update(x[1]) for x in results]
        
        primers['rv_count']=[sum(x.values()) for x in primers['rv_kmer']]
        reads=SeqIO.parse(open(readsFolder+readsFN+'.fastq','rU'),'fastq')
    #    primerReadsDic=dict((x,[]) for x in primers.index)
        primerReadsDics=[dict((x,[]) for x in primers.index) for i in xrange(len(cutoff_list))]
        scores=[]
        scores_diff=[]
        for read in reads:
            read_endL=read[:endN]
            read_endR=read[-endN:]
            readLD,_ =kmerDistr(str(read_endL.seq),k)
            readRD,_ =kmerDistr(str(read_endR.seq),k)
            primer,score,score_diff=matchReadPrimer(readLD,readRD,primers,h_rc)
            if primer=='':
                continue
            else:
                if cutoff_option=='score':
                    for idx in np.where(cutoff-score<0)[0]:
                        primerReadsDics[idx][primer].append(read)
                elif cutoff_option=='score_diff':
                    for idx in np.where(cutoff-score_diff<0)[0]:
                        primerReadsDics[idx][primer].append(read)
                else: 
                    raise ValueError('cutoff_option not valid')
            scores.append(score)
            scores_diff.append(score_diff)
    
        meanDics=[dict((x,0.) for x in primers.index) for i in xrange(len(cutoff_list))]
        stdDics=[dict((x,0.) for x in primers.index) for i in xrange(len(cutoff_list))]
        for idx in xrange(len(primerReadsDics)):
            for key,val in primerReadsDics[idx].iteritems():
                len_list=map(len,val)
                meanDics[idx][key]=np.mean(len_list)
                stdDics[idx][key]=np.std(len_list)
        [SeqIO.write(primerReadsDics[idx][key],open(readsFolder+readsFN+'_'+str(key)+'_k'+str(k)+'_endN'+str(endN)+'_cutoff'+str(cutoff[idx])+'.fastq','wr'),'fastq') for key in primers.index for idx in xrange(len(cutoff))]
        mean_list=[x[ky] for x in meanDics for ky in primers.index]
        std_list=[x[ky] for x in stdDics for ky in primers.index]
        [pickle.dump({'score':score,'score_diff':score_diff,'mean':mean,'std':std},open(readsFolder+fN+'.p','wr')) for score,score_diff,fN,mean,std in zip(scores,scores_diff,fN_list,mean_list,std_list)]
    return scores,scores_diff,fN_list,mean_list,std_list

    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="split reads based on k-mer matching with primers")
    parser.add_argument('readsFN',type=str,help='File name for the reads. Fastq format. Do not include .fastq')
    parser.add_argument('primerFN',type=str,help='File name for the primers. csv format. See example included. Do not include .csv')
    parser.add_argument('-k',nargs='?',default=7,type=int,help='the size of the k-mers used to compare with primer')
    parser.add_argument('--cutoff',nargs='+',default=[0.2],type=float,help='the score cutoff used to delete reads with low pairing with any primer. Between 0 and 1.')
    parser.add_argument('-N','--endN',nargs='?',default=100,type=int,help='the number of bp at the ends of each read to consider when comparing with primer.')
    aa=parser.parse_args()
#    print(args)
    pdb.set_trace()
    splitReadsByPrimer_kmer(aa.readsFN,aa.primerFN,aa.cutoff,k=aa.k,endN=aa.endN)
#    

