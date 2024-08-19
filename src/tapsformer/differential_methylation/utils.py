import pysam
from tapsformer.constants import * 
import numpy as np

valid_base_for_read ={"A", "C", "G", "T","N","D"}

def load_next_pileup_region2(bam, c):
    '''
    load a pileup of the whole chromosome for easy iteration    
    '''
    return bam.pileup(
        c,
        start=0,
        truncate=True,
        stepper="all",
        ignore_overlaps=False,
        compute_baq=False,
        flag_require=3,
    )


def align_read_to_reference(chr, read, genome_handle):
    '''
    returns a read and refernece with N inserted to account for indels and soft skips
    also returns the preceding and following reference base. 
    '''
    reconstructed_read = ""
    reconstructed_ref = ""
    ref_pos = read.reference_start
    cigar_tuples=read.cigartuples
    sequence = read.query_sequence
    ref_sequence=read.get_reference_sequence()
    ref_start = ref_pos
    for i, (cigar_code, length) in enumerate(cigar_tuples):
        if cigar_code == 0:  # match 
            # cut the read and reference sequences as they are
            reconstructed_read += sequence[:length]
            sequence = sequence[length:]
            reconstructed_ref += ref_sequence[:length]
            ref_sequence = ref_sequence[length:]
            ref_pos += length
        elif cigar_code == 1:  # insertion
            # cut the read as it is 
            reconstructed_read += sequence[:length]
            sequence = sequence[length:]
            # insert N*length into the reference
            reconstructed_ref += "N" * length
        elif cigar_code == 2:  # deletion
            # insert N*length to the read 
            reconstructed_read += "N" * length
            # take the reference as is
            reconstructed_ref += ref_sequence[:length]
            ref_sequence = ref_sequence[length:]
            ref_pos += length
        elif cigar_code == 4:  # soft clipping
            sequence = sequence[length:]
            reconstructed_read = reconstructed_read + "N" * length
            if i==0:
                # if in the beginning fetch the reference for the soft skip
                reconstructed_ref += genome_handle.fetch(chr, ref_pos-length, ref_pos)
                ref_start -= length
            else:
                reconstructed_ref += genome_handle.fetch(chr, ref_pos, ref_pos+length)
                ref_pos += length
    reconstructed_read=reconstructed_read.upper()
    reconstructed_ref=reconstructed_ref.upper()
    if reconstructed_read is None:
        return None, None, None,None
    for c in set(reconstructed_ref):
        if not c in valid_base_for_read:
            print("ERROR: invalid base", c)
            return None, None, None,None
    ref_prev_base = genome_handle.fetch(chr, ref_start-1, ref_start)
    ref_following_base = genome_handle.fetch(chr, ref_pos, ref_pos+1)
    return reconstructed_read, reconstructed_ref, ref_prev_base, ref_following_base


def count_cpgs(ref,ref_prev_base, ref_following_base):
    '''
    Count the number of CpGs in the reference.
    '''
    cpg_count = 0
    for iii in range(len(ref)-1):
        if ref[iii:iii+2]=="CG":
            cpg_count+=1
    if ref[0]=="G" and ref_prev_base=="C":
        cpg_count+=1  
    if ref[-1]=="C" and ref_following_base=="G":
        cpg_count+=1   
    return cpg_count          


def methylation_status(read, ref, flag, following_ref_base):
    if len(ref) != len(read):
        raise ValueError("Both sequences must be of the same length")
    
    meth_status = len(read)*[2]
    for i in range(len(ref) - 1):
        if ref[i:i+2] == "CG":
            if ((flag & 99 == 99 or flag & 147 == 147) and read[i:i+2] == "TG"): 
                meth_status[i] = 1                
            elif ((flag & 83 == 83 or flag & 163 == 163) and read[i:i+2]== "CA"):
                meth_status[i] = 1 
            else:
                meth_status[i] = 0 
    if following_ref_base=='G' and ref[-1]=='C' and (flag & 99 == 99 or flag & 147 == 147):
        if read[-1]=='T':
            meth_status[-1]=1
        elif read[-1]=='C':
            meth_status[-1]=0     
    return meth_status


