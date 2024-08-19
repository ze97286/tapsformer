import argparse
import pandas as pd
import pysam
import os
from tapsformer.file_constants import *
from tapsformer.differential_methylation.utils import *

def process_chr(chr, bam_file, tumour_tissue, sample, dmr_dir):
    bam = pysam.AlignmentFile(bam_file)
    pileup_region = load_next_pileup_region2(bam, chr)
    df = pd.read_csv(os.path.join(TAPSFORMER_DATA_BY_CPG, dmr_dir, "hypomethylated_dmrs.bed"), sep='\t')
    df['dmr_id']=len(df)-df.index
    # filter only very strong hypomethylation regions
    df = df[df.hypomethylation_strength=='Very Strong']    
    df = df[(df.chr == chr)]
    df = df.sort_values(['start'])
    print(f"{len(df)} DMRs with hypomethylation")
    if len(df)==0:
        return        
    
    # Initialize lists to store data
    cpgs_count, mods, unmods, reads, refs, mstates = [], [], [], [], [], []
    prevs, followings, region_starts, region_ends, starts, ends, flags, qn, first_in_pair = [], [], [], [], [], [], [], [], []
    dmr_ids,mapqs, alignment_scores = [],[],[]

    ref_genome = pysam.FastaFile(GENOME_FILE)
    pileupcolumn = next(pileup_region, None)

    seen_id = set()
    
    for _, r in df.iterrows():
        start, end = r['start'], r['end']
        dmr_id = r['dmr_id']
        print(f"Processing DMR {start}-{end}, so far seen {len(reads)} reads")
        
        if pileupcolumn is None:
            break  # End of file
        
        while pileupcolumn is not None and pileupcolumn.pos <= end:
            while pileupcolumn.pos < start:
                pileupcolumn = next(pileup_region, None)
                if pileupcolumn is None or pileupcolumn.pos > end:
                    break
            
            if pileupcolumn is None or pileupcolumn.pos > end:
                break
            
            for pileupread in pileupcolumn.pileups:
                read = pileupread.alignment
                id = read.query_name + "_" + str(read.is_read1)
                if id in seen_id:
                    continue 
                
                if "I" in read.cigarstring or "D" in read.cigarstring:
                    continue 

                reconstructed_read, reconstructed_ref, ref_prev_base, ref_following_base = align_read_to_reference(chr, read, ref_genome)
                if reconstructed_read is None:
                    continue
                
                seen_id.add(id)

                skips_start = 0
                skips_end = 0
                if read.cigartuples[0][0]==4:
                    skips_start=int(read.cigartuples[0][1])
                if read.cigartuples[-1][0]==4:
                    skips_end=int(read.cigartuples[-1][1])

                cpg_count = count_cpgs(reconstructed_ref, ref_prev_base, ref_following_base)
                if cpg_count <=0:
                    continue

                mstatus = methylation_status(reconstructed_read, reconstructed_ref, read.flag, ref_following_base)
                
                print("skips_start",skips_start)
                print("skips_end",skips_end)
                print("len(reconstructed_ref)-skips_end",len(reconstructed_ref)-skips_end)
                
                reconstructed_ref = reconstructed_ref[skips_start,len(reconstructed_ref)-skips_end]
                mstatus=mstatus[skips_start,len(mstatus)-skips_end]
                reconstructed_read = reconstructed_read[skips_start,len(reconstructed_read)-skips_end]
                alignment_scores.append(read.get_tag("AS"))
                dmr_ids.append(dmr_id)
                mods.append(mstatus.count(1))
                unmods.append(mstatus.count(0))
                mapqs.append(read.mapping_quality)
                reads.append(reconstructed_read)
                prevs.append(ref_prev_base)
                followings.append(ref_following_base)
                refs.append(reconstructed_ref)
                cpgs_count.append(cpg_count)
                region_starts.append(start)
                region_ends.append(end)
                starts.append(read.reference_start)
                ends.append(read.reference_end)
                flags.append(read.flag)
                qn.append(read.query_name)
                first_in_pair.append(read.is_read1)
            
            pileupcolumn = next(pileup_region, None)
            if pileupcolumn is None or pileupcolumn.pos > end:
                break
    
    data_dict = {
        "dmr_id":dmr_ids,
        "methylation_state":mstates,
        "read": reads, 
        "flag": flags,
        "reference": refs, 
        "mod": mods,
        "unmod": unmods, 
        "informative_cpg": starts, 
        "ref_prev_base": prevs,
        "ref_following_base": followings,
        "region_start": region_starts, 
        "region_ends": region_ends,
        "cpgs": cpgs_count,
        "qn": qn,
        "first_in_pair": first_in_pair,
        "mapq":mapqs,
        "alignment_score":alignment_scores,        
    }
    
    df = pd.DataFrame(data_dict)
    df['label'] = 1 if tumour_tissue else 0
    df = df.drop_duplicates(['qn', 'first_in_pair'])
    
    out_dir = os.path.join(METHYLATION_SAMPLES_OUT, dmr_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    df.to_parquet(os.path.join(out_dir, f"{sample}_{chr}_reads.parquet"), index=False)
    print(f"Finished processing chromosome {chr} for sample {sample} with {len(df)} reads")
    bam.close()
    return df

def main():
    """
    Extract all reads from hypomethylated DMRs with 3 or more CpGs for the given chromosome and bin size.
    """
    parser = argparse.ArgumentParser(description="TAPSFORMER")
    parser.add_argument("--bam", dest="bam", help="Input BAM file", required=True)
    parser.add_argument("--chr", dest="chr", help="Chromosome", required=True)
    parser.add_argument("--sample", dest="sample", help="Output sample", required=True)
    parser.add_argument("--true-labels", dest="tissue", action=argparse.BooleanOptionalAction, help="Generate true labels", required=False, default=False)
    parser.add_argument("--dmrs", dest="dmr_dir", help="DMR directory", required=True)
    
    inputs = parser.parse_args()
    tumour_tissue = inputs.tissue
    bam_file = inputs.bam
    sample = inputs.sample
    chr = inputs.chr
    dmr_dir = inputs.dmr_dir
    
    process_chr(chr, bam_file, tumour_tissue, sample, dmr_dir)

if __name__ == "__main__":
    main()
