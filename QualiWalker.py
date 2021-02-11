# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 09:48:34 2015

@author: hilmar

QualiWalker
Hilmar Berger 
(C) Max Planck Institute for Infection Biology, 2015
"""

import os, string
from WindowWalker import WindowWalker, RefSeqStore
import pysam
from VcfReader import VcfReader
from vcf_template import vcf_template
import numpy as np
import pandas

def generate_per_window_characteristics(inp_file, ref_file, window_size, regions=None,  known_polymorphic_positions=None, call_variants = True, max_reads_in_bin = 1e5):

    ref_seq_segment_size = 1e6
    
    ff = pysam.AlignmentFile(inp_file,"rb")
    ref_seq_store = RefSeqStore(ref_file, ref_seq_segment_size)
    
    # loop through all contigs provided in the BAM file
    if len(regions)==0:
        regions = [(s,None,None) for s in ff.references]

    all_window_descriptions = None
    all_variants = {}
    total_read_cnt = 0
    accepted_read_cnt = 0
    excluded_read_cnt = 0
    bin_cnt = 0

    window_size_default = window_size
    
    # loop through all regions. Start and end positions are 1-based and inclusive. 
    for curr_region in regions:
        seq = curr_region[0]
        start_pos = curr_region[1]
        end_pos = curr_region[2]

        if start_pos is None:
            start_pos = 1
            end_pos = ref_seq_store.get_ref_seq_length(seq)

        else:
            if start_pos < 1 or start_pos > ref_seq_store.get_ref_seq_length(seq):
                raise Exception, "Invalid start position [%d]" % start_pos
            if end_pos < start_pos or end_pos > ref_seq_store.get_ref_seq_length(seq):
                raise Exception, "Invalid end position [%d]" % end_pos

        window_size = min(window_size_default, end_pos - start_pos + 1)        

        bin_starts = np.arange(start_pos, end_pos+1, window_size)
        bin_ends = np.minimum(end_pos, bin_starts + window_size-1)
        bin_num = len(bin_starts)
        result_dict = {"seq":np.array([seq]*bin_num,dtype=object), "window_start":bin_starts, "window_end":bin_ends}
        result_dict["secondary_alignments"] = np.array([np.nan] * bin_num, dtype=float)
        result_dict["chimeric_alignments"] = np.array([np.nan] * bin_num, dtype=float)

        w_ind = 0
        #window_start, window_end = result_df.iloc[w_ind, 1:3]

        window_start = start_pos
        window_end = min(window_start + window_size-1, end_pos)

        previous_win_store = WindowWalker(window_size, 0, 10, "NNNNNNNNNN")

        # loop through all genomic windows
        while window_start < end_pos:
            secondary_mappings = 0
            chimeric_read_mappings = 0
            
            #print "Processing window: %s\t%d\t%d" % (seq, window_start, window_end)
            # get ref segment for window
            ref_seq_segment = ref_seq_store.get_segment(seq, window_start, window_end)
            
            snp_pos_list = {}
            # populate list of known SNP/Indel positions
            if not known_polymorphic_positions is None:
                snp_pos = known_polymorphic_positions.get_block(seq, window_start, window_end)
                for e in snp_pos:
                    # positions read from VCF are 0 based - we have to convert to 1-based here
                    # NOTE: Do not expect the ref allele from e.g. dbSNP or other genotype files to correspond to the genome reference sequence
                    # SNV are not stored in that way, so the only thing you can reasonably expect is that the ref sequence matches either ref or alt
                    # FIXME: This does not handle InDels correctly
                    snp_pos_list[e.pos+1] = e
            # initialize window store
            win_store =  WindowWalker(window_start, window_end, ref_seq_segment, known_polymorphic_positions = snp_pos_list)
            bin_cnt += 1
            
            # add all reads in window to store
            # note that pysam API uses 0-based half-open coordinates
            reads_in_bin = 0
            rr=ff.fetch(seq, window_start-1, window_end)
            for r in rr:
                total_read_cnt += 1
                reads_in_bin += 1
                if reads_in_bin > max_reads_in_bin:
                    break
                
                # FIXME: Secondary hits might be informative regarding quality of mappability in a window. Consider including them accordingly.
                if not (r.is_unmapped or r.is_secondary or r.is_supplementary):
                    try:
                        if r.reference_start+1 < window_start:
                            read_vec, read_id = previous_win_store.get_reads_from_store(r)
                            if read_id is not None:
                                win_store.insert_read_vector(read_id, read_vec)
                            else: 
                                win_store.add_read(r)
                                accepted_read_cnt += 1

                        else:
                            win_store.add_read(r)
                            accepted_read_cnt += 1
                            
                    except Exception, e:
                        print "Could not add read [%s], read # %d" % (r.query_name, total_read_cnt)
                        raise e
                else:
                    if (r.is_secondary):
                        secondary_mappings += 1
                    if (r.is_supplementary):
                        chimeric_read_mappings += 1
                        
                    excluded_read_cnt += 1
                    
            # get metrics for this window
            tmp = win_store.describe_window()
            for k in tmp.keys():
                try:
                    result_dict[k][w_ind] = tmp[k]
                except KeyError:
                    result_dict[k]= np.array([np.nan] * bin_num, dtype=float)
                    result_dict[k][w_ind] = tmp[k]

            result_dict["secondary_alignments"][w_ind] = secondary_mappings
            result_dict["chimeric_alignments"][w_ind] = chimeric_read_mappings

            if call_variants:
                win_store.call_variants()
                all_variants[(seq, window_start)] = win_store.get_variants()

            window_start = window_start + window_size
            window_end = min(window_start + window_size, end_pos)
            #print "->:", window_start, window_end
            previous_win_store = win_store
            del win_store
            w_ind += 1

        # FIXME: make sure that final table will have all columns in case a sub-table has less columns
        result_df = pandas.DataFrame(result_dict)
        #result_df.set_index(['seq','window_start'], inplace=True, drop=False)
        if(all_window_descriptions is None):
            all_window_descriptions = result_df
        else:
           all_window_descriptions = all_window_descriptions.append(result_df)

    print "Total reads seen: %d\tUnique reads: %d\tExcluded secondary alignments: %d" % (total_read_cnt, accepted_read_cnt, excluded_read_cnt)

    all_window_descriptions.set_index(['seq','window_start','window_end'], inplace=True, drop=False)
    
    return all_window_descriptions, all_variants

###############################################################################    
if __name__ == '__main__':

    import argparse, sys
#    import time
    from operator import itemgetter
    
    parser = argparse.ArgumentParser(description='Generate overview in genomic bins')
    parser.add_argument('input_file', type=str, help='Input BAM file')

    parser.add_argument('-f','--reference-fasta', action='store', type=str, default=None,
                       help='Path to reference genome FASTA file')

    parser.add_argument('-w','--window-size', action='store', type=int, default=300,
                       help='Genomic window size. Defaults to 300.')

    parser.add_argument('-r','--region', action='store', type=str, default=None,
                       help='Restrict analysis to genomic region (chr:start-end')

    parser.add_argument('-R','--region-file', action='store', type=str, default=None,
                       help='Restrict analysis to genomic regions in bed format (chr start end) file')

    parser.add_argument('-o','--output-file', action='store', type=str, default=None,
                       help='Output file. Default: STDOUT')

    parser.add_argument('-s','--snp-positions', action='store', type=str, default=None,
                       help='Indexed VCF holding positions of known polymorphic sites (e.g. dbSNP).')
                       
    parser.add_argument('-c','--called_variants_file', action='store', type=str, default=None,
                        help='File name for variants detected in sequences after filtering out bad reads and known SNP positions. If not given, variants will not be called.')

    parser.add_argument('-m','--max_reads_per_bin', action='store', type=int, default=1e6,
                        help='Maximum number of reads to be considered per bin. Excess reads in affected bins will be ignored. Default 1e6')

    args = parser.parse_args()

    inp_file = args.input_file
    ref_file = args.reference_fasta

    known_polymorphic_positions = None
    if not args.snp_positions is None:
        known_polymorphic_positions = VcfReader(args.snp_positions)

    call_variants = True
    if args.called_variants_file is None:
        call_variants = False

    regions = []
    # Regions from BED file    
    # FIXME: We do not check here if regions are non-overlapping (which they should; otherwise window routines might fail)
    if not args.region_file is None:
        with open(args.region_file, "r") as region_file:
            for lr in region_file:
                fields = lr.strip("\n\r").split("\t")
                if len(fields) < 3: continue
                seq = fields[0]
                start = int(fields[1]) + 1
                end = int(fields[2])
                regions.append((seq, start, end))
    else: 
        if not args.region is None: # Regions from command line
            f1 = args.region.split(":")
            seq = f1[0]
            if len(f1) > 1:
                start, end = [int(e) for e in f1[1].split("-")]
                regions.append((seq, start, end))
            else:
                start = None
                end = None
                regions.append((seq, start, end))

    max_reads_per_bin = 1 if args.max_reads_per_bin < 1 else args.max_reads_per_bin

    # generate results
    res, variants = generate_per_window_characteristics(inp_file, ref_file, args.window_size, regions, known_polymorphic_positions, call_variants, max_reads_per_bin )
    columns = res.columns

    bin_index_cols = ["seq", "window_start","window_end"]
    
#    # Write results
#    if args.output_file is None:
#        ofile = sys.stdout
#    else:
#        # make sure that we have a sufficiently large output buffer to improve perforrmance when writing over NFS or similar
#        ofile = open(args.output_file, "w", 10**7)
#
    columns = sorted(columns)
    columns = [e for e in columns if not e in bin_index_cols]
#    first = True
#    
#    for k in res.index:
#        if first:
#            first = False
#            headers = bin_index_cols + columns
#            ofile.write( "\t".join(headers) + "\n")
#        
#        ofile.write("%s\t%d\t%d" % k)
#        if res.loc[k,"read_number"]>0:
#            for c in columns:        
#                ofile.write("\t%s" % str(res.at[k,c]))
#            ofile.write("\n")
#        else:
#            for c in columns:
#                if c=="read_number":
#                    ofile.write("\t0")
#                else:
#                    ofile.write("\tNA")
#            ofile.write("\n")
#            
#    if args.output_file is not None:
#        ofile.close()

    res.reindex(sorted(columns), axis=1)
    res.to_csv(args.output_file, sep="\t", index=True, na_rep="NA", columns = columns)

    list_to_str = lambda x: ",".join([str(e) for e in x])

    # write variants to VCF
    if call_variants:       
        vcf_ofile = open(args.called_variants_file+".vcf", 'w')
        vcf_ofile.write(string.Template(vcf_template).substitute({"sample_name":os.path.basename(inp_file), "ref_file":ref_file}))
        format_str="DP:DV:DPR:RFSP:R1FSP:MQ:BQ:POS"
        for rr in sorted(variants, key=itemgetter(0,1)):
            for vpos in sorted(variants[rr].keys()):
                vv = variants[rr][vpos]
                var_alleles = vv["var"].keys()
                ref = vv["ref"]
                allele_cnts = [vv["var"][e]["allele_cnt"] for e in var_alleles]
                r1fsp = [vv["var"][e]["read_plus_prop"] for e in var_alleles]
                r1p = [vv["var"][e]["read1_plus_prop"] for e in var_alleles]
                mq = [vv["var"][e]["mq_median"] for e in var_alleles]
                bq = [vv["var"][e]["bq_median"] for e in var_alleles]
                pos = [vv["var"][e]["pos_median"] for e in var_alleles]
                info_str = "DP=%d" % ref["cov_total"]
                non_ref_cnt = sum(allele_cnts)
                total_cov_hq = non_ref_cnt + vv["ref"]["allele_cnt"]
                sample_str = "%d:%d:%s:%s:%s:%s:%s:%s" % (total_cov_hq,non_ref_cnt, list_to_str([total_cov_hq-non_ref_cnt] + allele_cnts), list_to_str(r1fsp), list_to_str(r1p), list_to_str(mq), list_to_str(bq), list_to_str(pos))
                output_str="%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n" % (rr[0],vpos, ".", ref["base"], list_to_str(var_alleles),0,".",info_str, format_str, sample_str )
                vcf_ofile.write(output_str)
        vcf_ofile.close()
        
