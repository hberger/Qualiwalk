# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 09:48:34 2015

@author: hilmar
(C) Max Planck Institute for Infection Biology, 2015
"""

import pysam
import numpy as np
import string 
import bottleneck as bn

from seq_to_vector import seq2num_vec, set_translation

set_translation()

nuc_to_code_trans = string.maketrans("ACGTNUWSMKRYBDHV","1234545555555555")
code_to_nuc_trans = string.maketrans("12345","ACGTN")

class WindowWalker:
    """stores read information from a single genomic window and computes per window metrics
       genomic window is modeled as two [total_read_number, window_length] integer arrays (one for seq, one for quality)
       # Code alignment result as numeric vector in the window
       # -1 - (soft-) clipped 
       # 0 - not aligned
       # 1..4 - aligned, read base (1-A, 2-C, 3-G, 4-T)
       # 5 ... aligned, read base N
       # 6 - aligned, Del
       # for each read also store strand, R1/R2, number and position (relativ to read and genome) of inserts    
       
       start and end are both inclusive and refer to a 1-based coordinate system for all functions. 
       Note that pysam API calls operate on 0-based half-open intervals instead.
    """
    # FIXME: Define start/end interval correctly (i.e. [] vs.  [) )
    def __init__(self, start, end, ref_seq, alignment_len_default = 200, known_polymorphic_positions = None, insert_size_mean = 300, insert_size_sd=50):
        self.__read_hash = {}
        self.window_length = end-start+1
        self.window_start = start
        self.window_end = end
        #print start, end, self.window_length
        self.window_ref_seq = ref_seq
        self.alignment_len_default = alignment_len_default
        self.window_array_seq = None
        self.window_array_qual = None
        self.read_characteristics = None
        self.window_description = {}
        self.window_positional_stats = {}
        self.trans_tab = nuc_to_code_trans
        # FIXME: Currently unused. Probably we wanted to enable an a-priori filtering of reads with abnormal insert sizes. 
        self.insert_size_min_thr = insert_size_mean - 3 * insert_size_sd
        self.insert_size_max_thr = insert_size_mean + 3 * insert_size_sd
        if not known_polymorphic_positions is None:
            self.known_snp_pos = np.array([e-start for e in known_polymorphic_positions.keys()], dtype=int)
        else:
            self.known_snp_pos = None
        self.variants = {}
        
    def add_read(self, read):
        """add the converted sequence, quality and additional read characteristics to the store"""
        #vv = self.__seq2num_vec(read.query_sequence, read.query_qualities, read.reference_start+1, read.cigartuples   )
        #print ">%s<" % read.query_sequence
        try:
            vv = seq2num_vec(read.query_sequence, read.query_qualities, read.reference_start+1, read.cigartuples, self.alignment_len_default )
        except Exception, e:
            print
            raise e
            
        if vv is None:
            raise Exception, "Read alignment invalid"
            
        mate_num = 1 if read.is_read1 else 2
        read_strand = "+" if not read.is_reverse else "-"
        first_read_strand = None
        if mate_num == 1:
            first_read_strand = read_strand
        else:
            first_read_strand = "+" if not read.mate_is_reverse else "-"
            
        proper_pair = True if read.is_paired and read.is_proper_pair else False
        # mates are considered invalid if unmapped or on a different contig/chromosome
        valid_chrom_mate = False if (read.mate_is_unmapped or (read.next_reference_id != read.reference_id)) else True
        template_length = read.template_length
        mq = read.mapping_quality
        read_id = (read.query_name, mate_num)
        real_start = read.reference_start+1-vv[4] if read_strand=="+" else read.reference_start+1-vv[4]+len(read.query_sequence)
        read_characteristics = { "read_id": read_id, "cigartuples": read.cigartuples, "mate_num": mate_num, "read_strand": read_strand, 
                                "first_read_strand": first_read_strand, "proper_pair": proper_pair, "read_start": read.reference_start+1, 
                                "indels": vv[2], "alignment_length": vv[3], "start_position_offset": vv[4], "mq": mq , "5prime_start": real_start,
                                "insert_size":template_length, "valid_chrom_mate":valid_chrom_mate }
                                
        self.__read_hash[read_id] = (read_characteristics, vv[0], vv[1])

    def insert_read_vector(self, read_id, read_vec):
        """This method allows to reuse data from reads overlapping the window boundary from a previous window"""
        self.__read_hash[read_id] = read_vec

    def get_reads_from_store(self, read=None):
        if read is None:
            return self.__read_hash
        mate_num = 1 if read.is_read1 else 2
        read_id = (read.query_name, mate_num)
        if read_id in self.__read_hash:
            return self.__read_hash[read_id], read_id
        else:
            return None, None
        
    def get_window_array(self):
        if self.window_array_seq is None:
            self.__generate_window_arrays()
        return (self.window_array_seq, self.window_array_qual, self.read_characteristics)


    def prepare_window_metrics(self, bq_threshold=20, low_qual_bq_thr = 10):
        if self.window_array_seq is None:
            self.__generate_window_arrays()

        #tmp = np.hstack((self.window_array_qual.transpose(), self.window_array_seq.transpose()) )
        #np.savetxt("qual_seq_array.txt", tmp, delimiter='\t', fmt='%d')
        
        read_num = self.window_array_seq.shape[0]

        if read_num != 0:
            # all bases in reads, including soft-clipping
            self.all_bases = (self.window_array_seq != 0)
            # Coverage: Mask unmapped, N and clipped bases, count everything else as aligned
            self.aligned_mask = np.logical_or(self.window_array_seq < 1, self.window_array_seq ==5)
            # Exclude dels from SNV 
            self.snv_mask = np.logical_or(self.aligned_mask, self.window_array_seq == 6)
            # Exclude low quality bases
            self.quality_mask = self.window_array_qual < bq_threshold
            # translate ref seq to numeric code
            self.ref_code = np.array([int(e) for e in string.translate(self.window_ref_seq, self.trans_tab)], dtype=int)
            # identify positions in ref seq that are not determined (i.e. N bases in ref) and exclude them
            self.genomic_N_pos = np.zeros_like(self.window_array_seq)
            self.genomic_N_pos[:,(self.ref_code == 5)] = 1
            # define the array of non-N-ref positions with valid (non-del, aligned) bases
            self.valid_ref_snv_mask = np.logical_or(self.snv_mask, self.genomic_N_pos)

            self.valid_pos_seq_array = np.copy(self.window_array_seq)
            self.valid_pos_qual_array = np.copy(self.window_array_qual)

            # determine differences between ref seq and valid positions
            self.snv = np.array((self.valid_pos_seq_array.transpose() - self.ref_code[:, np.newaxis]).transpose() != 0, dtype=float)
            self.snv[self.valid_ref_snv_mask] = np.NaN

            self.valid_pos_seq_array[self.valid_ref_snv_mask] = np.NaN
            self.valid_pos_qual_array[self.valid_ref_snv_mask] = np.NaN
            
            # compute coverage per position in the window
            self.coverage_per_pos = bn.nansum(np.logical_not(self.aligned_mask).astype(int), axis=0)
            # some reads have higher number of mismatches - identify them
            self.mismatch_reads = bn.nansum(self.snv, axis=1)>2
            self.mq_per_read = np.array([e["mq"] for e in self.read_characteristics], dtype=int)
            
            self.low_qual_bases =  (self.window_array_qual < low_qual_bq_thr)
            
            # quality values 
            #self.read_median_bq = bn.nanmedian(self.valid_pos_qual_array, axis=1)
            self.read_mean_bq = bn.nanmean(self.valid_pos_qual_array, axis=1)
            self.read_min_bq = bn.nanmin(self.valid_pos_qual_array, axis=1)
            self.read_lowqual_bases = np.logical_and(np.logical_not(self.valid_ref_snv_mask), self.low_qual_bases)
            self.read_lowqual_base_count = bn.nansum(self.read_lowqual_bases.astype(int), axis=1)
            # identify positions that contain mismatches
            self.snv_pos = bn.nansum(self.snv, axis=0)>0

            # define array of high-quality bases and compute SNV on them
            self.valid_ref_snv_good_qual_mask = np.logical_or(self.valid_ref_snv_mask, self.quality_mask)
            self.snv_bq_filtered = np.copy(self.snv)
            self.snv_bq_filtered[self.valid_ref_snv_good_qual_mask] = np.NaN
            self.snv_pos_bq_filtered = bn.nansum(self.snv_bq_filtered, axis=0)>0
            
            # quality of aligned positions only
            self.qual_aligned_mask = np.logical_or(self.window_array_seq==0, self.genomic_N_pos)
            self.qual_aligned_lowqual_bases = np.logical_and(np.logical_not(self.qual_aligned_mask), self.low_qual_bases)
            
            # from array of high-quality bases also compute "bad" reads and compute SNV on the rest. 
            self.good_read_mask = np.zeros_like(self.window_array_seq)
            self.good_read_mask[np.logical_or(np.logical_or(self.read_lowqual_base_count > 2,self.mismatch_reads), self.mq_per_read==0),:] = 1
            
            # SNV from high-quality aligned bases and high-quality reads
            self.valid_ref_snv_good_qual_good_reads_mask = np.logical_or(self.good_read_mask, self.valid_ref_snv_good_qual_mask )
            self.snv_bq_filtered_good_reads = np.copy(self.snv)
            self.snv_bq_filtered_good_reads[self.valid_ref_snv_good_qual_good_reads_mask] = np.NaN
            self.snv_pos_bq_filtered_good_reads = bn.nansum(self.snv_bq_filtered_good_reads, axis=0)>0
            
            # SNV from high-quality aligned bases and high-quality reads with known SNP positions excluded
            self.known_snp_mask = np.zeros_like(self.window_array_seq)
            self.known_snp_mask[:,self.known_snp_pos] = 1

            self.valid_ref_snv_good_qual_good_reads_no_known_site_mask = np.logical_or(self.valid_ref_snv_good_qual_good_reads_mask, self.known_snp_mask)
            self.snv_bq_filtered_good_reads_no_snp = np.copy(self.snv)
            self.snv_bq_filtered_good_reads_no_snp[self.valid_ref_snv_good_qual_good_reads_no_known_site_mask] = np.NaN
            self.snv_pos_bq_filtered_good_reads_no_snp = bn.nansum(self.snv_bq_filtered_good_reads_no_snp, axis=0)>0

            self.pos_median_bq = bn.nanmedian(self.valid_pos_qual_array, axis=0)
            self.pos_min_bq = bn.nanmin(self.valid_pos_qual_array, axis=0)
            self.pos_lowqual_base_count = bn.nansum(self.read_lowqual_bases.astype(int), axis=0)

            # note that for now we only count start positions
            # FIXME: CHeck if we can code this in the windows sequence matrix, possibly using an bit-overlay matrix or similar (i.e. code insertion flanking nucleotides with different higher bits)
            self.window_ins_positions = {}
            self.window_del_positions = {}
            tmp_insert_sizes = {}
            
            for e in self.read_characteristics:
                insertions = e["indels"]["I"]
                deletions = e["indels"]["D"]
                rstart = e["read_start"]+e["start_position_offset"]

                for ii in insertions.keys():
                    genomic_pos = ii+rstart
                    if genomic_pos >= self.window_start and genomic_pos <= self.window_end:
                        try:
                            self.window_ins_positions[genomic_pos] += 1
                        except KeyError:
                            self.window_ins_positions[genomic_pos] = 1

                for ii in deletions.keys():
                    genomic_pos = ii+rstart
                    if genomic_pos >= self.window_start and genomic_pos <= self.window_end:
                        try:
                            self.window_del_positions[genomic_pos] += 1
                        except KeyError:
                            self.window_del_positions[genomic_pos] = 1
                
                # make sure we consider insert sizes only once per fragment in the same window
                tmp_insert_sizes[e["read_id"][0]] = np.abs(e["insert_size"]) if e["insert_size"]<>0 else np.nan

            self.insert_sizes_unique = np.array(tmp_insert_sizes.values(),dtype=float)
            
        self.window_description["read_number"] = read_num
#        print "Mean BQ by read %s" % str(bn.nanmean(self.valid_pos_qual_array, axis=1))
#        print "Mean BQ %s" % str(bn.nanmean(bn.nanmean(self.valid_pos_qual_array, axis=1)))
#        print "Mean BQ v2 %s" % str(bn.nanmean(self.valid_pos_qual_array))

    def describe_window(self, call_variants = True, bq_threshold=20, low_qual_bq_thr = 10):
        # prepare window metrics    
        if not "read_number" in self.window_description:
            self.prepare_window_metrics(bq_threshold, low_qual_bq_thr)
        read_num = self.window_description["read_number"]
        
        if read_num != 0:
            self.window_description["cov_median"] = bn.nanmedian(self.coverage_per_pos)
            self.window_description["cov_mean"] = bn.nanmean(self.coverage_per_pos)
            self.window_description["cov_q1"] = np.percentile(self.coverage_per_pos, 25)
            self.window_description["cov_p10"] = np.percentile(self.coverage_per_pos, 10)
            self.window_description["cov_0"] = bn.nansum(self.coverage_per_pos==0)
            self.window_description["snv_total"] = bn.nansum(self.snv)
            self.window_description["snv_positions"] = bn.nansum(self.snv_pos.astype(int))
            self.window_description["snv_total_bq_filtered"] = bn.nansum(self.snv_bq_filtered)
            self.window_description["snv_positions_bq_filtered"] = bn.nansum(self.snv_pos_bq_filtered.astype(int))
            # FIXME: The name of this metric is somewhat misleading, since we completely ignore SNV outside the current window. This should be documented
            self.window_description["reads_gt2_snv"] = bn.nansum(self.mismatch_reads.astype(int))
            self.window_description["mq_median"] = bn.nanmedian(self.mq_per_read)
            self.window_description["mq_mean"] = bn.nanmean(self.mq_per_read)
            self.window_description["read_plus_strand_proportion"] = float(bn.nansum(np.array([e["read_strand"]=="+" for e in self.read_characteristics], dtype=int)))/read_num
            self.window_description["read1_plus_strand_proportion"] = float(bn.nansum(np.array([e["first_read_strand"]=="+" for e in self.read_characteristics], dtype=int)))/read_num
            self.window_description["R1_proportion"] = float(bn.nansum(np.array([e["mate_num"]==1 for e in self.read_characteristics], dtype=int)))/read_num
            self.window_description["insertion_num"] = bn.nansum(np.array(self.window_ins_positions.values(), dtype=int))
            self.window_description["deletion_num"] = bn.nansum(np.array(self.window_del_positions.values(), dtype=int))
            self.window_description["ins_positions"] = len(self.window_ins_positions.keys())
            self.window_description["del_positions"] = len(self.window_del_positions.keys())
            self.window_description["LowQual_bases"] = bn.nansum(self.qual_aligned_lowqual_bases.astype(int))
            self.window_description["LowQual_gt2_reads"] = bn.nansum((self.read_lowqual_base_count > 2).astype(int))
            self.window_description["median_bq"] = bn.nanmedian( self.valid_pos_qual_array )
            self.window_description["mean_bq"] = bn.nanmean( self.read_mean_bq )
            self.window_description["avg_read_min_bq"] = bn.nanmean( self.read_min_bq )
            self.window_description["snv_total_bq_filtered_good_reads"] = bn.nansum(self.snv_bq_filtered_good_reads)
            self.window_description["snv_positions_bq_filtered_good_reads"] = bn.nansum(self.snv_pos_bq_filtered_good_reads.astype(int))
            self.window_description["snv_total_bq_filtered_good_reads_no_snp"] = bn.nansum(self.snv_bq_filtered_good_reads_no_snp)
            self.window_description["snv_positions_bq_filtered_good_reads_no_snp"] = bn.nansum(self.snv_pos_bq_filtered_good_reads_no_snp.astype(int))
            self.window_description["bases_aligned"] = bn.nansum(np.logical_not(self.snv_mask).astype(int))
            self.window_description["soft_clipped_bases"] = bn.nansum((self.window_array_seq==-1).astype(int))
            self.window_description["total_base_cnt"] = bn.nansum((self.all_bases).astype(int))
            self.window_description["improper_pairs"] = bn.nansum(np.array([e["proper_pair"]==False for e in self.read_characteristics], dtype=int))
            self.window_description["median_insert_size"] = bn.nanmedian(self.insert_sizes_unique) 
            self.window_description["mean_insert_size"] = bn.nanmean(self.insert_sizes_unique) 
            self.window_description["invalid_chrom_mates"] = bn.nansum(np.array([e["valid_chrom_mate"]==False for e in self.read_characteristics], dtype=int))

        return self.window_description
        
    def get_variants(self):
        return self.variants

    def __count_nuc_per_position(self, mat):
        nuc_counts = {}
        nuc_counts["A"] = bn.nansum((mat==1).astype(int), axis=0)
        nuc_counts["C"] = bn.nansum((mat==2).astype(int), axis=0)
        nuc_counts["G"] = bn.nansum((mat==3).astype(int), axis=0)
        nuc_counts["T"] = bn.nansum((mat==4).astype(int), axis=0)
        return nuc_counts
        
    def describe_positions(self, bq_threshold=20, low_qual_bq_thr = 10):
        """
        extract per-position metrics for all base positions in the window
        Metrics:
          - total coverage (all aligned bases)
          - total coverage (all bases, including soft-clipped)
          - total cov of high quality bases
          - median BQ per pos
          - min BQ per pos
          - number of low BQ bases per pos
          - counts for A,T,C,G total per pos
          - counts for A,T,C,G in HQ bases per pos
          - N count per pos
          - insertion count in adjacent bases (+/-1) per pos
          - del count per pos
        """
        
        # prepare metrics
        if not "read_number" in self.window_description:
            self.prepare_window_metrics(bq_threshold, low_qual_bq_thr)
        read_num = self.window_description["read_number"]

        self.window_positional_stats["position"] = xrange(self.window_start, self.window_end+1)
        if read_num != 0:
            self.window_positional_stats["cov_aligned"] = self.coverage_per_pos
            self.window_positional_stats["cov_total"] = bn.nansum(self.all_bases.astype(int), axis=0)
            #FIXME: We here also exclude N positions in genome and deletions - check if that is what we really want
            self.window_positional_stats["cov_good_bases"] = bn.nansum(np.logical_not(self.valid_ref_snv_good_qual_mask).astype(int),axis=0) # does not count DELs 
            self.window_positional_stats["median_BQ"] = bn.nanmedian(self.valid_pos_qual_array, axis=0) # this does exclude DELs and non-aligned bases
            self.window_positional_stats["min_BQ"] = bn.nanmedian(self.valid_pos_qual_array, axis=0) # this does exclude DELs and non-aligned bases]
            self.window_positional_stats["low_BQ_count"] = bn.nansum(self.read_lowqual_bases.astype(int), axis=0)
            positional_nuc_cnts = self.__count_nuc_per_position(self.valid_pos_seq_array)
            self.window_positional_stats["Nuc_A_cnt_all"] = positional_nuc_cnts["A"]
            self.window_positional_stats["Nuc_T_cnt_all"] = positional_nuc_cnts["T"]
            self.window_positional_stats["Nuc_C_cnt_all"] = positional_nuc_cnts["C"]
            self.window_positional_stats["Nuc_G_cnt_all"] = positional_nuc_cnts["G"]
            tmp_seq_array = np.copy(self.valid_pos_seq_array)
            tmp_seq_array[self.valid_ref_snv_good_qual_good_reads_mask] = np.NaN
            positional_hq_nuc_cnts = self.__count_nuc_per_position(tmp_seq_array)
            self.window_positional_stats["NucA_cnt_HQ"] = positional_hq_nuc_cnts["A"]
            self.window_positional_stats["NucT_cnt_HQ"] = positional_hq_nuc_cnts["T"]
            self.window_positional_stats["NucC_cnt_HQ"] = positional_hq_nuc_cnts["C"]
            self.window_positional_stats["NucG_cnt_HQ"] = positional_hq_nuc_cnts["G"]
            self.window_positional_stats["Del_Cnt"] =  bn.nansum((self.window_array_seq == 6).astype(int), axis=0)
            self.window_positional_stats["N_Cnt"] =  bn.nansum((self.window_array_seq == 5).astype(int), axis=0)
            
            self.window_insertions = np.zeros((read_num, self.window_length))
            row_index = 0
            for rc in self.read_characteristics:
                # For soft-clipped reads, the read start defines the start of the matching sequence, i.e. after the leftmost
                # clipped sequence. We therefore have to adjust for this difference and take the clipped sequence into account 
                # in the same way as any other CIGAR op
                read_start = rc["read_start"] - rc["start_position_offset"] # adjust offset for soft-clipped reads
                w_offset = max(0,read_start-self.window_start) 
                r_offset = max(0,self.window_start - read_start)
    
                ins_events = rc["indels"]["I"]
                for qpos in ins_events.keys():
                    window_ins_pos = qpos - r_offset + w_offset
                    self.window_array_seq[row_index, max(0,window_ins_pos-1):min(self.window_length,window_ins_pos+1)] = 1                    
                row_index += 1
            
            self.window_positional_stats["Ins_Cnt"] = bn.nansum(self.window_insertions, axis=0)

        else:
            self.window_positional_stats["cov_aligned"] = np.zeros(self.window_length)
        
        self.window_positional_stats["read_number"] = read_num
            
        return self.window_positional_stats
          
    def call_variants(self, bq_threshold=20, low_qual_bq_thr = 10):
        if not "read_number" in self.window_description:
            self.prepare_window_metrics(bq_threshold, low_qual_bq_thr)
        read_num = self.window_description["read_number"]

        if read_num != 0:
            # loop through all positions with a SNV
            for vp in self.snv_pos_bq_filtered_good_reads.nonzero()[0]:
                # generate pileup for position
                variants = self.snv_bq_filtered_good_reads[:,vp]
                valid_read_mask = np.logical_not(np.isnan(variants))
                valid_read_snv_mask = variants[valid_read_mask]
                orig_calls = self.window_array_seq[valid_read_mask,vp]
                var_quals = self.window_array_qual[valid_read_mask,vp]
                
                # get read characteristics for reads overlapping this columns and contributing to pileup (i.e. without low qual bases and reads)
                read_characteristics = np.array(self.read_characteristics, dtype=np.dtype(object))[valid_read_mask]
                ref_base_code = self.ref_code[vp]
                ref_calls = np.sum(valid_read_snv_mask==0)
                
                # identify reads supporting any variant at this site
                variant_alleles = {1:None,2:None,3:None,4:None}
                for vv in variant_alleles.keys():
                    vv_ind = (orig_calls==vv).nonzero()[0]
                    if len(vv_ind) > 0:
                        variant_alleles[vv] = vv_ind
                
                # get reference base at site
                var_genomic_pos = vp+self.window_start
                ref_base = string.translate(str(ref_base_code), code_to_nuc_trans)

                # loop through results 
                self.variants[var_genomic_pos] = {"ref": {"base": ref_base, "allele_cnt": ref_calls, "cov_total": self.coverage_per_pos[vp]},"var":{}}
                for vv in variant_alleles.keys():
                    if variant_alleles[vv] is not None:
                        if vv == ref_base_code:
                            continue
                        supporting_reads_ind = variant_alleles[vv]
                        supporting_reads_num = len(supporting_reads_ind)
                        v_base = string.translate(str(vv), code_to_nuc_trans)
                        rc_var = [read_characteristics[i] for i in supporting_reads_ind]
                        rel_read_starts = [e["5prime_start"]-self.window_start for e in rc_var]
                        
                        callset = { "allele_cnt": supporting_reads_num,
                                    "read_plus_prop": float(bn.nansum(np.array([e["read_strand"]=="+" for e in rc_var], dtype=int)))/supporting_reads_num,
                                    "read1_plus_prop": float(bn.nansum(np.array([e["first_read_strand"]=="+" for e in rc_var], dtype=int)))/supporting_reads_num,
                                    "R1_prop": float(bn.nansum(np.array([e["mate_num"]==1 for e in rc_var], dtype=int)))/supporting_reads_num,
                                    "mq_median": bn.nanmedian(self.mq_per_read[valid_read_mask][supporting_reads_ind]),
                                    "bq_median": bn.nanmedian(var_quals[supporting_reads_ind]),
                                    "pos_median": bn.nanmedian(np.abs(np.array([rel_read_starts[i]-vp if rc_var[i]["read_strand"]=="+" else vp-rel_read_starts[i] for i in xrange(supporting_reads_num)], dtype=int)))
                                    }
                        self.variants[var_genomic_pos]["var"][v_base] = callset
    
    # FIXME: seq_startspos is not required here - remove on refactoring
    def __seq2num_vec(self, sequence, qual, seq_startpos, cigar_tuples):
        """
        Convert read alignment to numeric vector with length equal to alignment length (plus soft clipped bases)
        """
        
        cur_alignment_vec_len = self.alignment_len_default
        cig_pos = 0 # this is the position after the last inserted cigar operation relative to the leftmost alignment start on the ref seq
        q_pos = 0 # this is the corresponding position in the query sequence (i.e. the read sequence)
        seq_vec = np.zeros(cur_alignment_vec_len, dtype=int)
        qual_vec = np.zeros(cur_alignment_vec_len, dtype=int)
        ins_list = {"I": {},"D":{}} # list of indel operations in the alignment
        leftmost_match_seen = False
        start_position_offset = 0 # if the alignment starts with a soft clip then the start position of the read refers to the first really aligned base
        for cig_op in cigar_tuples:
            ll = cig_op[1]
            # check if the alignment vector is still long enough to store the remaining cigar op. If not, resize it
            if ll+cig_pos > cur_alignment_vec_len:
                seq_vec = np.append(seq_vec, np.zeros(cur_alignment_vec_len+ll, dtype=int))
                qual_vec = np.append(qual_vec, np.zeros(cur_alignment_vec_len+ll, dtype=int))                
                cur_alignment_vec_len += (cur_alignment_vec_len+ll)
            
            op_code = cig_op[0]

            if op_code in (0,7,8): # match
                seq_vec[cig_pos:(cig_pos+ll)] = [int(e) for e in string.translate(sequence[q_pos:(q_pos+ll)],self.trans_tab)]
                qual_vec[cig_pos:(cig_pos+ll)] = qual[q_pos:(q_pos+ll)]
                cig_pos += ll
                q_pos += ll
                leftmost_match_seen = True
            elif op_code == 4: # soft clipping
                seq_vec[cig_pos:(cig_pos+ll)] = -1
                qual_vec[cig_pos:(cig_pos+ll)] = qual[q_pos:(q_pos+ll)]
                cig_pos += ll
                q_pos += ll
                if not leftmost_match_seen:
                    start_position_offset += ll
            elif op_code == 2: # deletion
                seq_vec[cig_pos:(cig_pos+ll)] = 6
                ins_list["D"][q_pos] = ( ll, sequence[cig_pos:(q_pos+ll)] ) 
                cig_pos += ll
                leftmost_match_seen = True # It is not clear from the SAM spec if we could have xSyD and how that would affect start position. Lets assume that Del counts as alignment
            elif op_code == 1: # insertion
                ins_list["I"][q_pos] = ( ll, sequence[cig_pos:(q_pos+ll)] ) 
                q_pos += ll
                leftmost_match_seen = True # It is not clear from the SAM spec if we could have xSyI and how that would affect start position. Lets assume that Ins counts as alignment
            elif op_code in (5,6): # hard clipping, padding
                continue
        return (seq_vec, qual_vec, ins_list, cig_pos, start_position_offset)
         
    def __generate_window_arrays(self):
        """Build 2D arrays for seq and qual from individual reads, clipping portions that lie outside of the window
        and adjusting all reads to the same positions. Any column will then correspond to the same genomic position"""
        read_num = len(self.__read_hash.keys())
        self.window_array_seq = np.zeros((read_num, self.window_length))
        self.window_array_qual = np.zeros((read_num, self.window_length))
        self.read_characteristics = [None] * read_num
        row_index = 0
        for r in self.__read_hash.keys():
            data = self.__read_hash[r]
            rc = data[0]
            read_start = rc["read_start"] - rc["start_position_offset"] # adjust offset for soft-clipped reads
            w_offset = max(0,read_start-self.window_start)
            r_offset = max(0,self.window_start - read_start)
            r_end = min(self.window_length-w_offset+r_offset, rc["alignment_length"])
            w_end = min(self.window_length, (r_end - r_offset) + w_offset)

            try:
                self.window_array_seq[row_index, w_offset:w_end] = data[1][r_offset:r_end]
                self.window_array_qual[row_index, w_offset:w_end] = data[2][r_offset:r_end]
                self.read_characteristics[row_index] = rc
            except Exception, e:
                print "Read hash key: %s\tRead details: %s\nw_offset %d\tw_end %d\tr_offset %d\tr_end %d\tread_start %d\twindow_start %d\twindow_len %d\talignment_length %d\trow_index %d\n" % (r, self.__read_hash[r][0], w_offset, w_end, r_offset, r_end, rc["read_start"], self.window_start, self.window_length, rc["alignment_length"], row_index)
                raise e
            row_index += 1
        #print self.window_array_seq.shape
        
class RefSeqStore:
    """Fetches and keeps a larger chunk of the reference sequence in memory"""
    def __init__(self, ref_file, ref_segment_size=1e6):
        self.fasta_file = pysam.FastaFile(ref_file)
        self.known_seqs = self.fasta_file.references
        self.ref_segment_size = ref_segment_size
        self.__store = None
        self.__store_start = None
        self.__store_end = None
        self.__store_seq = None
        
    def get_ref_seq_names(self):
        return self.known_seqs
        
    def get_ref_seq_length(self, seq):
        return self.fasta_file.get_reference_length(seq)
    
    def get_segment(self, seq, start, end):
        if (not self.__store is None) and (seq == self.__store_seq and start >= self.__store_start and end <= self.__store_end):
            return self.__store[(start - self.__store_start):(end - self.__store_start+1)]
        else:
            if not seq in self.known_seqs:
                raise Exception, "Sequence [%s] not found in reference" % seq
            rl = self.fasta_file.get_reference_length(seq)
            if start < 1 or start > rl or end < start or end > rl:
                raise Exception, "Invalid start/end for reference: [%s:%d-%d]" % (seq, start, end)
            self.__store_seq = seq
            self.__store_start = start
            self.__store_end = min(start+self.ref_segment_size, rl)
            self.__store = self.fasta_file.fetch(seq, self.__store_start-1, self.__store_end )
            
        return self.__store[(start - self.__store_start):(end - self.__store_start+1)]

