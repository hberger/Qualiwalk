#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:33:56 2017

@author: hilmar
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas
import bottleneck as bn
from collections import OrderedDict

pos_index_columns = ["seq","window_start","window_end"]

def spearmanr(array2d):
    """
    Spearman correlation coefficient on a matrix columns with missing values and ties
    Should give the same result as cor(x, method="spearman", use="pairwise") in R
    """
   
    ra = np.ma.masked_invalid(bn.rankdata(array2d, axis=0))
    
    ncols = ra.shape[1]
    cor_mat = np.empty((ncols, ncols), dtype=float)
    for j in xrange(ncols):
        x = ra[:,j]
        xm = ra.mask[:,j]
        
        for k in xrange(ncols):
            if j==k:
                cor_mat[j,k] = 1
                continue
            if k>j:
                continue
            y = ra[:,k]
            ym = ra.mask[:,k]
            both_valid = np.logical_not(np.logical_or(xm, ym))
            r = np.corrcoef(x[both_valid], y[both_valid])[1,0]
            cor_mat[j,k] = r
            cor_mat[k,j] = r
    
    return(cor_mat)           


def mad_fun(x, med):
    return(bn.nanmedian(np.abs(med-x)))

def mad_masked(array2d, constant=1.4826):
    """
    Compute MAD on columns of a masked array
    """
    med = bn.nanmedian(array2d, axis=0)
    
    ncol = array2d.shape[1]
    res = np.empty(ncol, dtype=float)
    for i in xrange(ncol):
        res[i] = mad_fun(array2d[:,i], med[i])
    return(res * constant)
    
def GetMetricNonconstantVariables(res):
    """
    Determine which of the columns of the result data frame should be kept for statistical analysis
    """
    all_columns = res.columns.values.tolist()
    cols_to_exclude = pos_index_columns
    accepted_cols = [e for e in all_columns if not e in cols_to_exclude]
    # exclude columns with zero variation
    ranges = np.ptp(res[accepted_cols].values, axis=0)
    sel_columns = [accepted_cols[i] for i in xrange(len(accepted_cols)) if ranges[i] != 0]
    
    return(sel_columns)

def SelectNonredundantSet(res, columns, max_cor = 0.5, mandatory_cols = ["LowQual_bases", "total_base_cnt"]):
    """
    Takes a pandas data frame where each column is one variable measured across tiles covering the whole target sequence and
    computes the spearman correlation between all variables. 
    """
    
    # convert data frame to FLOAT matrix to make sure fast BLAS methods are used
    # 
    cor_mat = spearmanr(res[columns].values.astype(float)) 
    excluded_col_ind = []
    tmp = np.copy(cor_mat)
    tmp[np.triu_indices_from(tmp, k=0)] = np.nan
    tmp = np.ma.masked_invalid(tmp)

    mandatory_column_ind = [columns.index(e) for e in mandatory_cols if e in mandatory_cols]  
    
    
    while np.ma.max(tmp) > max_cor:
        mc = np.ma.max(np.ma.masked_invalid(np.ma.abs(tmp)))
        highest_cor_ind = np.ma.nonzero(np.ma.abs(tmp) == mc)
        if (len(highest_cor_ind[0])==0): 
            break
        
        highest_cor_metric_x_ind = highest_cor_ind[0][0]
        highest_cor_metric_y_ind = highest_cor_ind[1][0]
        
        print("X: %s\tY: %s" % (columns[highest_cor_metric_x_ind], columns[highest_cor_metric_y_ind]))
        
        if (highest_cor_metric_x_ind in mandatory_column_ind) and (highest_cor_metric_y_ind in mandatory_column_ind):
            tmp[highest_cor_metric_x_ind, highest_cor_metric_y_ind] = 0
            continue
        else:
            highest_cor_metric_ind = highest_cor_metric_y_ind if highest_cor_metric_x_ind in mandatory_column_ind else highest_cor_metric_x_ind
        
        tmp[highest_cor_metric_ind, :] = np.nan
        tmp[:, highest_cor_metric_ind] = np.nan
        tmp = np.ma.masked_invalid(tmp)
        excluded_col_ind.append(highest_cor_metric_ind)
        print("Excluded: %s" % columns[highest_cor_metric_ind])
    
    sel_cols = [columns[i] for i in xrange(len(columns)) if not i in excluded_col_ind]
    return(cor_mat, sel_cols)


    
def ScaleData(res, columns):
    ma = np.ma.masked_invalid(res[columns])
    mads = mad_masked(ma)
    medians = np.ma.median(ma, axis=0)
    mads_fixed = [m if m>0 else 2.5 for m in mads]
    
    new_df = res[pos_index_columns].copy(deep=True)
    for i in xrange(len(columns)):
        cc = columns[i]
        new_df[cc] = (res[cc]-medians[i])/mads_fixed[i]
    
    return(new_df)


def compute_outlier_num(df, columns, threshold=3):
    
    outlier_matrix = np.ma.masked_invalid(df[columns].values)
    outlier_num = np.ma.sum(np.abs(outlier_matrix) > threshold, axis=1)
    
    outlier_df = df[pos_index_columns].copy(deep=True)
    outlier_df["outlier_number"] = outlier_num
    
    return(outlier_df)

def reduce_regions(df):
    """
    Collapse neighboring regions
    """
    dd = df.values
    current = dd[0]
        
    all_regions = []
    for new_reg in dd[1:len(dd)]:
        
        if new_reg[0] == current[0]:
            if current[2]+1 >= new_reg[1]:
                current[2] = new_reg[2]
                continue
    
        new_entry = current
        all_regions.append(new_entry)
        current = new_reg
    
    new_entry = current
    all_regions.append(new_entry)
    res = pandas.DataFrame(all_regions)
    res.columns = ['#seq','start','end','outlier_cnt']
    return(res)

def DefineAcceptedRegions(df, cnt_column, max_outliers = 1):
    """
    """
    outlier_rows = df.query("%s <= %d" % (cnt_column, max_outliers) )
    
    return((reduce_regions(outlier_rows),None)) 
 
###############################################################################    
if __name__ == '__main__':
        
    import argparse, os
    from matplotlib import pyplot as plt
    import feather
#    import time
    
    parser = argparse.ArgumentParser(description='Generate accepted regions from QualiWalker output.')
    parser.add_argument('input_file', type=str, help='Input BAM file')

    parser.add_argument('-o','--output-file', action='store', type=str, default=None,
                       help='Output file prefix')

    parser.add_argument('-c','--max-correlation', action='store', type=float, default=0.5,
                       help='Maximum correlation up to which parameters will be removed from non-redundant set.')

    args = parser.parse_args()

    inp_file = args.input_file
    
    if args.output_file is None:
        output_file_prefix = os.path.splitext(inp_file)[0]
    else:
        output_file_prefix = args.output_file
    
    header_lines = pandas.read_csv(inp_file, sep="\t", nrows=10)
    index_cols = pos_index_columns
    value_cols = [e for e in header_lines.columns.values if not e in index_cols]
    dtypes={'seq': str, 'window_start': np.int64, 'window_end': np.int64}
    for v in value_cols:
        dtypes[v] = float
        
    all_window_descriptions = pandas.read_csv(inp_file, sep="\t", dtype=dtypes)
    #all_window_descriptions['seq'].astype(str)
    
    # set index to seq, start, end
    all_window_descriptions.set_index(index_cols, inplace=True, drop=False)
    
    # determine windows without coverage and set all parameter values to missing
    uncovered_windows = (all_window_descriptions["read_number"]==0).astype(bool)
    value_col_flags = np.array([True if e in value_cols else False for e in all_window_descriptions.columns.values], dtype=bool)
    all_window_descriptions.loc[uncovered_windows,value_cols] = np.nan
    
    # Determine which columns are numeric and not constant
    numeric_var_cols = GetMetricNonconstantVariables(all_window_descriptions)
    cols_to_be_excluded = ["soft_clipped_bases"]
    
    # Scale numeric data by Median and MAD
    scaled_data = ScaleData(all_window_descriptions, numeric_var_cols)
 
    # write scaled data
    output_file_scaled_data = output_file_prefix+"_scaled.feather"
    # Solutions below produce either larger files or a slower for exporting data
    #scaled_data.to_csv(output_file_scaled_data +".txt", sep="\t", index=False)
    #scaled_data.to_hdf(output_file_scaled_data, 'scaled_table', complevel=5, complib="zlib")
    #store_export = pandas.HDFStore(output_file_scaled_data, complevel = 5, complib="zlib")
    #store_export.append('scaled_metrics', scaled_data, data_columns=scaled_data.columns)
    #store_export.close()
    feather.write_dataframe(scaled_data, output_file_scaled_data)
    
    metrics_cor, sel_cols = SelectNonredundantSet(all_window_descriptions, [e for e in numeric_var_cols if not e in cols_to_be_excluded], max_cor=0.5)
#    plt.imshow(metrics_cor, interpolation='none')
#    plt.show()
    print sel_cols
    
    outl_df = compute_outlier_num(scaled_data, sel_cols)
    outl_df.loc[np.isnan(scaled_data['read_number'].values),"outlier_number"] = len(sel_cols)+1

    output_file_outliers = output_file_prefix+"_outlier_num_all.feather"
    # write outlier numbers per window
    #outl_df.to_csv(output_file_outliers+".txt", sep="\t", index=False)
    #outl_df.to_hdf(output_file_outliers, 'outlier_table', complevel=5, complib="zlib")
    feather.write_dataframe(outl_df, output_file_outliers)

    # define accepted regions    
    output_file_accepted_regions =  output_file_prefix+"_accepted_regions.bed"
    accepted_regions, rejected_regions = DefineAcceptedRegions(outl_df, "outlier_number")
    accepted_regions.to_csv(output_file_accepted_regions, sep="\t", index=False)
    