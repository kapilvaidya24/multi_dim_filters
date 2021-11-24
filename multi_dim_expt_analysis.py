


import json,os, random, sys
import matplotlib
matplotlib.use('agg') #this is to run this code in azure server that does not support plt.show(), we can still save the plot
import matplotlib.pyplot as plt
from os import path
import numpy as np
from scipy.stats import trim_mean
from scipy.stats.mstats import gmean
from tqdm import tqdm
import math
import copy

from matplotlib import rc, rcParams


def expt_1():

    query_subset=[]
    query_actual_block_read=[]
    query_single_dim_blocks_read=[]

    with open('subset_dmv.out') as f:
        lines = f.readlines()
        curr_query_type=""
        
        for line in lines:
            if "Stats for columns:" in line:
                curr_query_type=line[line.find('['):-1]
                query_subset.append(curr_query_type)
                continue

            if "Actual Number of blocks to be touched per query:" in line:
                line_split=line.split(" ")
                query_actual_block_read.append(float(line_split[-1]))
                continue

            if "Single Dim Filter blocks to be touched per query:" in line:
                line_split=line.split(" ")
                query_single_dim_blocks_read.append(float(line_split[-1]))
                continue    


    plt.figure(figsize=(100,6))

    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}

    matplotlib.rc('font', **font)
    color_list=['blue','green','red',"orange","black"]

    x_list=[]

    for i in range(0,len(query_single_dim_blocks_read)):
        x_list.append(i+1)

    # print("plotting: ",space_dict[rosetta_key]," ",fpr_dict[rosetta_key])
    plt.plot(x_list, query_actual_block_read,marker='o',markersize=15, lw=4, label="Actual Blocks", color=color_list[0])
    plt.plot(x_list, query_single_dim_blocks_read,marker='v',markersize=15, lw=4, label="Blocks Read", color=color_list[1])
    # plt.plot(space_dict[gcs_key], fpr_dict[gcs_key],marker='D',markersize=15, lw=4, label="SNARF", color=color_list[2])
    # plt.plot(space_dict[ef_key], fpr_dict[ef_key],marker='D',markersize=15, lw=4, label="SNARF with Elias Fanno coding", color=color_list[3])
    # plt.plot(cuckoo_space, cuckoo_fpr,marker='+',markersize=15, lw=4, label="Cuckoo Filter", color=color_list[4])

    # plt.xlabel('Space Used (bits per key)',fontsize=20,weight='bold')
    plt.legend(fontsize=16)
    #plt.xlim(0.9, 100000)
    # plt.xlabel( _metric + ' ' + subopt_type)
    plt.yscale('log')
    # plt.ylabel('FPR Acheived',fontsize=20,weight='bold')
    # plt.ylabel('95th Percentile Suboptimality')
    #plt.ylim(0.01, 100000)
    # plt.ylim(0.00001,1.5)
    plt.yticks(fontsize=40,weight='bold')
    plt.xticks(fontsize=20,weight='bold')
    plt.xticks(x_list,query_subset)
    # plt.xscale('log')
    # template_name=all_QTs[0].split('/')[-1]
    savefilename = "multi_dim_query.png"
    # savepath = path.join(base_path, 'charts/' + savefilename)
    # plt.title("Geometric mean Suboptimality on a sequence of queries")
    # plt.title("Tail Suboptimality on a sequence of queries")
    plt.tight_layout()
    plt.savefig("figures/"+savefilename)
    plt.close()
    plt.clf()







expt_1()