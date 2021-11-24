import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from sklearn.svm import SVR
# from sklearn.svm import SVC
# from sklearn.ensemble import RandomForestRegressor 
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import confusion_matrix
# from statistics import mode
import statistics
import pickle
import sys
import random
from collections import Counter
import scipy.stats as ss
import copy

import math
import csv
from tqdm import tqdm

from collections import Counter
from itertools import repeat, chain
import itertools
from scipy.stats import gmean

from functools import cmp_to_key

# from pybloomfilter import BloomFilter


import seaborn, time
seaborn.set_style('whitegrid')

# from sklearn.linear_model import LinearRegression

# from pomegranate import BayesianNetwork

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

# def return_size(start,end,matrix,cola

class discrete_correl:
  exception_list_0=[]
  exception_list_1=[]
  exception_list_not_one=[]
  factor_0_to_1=0

def encode_discrete(matrix,df,i,j):


  a=np.array(matrix)
  temp_matrix=unique_rows(a)

  print("Sanity Check unique rows before",len(temp_matrix))

  col_name=list(df.columns)
  correl_data_struct=discrete_correl()

  print("\n")  
  print("columns",col_name[i],col_name[j])

  prim_sec_map={}
  sec_prim_map={}

  for t in range(0,len(matrix)):
    prim_sec_map[matrix[t][i]]=set([])
    sec_prim_map[matrix[t][j]]=set([])

  for t in range(0,len(matrix)):
    prim_sec_map[matrix[t][i]].add(matrix[t][j])
    sec_prim_map[matrix[t][j]].add(matrix[t][i])

  factor=0
  factor_list=[]
  temp_exception_list_0=set([])
  temp_exception_list_not_one=set([])
  
  for t in range(0,len(matrix)):
    if not(len(prim_sec_map[matrix[t][i]])==1 and len(sec_prim_map[matrix[t][j]])==1):
      temp_exception_list_not_one.add(matrix[t][i])

    if len(prim_sec_map[matrix[t][i]])==1:
      temp=0
      val=next(iter(prim_sec_map[matrix[t][i]]))
      factor=max(len(sec_prim_map[val]),factor) 
      factor_list.append(len(sec_prim_map[val]))
    else:
      temp_exception_list_0.add(matrix[t][i])  

  factor=100    
  correl_data_struct.factor_0_to_1=factor


  encoding_map={}    
  count_map={}
  count_exception=0
  max_col_1=0

  for t in range(0,len(matrix)):
    max_col_1=max(max_col_1,matrix[t][j])

  for t in range(0,len(matrix)):

    if matrix[t][i] in encoding_map:
      continue

    if matrix[t][i] in temp_exception_list_0:
      encoding_map[matrix[t][i]]=math.floor(factor*(max_col_1+4)+count_exception)
      count_exception+=1
    else:
      if matrix[t][j] in count_map:
        encoding_map[matrix[t][i]]=math.floor(matrix[t][j]*factor+count_map[matrix[t][j]])
        count_map[matrix[t][j]]+=1
      else:
        count_map[matrix[t][j]]=0
        encoding_map[matrix[t][i]]=math.floor(matrix[t][j]*factor+count_map[matrix[t][j]])
        count_map[matrix[t][j]]+=1

  one_one_val=0
  for key in count_map.keys():      
    if len(sec_prim_map[key])==1:
      one_one_val+=1

      

  for t in range(0,len(matrix)):
    matrix[t][i]=encoding_map[matrix[t][i]]

  for t in temp_exception_list_0:
    correl_data_struct.exception_list_0.append(encoding_map[t])  

  for t in temp_exception_list_not_one:
    correl_data_struct.exception_list_not_one.append(encoding_map[t])   

  print("one one mappings are:",one_one_val,"proportion:",one_one_val*1.00/len(prim_sec_map.keys()))
  
  print("many to many vals:",len(correl_data_struct.exception_list_0)*1.00/len(prim_sec_map.keys()))
  print("one to many vals:",(len(correl_data_struct.exception_list_not_one)-len(correl_data_struct.exception_list_0))*1.00/len(prim_sec_map.keys()))  
  print("one to one vals:",1.00-(len(correl_data_struct.exception_list_not_one)*1.00/len(prim_sec_map.keys())))   

  a=np.array(matrix)
  temp_matrix=unique_rows(a)

  print("Sanity Check unique rows after",len(temp_matrix))  
  print("\n\n")

  return correl_data_struct  

def analyse_fpr(matrix,df,i,j,correl_data_struct,target_fpr,block_size):

  num_blocks=math.floor(len(matrix)/block_size)

  print("num blocks:",num_blocks)


  many_many_elements=set(correl_data_struct.exception_list_0)
  one_many_elements=set(correl_data_struct.exception_list_not_one)


  size_correl=0.0
  size_normal=0.0

  block_bloom_list_0_normal=[]
  block_bloom_list_0_correl=[]
  block_bloom_list_1=[]

  block_set_0=[]
  block_set_1=[]


  for t in range(0,num_blocks):
    block_set_0.append(set([]))
    block_set_1.append(set([]))

  for t in range(0,int(block_size*num_blocks)):
    ind=math.floor(t/block_size)
    block_set_0[ind].add(matrix[t][i])
    block_set_1[ind].add(matrix[t][j])

    
  for t in range(0,num_blocks):

    count_to_add=0

    for item in block_set_0[t]:
      if item in one_many_elements:
        count_to_add+=1

    block_bloom_list_0_correl.append(BloomFilter(count_to_add, target_fpr))
    block_bloom_list_0_normal.append(BloomFilter(len(block_set_0[t]), target_fpr))
    block_bloom_list_1.append(BloomFilter(len(block_set_1[t]), target_fpr)) 

    for item in block_set_0[t]:
      block_bloom_list_0_normal[-1].add(item)
      if item in one_many_elements:
        block_bloom_list_0_correl[-1].add(item)

    # print("perecentage used:",count_to_add*1.00/len(block_set_0[t]))    

    for item in block_set_1[t]:
      block_bloom_list_1[-1].add(item)  

    size_normal+=1.44*math.log(1.00/target_fpr,2)*len(block_set_0[t])  
    size_correl+=1.44*math.log(1.00/target_fpr,2)*count_to_add

  print("Size Ratio:",size_correl*1.00/size_normal)
  # correl_bf=BloomFilter(len(correl_data_struct.exception_list_0), 0.01)
  # for item in correl_data_struct.exception_list_0:
  #   correl_bf.add(item)
  #   # print(item)

  # correl_bf_not_one=BloomFilter(len(correl_data_struct.exception_list_not_one), 0.01)
  # for item in correl_data_struct.exception_list_not_one:
  #   correl_bf_not_one.add(item)  

  # size_correl=size_normal
  # size_correl+=1.44*math.log(1.00/0.01,2)*len(correl_data_struct.exception_list_0)
  # size_correl+=1.44*math.log(1.00/0.01,2)*len(correl_data_struct.exception_list_not_one)

  num_queries_per_block=1000

  total_negatives=0
  total_false_positives_normal=0
  total_false_positives_correl=0



  


  for curr_block in tqdm(range(0,num_blocks)):
    rand_list=np.random.uniform(0,1.0,num_queries_per_block)

    for t in range(0,num_queries_per_block):
      ind=math.floor(rand_list[t]*num_blocks*block_size)

      if matrix[ind][i] in block_set_0[curr_block]:
        if matrix[ind][i] not in many_many_elements:
          val=math.floor(matrix[ind][i]/correl_data_struct.factor_0_to_1)
          if val not in block_bloom_list_1[curr_block] or val not in block_set_1[curr_block]:
            while(True):
              print("ERROR",val,matrix[ind][i],matrix[ind][j])
        continue

      total_negatives+=1


      
      if matrix[ind][i] in block_bloom_list_0_normal[curr_block]:
        total_false_positives_normal+=1
        

      if matrix[ind][i] in many_many_elements:
        if matrix[ind][i] in block_bloom_list_0_correl[curr_block]:
          total_false_positives_correl+=1
      else:
        val=math.floor(matrix[ind][i]/correl_data_struct.factor_0_to_1)
        if matrix[ind][i] in one_many_elements:
          if matrix[ind][i] in block_bloom_list_0_correl[curr_block] and val in block_bloom_list_1[curr_block]:
            total_false_positives_correl+=1
        else:
          if val in block_bloom_list_1[curr_block]:
            total_false_positives_correl+=1      
        

  fpr_correl=total_false_positives_correl*1.00/total_negatives
  fpr_normal=total_false_positives_normal*1.00/total_negatives
  print("Normal False positive rate:",fpr_normal)
  print("Correl False positive rate:",fpr_correl)       



  print("\n\n")


  return  fpr_correl,size_correl,fpr_normal,size_normal        


def chop_blocks(matrix,block_size):

    num_blocks=math.ceil(len(matrix)*1.00/block_size)

    block_list=[]
    for i in range(0,num_blocks):
        block_list.append([])

    for i in range(0,len(matrix)):
        block_id=math.ceil(i*1.00/block_size)    
        block_list[block_id-1].append(matrix[i])

    print("chopped blocks; nuumblocks:",num_blocks," avg block_size: ",block_size)    

    return block_list    


def findsubsets(s, n):
    return list(itertools.combinations(s, n))

def findallsubsets(s):
    final_list=[]
    for i in range(0,len(s)):
      temp_list=list(itertools.combinations(s, i+1))
      final_list=final_list+temp_list
    return final_list


def print_stats(block_list,col_list,col_names):

    block_id_dict={}
    count_id_dict={}
    bloom_filt_list=[]
    total_rows=0

    for i in range(0,len(block_list)):
        bloom_filt_list.append([])
        for k in range(0,len(col_list)):
            bloom_filt_list[-1].append(set([]))

        for j in range(0,len(block_list[i])):
            key=""
            for k in range(0,len(col_list)):
                bloom_filt_list[-1][k].add(float(block_list[i][j][col_list[k]]))
                key+=str(block_list[i][j][col_list[k]])+"_"

            if key not in block_id_dict.keys():
                block_id_dict[key]=set([])

            block_id_dict[key].add(i)

            if key not in count_id_dict.keys():
                count_id_dict[key]=0

            count_id_dict[key]+=1
            total_rows+=1

    total_blocks=0        
    for i in block_id_dict.keys():
        total_blocks+=len(block_id_dict[i])

    total_filt_blocks=0

    total_single_blocks=0

    for i in block_id_dict.keys():
        total_single_blocks+=math.ceil(count_id_dict[i]*1.00/len(block_list[0]))

    bloom_distinct_elements=0
    for i in range(0,len(bloom_filt_list)):
      # bloom_distinct_elements+=len(bloom_filt_list[i][-1])
      for j in range(0,len(bloom_filt_list[i])):   
        bloom_distinct_elements+=len(bloom_filt_list[i][j])

      

    print("checking filter")
    for key in tqdm(block_id_dict.keys()):
        key_arr=key.split("_")
        for j in range(0,len(bloom_filt_list)):
            bool_check=True
            for k in range(0,len(bloom_filt_list[j])):
                if float(key_arr[k]) not in bloom_filt_list[j][k]:
                    bool_check=False
                    break
            if bool_check:
                total_filt_blocks+=1        

    # optimal_val_with_qdtree=total_rows*1.00/(len(block_list[0])*len(block_id_dict.keys()))
    optimal_val_with_qdtree=total_single_blocks*1.00/(len(block_id_dict.keys()))
    optimal_with_perfect_filtering=total_blocks*1.00/len(block_id_dict.keys())
    optimal_with_single_filters=total_filt_blocks*1.00/len(block_id_dict.keys())


    print("-------------------------------------------------------")
    print("Stats for columns:",str(col_names))
    print("Distinct Queries:",len(block_id_dict.keys()),"Total Number of Blocks",len(block_list),"Block Size:",len(block_list[0]))        
    print("Optimal number of blocks per query:",optimal_val_with_qdtree)
    print("Actual Number of blocks to be touched per query:",optimal_with_perfect_filtering)
    print("Single Dim Filter blocks to be touched per query:",optimal_with_single_filters)

    print("total distinct filter elements:",bloom_distinct_elements)   
    print("Bloom Filter Size(in MB): ",bloom_distinct_elements*10.00/(8.00*pow(10,6)))
    print("Uncompressed Data Size(in MB): ",total_rows*len(block_list[0][0])*64.00/(8.00*pow(10,6)))
    print("-------------------------------------------------------")

    return optimal_with_single_filters*1.00/optimal_with_perfect_filtering, optimal_with_perfect_filtering*1.00/optimal_val_with_qdtree



def generate_data(column_card,num_rows,num_cols):

  matrix=[]

  for i in range(0,num_rows):
    temp_list=[]
    for j in range(0,num_cols):
      temp_val=random.randint(1,column_card)
      temp_val-=1
      temp_list.append(temp_val)
    matrix.append(temp_list)  

  df = pd.DataFrame.from_records(matrix)

  return df  


def generate_workload(block_list,col_list,col_names,block_index_list):

  workload=[]
  num_queries=100000
  # chop=1000

  for i in range(0,num_queries):
    
    val=1
    # while (val<=1):  
    #   val=random.randint(2,32)
    #   val=math.ceil(math.log(val,2.0))
    #   val=6-val

    # val=2
    #   # break
    # val=6

    # val_subsets=findsubsets(col_list,val)

    # for subset in range(0,len(val_subsets)):
    #   val_subsets[subset].sort()


    # for subset_id in range(0,len(val_subsets)):
    #   for 


    query_col_list=random.sample(col_list, val)
    query_col_list.sort()

    query_col_list=[1,2,3,4]
    val_list=[]

    block_index=random.randint(0,len(block_list)-1)
    block_id=random.randint(0,len(block_list[block_index])-1)

    for j in range(0,len(query_col_list)):
      val_list.append(block_list[block_index][block_id][query_col_list[j]])

    workload.append([query_col_list,val_list,1.0])
    # print(i,val,workload[-1])

  all_subsets=findallsubsets(col_list)
  for i in range(0,len(all_subsets)):
    all_subsets[i]=list(all_subsets[i])
    all_subsets[i].sort()
  # print(all_subsets)

  print("building main filter")
  all_filter_dict={}
  for i in tqdm(range(0,len(all_subsets))):
    all_filter_dict[str(all_subsets[i])]={}
    for j in range(0,len(block_index_list)):
      all_filter_dict[str(all_subsets[i])][block_index_list[j]]=set([])
      for k in range(0,len(block_list[block_index_list[j]])):
        temp_list=[]
        for p in range(0,len(all_subsets[i])):
          temp_list.append(block_list[block_index_list[j]][k][all_subsets[i][p]])
        all_filter_dict[str(all_subsets[i])][block_index_list[j]].add(str(temp_list))

  print("checking true positives")
  true_positives=0
  for i in tqdm(range(0,len(workload))):
    workload[i].append([])
    for j in range(0,len(block_index_list)):
      # print(str(workload[i][1]),str(workload[i][0]),all_filter_dict[str(workload[i][0])][block_index_list[j]])
      if str(workload[i][1]) in all_filter_dict[str(workload[i][0])][block_index_list[j]]:
        workload[i][-1].append(block_index_list[j])  
        true_positives+=1    
        # print("Added")

  # for i in range(0,len(workload)):
  #   val=len(workload[i][-1])*1.00/len(block_index_list)
  #   if(val>1.2):
  #     workload[i][2]=0.0
  #   else:
  #     workload[i][2]=1000.00/(len(workload[i][-1])+1)  

    # print("prop hit: ",len(workload[i][-1])*1.00/len(block_index_list),str(workload[i][0]),str(workload[i][1]))   
    

  print("Prop of true postives:",true_positives*1.00/(len(workload)*len(block_index_list)))

  return workload    
    

def single_dim_filter_analysis(block_list,col_list,col_names,workload,block_index_to_test):

    bloom_filt_list=[]
    total_rows=0
    print("Analysis Started")

    print("filter buiding")
    for i in tqdm(range(0,len(block_list))):
        bloom_filt_list.append([])
        for k in range(0,len(col_list)):
            bloom_filt_list[-1].append(set([]))

        for j in range(0,len(block_list[i])):
            key=""
            for k in range(0,len(col_list)):
                bloom_filt_list[-1][k].add(float(block_list[i][j][col_list[k]]))
                key+=str(block_list[i][j][col_list[k]])+"_"

            total_rows+=1
    
    total_queries=0
    total_passed=0
    true_positives=0
    fpr_passed=[0,0,0,0,0,0]
    fpr_passed_weighted=[0,0,0,0,0,0]
    print("Identifying usefullness")

    # check_tp=0
    # for i in range(0,len(workload)):
    #   check_tp+=len(workload[i][-1])

    # print("check tp: ",check_tp*1.00/(len(block_index_to_test)*len(workload)))

    total_weight=0
    true_positive_weighted=0
    for i in tqdm(range(0,len(block_index_to_test))):
      for j in range(0,len(workload)):        
        bool_present=0
        total_queries+=1
        total_weight+=workload[j][2]
        if block_index_to_test[i] in workload[j][-1]:
          true_positives+=1
          true_positive_weighted+=workload[j][2]
          # continue
        
        for k in range(0,len(workload[j][0])):
          if(workload[j][1][k] not in bloom_filt_list[block_index_to_test[i]][workload[j][0][k]]):
            bool_present+=1
        fpr_passed[bool_present]+=1
        fpr_passed_weighted[bool_present]+=workload[j][2]

    
    total_distinct_elements=0
    total_space=0
    for i in tqdm(range(0,len(block_index_to_test))):
      total_space+=5*len(block_list[block_index_to_test[i]])*64
      for k in range(0,len(bloom_filt_list[block_index_to_test[i]])):
        total_distinct_elements+=len(bloom_filt_list[block_index_to_test[i]][k])        

    fpr_list=[0.75,0.5,0.25,0.1,0.05,0.01,0.001,0.00001]

    final_prop_space=[]
    final_prop_passed=[]
    final_prop_weighted=[]

    for i in range(0,len(fpr_list)):
      prop_passed=0
      prop_weighted=0
      for j in range(0,len(fpr_passed)):
        prop_passed+=fpr_passed[j]*pow(fpr_list[i],j)
        prop_weighted+=fpr_passed_weighted[j]*pow(fpr_list[i],j)
      
      final_prop_passed.append(prop_passed*1.00/total_queries)
      final_prop_weighted.append(prop_weighted*1.00/total_weight)
      final_prop_space.append(total_distinct_elements*1.00*(math.log(1.00/fpr_list[i],2)+1)/total_space) 

    print("------------Space Usage for Single Dim Bloom Filter---------")
    print("True Positives:",true_positives*1.00/total_queries)
    for i in range(0,len(final_prop_passed)):
      print("Prop Space Used: ",final_prop_space[i]," PROP BLOCKS Read: ",final_prop_passed[i])

    print("True Positives:",true_positive_weighted*1.00/total_weight)
    for i in range(0,len(final_prop_passed)):
      print("Prop Space Used: ",final_prop_space[i]," PROP BLOCKS Read WEIGHTED: ",final_prop_weighted[i])
    # print(final_prop_passed)
    # print(final_prop_space)   
    print("-------Space Usage for Single Dim Bloom Filter--------")

    return final_prop_passed,final_prop_space


def smart_dim_filter_analysis(block_list,col_list,col_names,workload,block_index_to_test):

  all_subsets=findallsubsets(col_list)
  for i in range(0,len(all_subsets)):
    all_subsets[i]=list(all_subsets[i])
    all_subsets[i].sort()

    
  print("building main filter")
  all_filter_dict={}
  for i in tqdm(range(0,len(all_subsets))):
    all_filter_dict[str(all_subsets[i])]={}
    for j in range(0,len(block_index_to_test)):
      all_filter_dict[str(all_subsets[i])][block_index_to_test[j]]=set([])
      for k in range(0,len(block_list[block_index_to_test[j]])):
        temp_list=[]
        for p in range(0,len(all_subsets[i])):
          temp_list.append(block_list[block_index_to_test[j]][k][all_subsets[i][p]])
        all_filter_dict[str(all_subsets[i])][block_index_to_test[j]].add(str(temp_list))

  steps=25
  step_size=2000.0


  print("smart solution taking steps")
  true_weight=0
  weight_per_step=[]
  for i in range(0,steps):
    weight_per_step.append(0.0)

  for j in tqdm(range(0,len(block_index_to_test))):
    workload_weight_list=[]
  
    for k in range(0,len(workload)):
      workload_weight_list.append(workload[k][2])

    true_weight=sum(workload_weight_list)

    for step_id in range(0,steps):
      benefit_list=[]
      for subset_id in range(0,len(all_subsets)):
        
        prob=random.randint(1,10)
        temp_benefit=0
        fpr=pow(0.5,step_size/len(all_filter_dict[str(all_subsets[subset_id])][block_index_to_test[j]]))
        # print("fpr is: ",fpr)
        if prob<=10:
          for k in range(0,len(workload)):     
            if set(all_subsets[subset_id])<=set(workload[k][0]):
              # check_bool=True
              temp_list=[]
              for t in range(0,len(all_subsets[subset_id])):
                temp_index=workload[k][0].index(all_subsets[subset_id][t])
                temp_list.append(workload[k][1][temp_index])
                
              if str(temp_list) in all_filter_dict[str(all_subsets[subset_id])][block_index_to_test[j]]:
                temp_benefit+=workload_weight_list[k]
              else:
                temp_benefit+=workload_weight_list[k]*fpr
            else:
              temp_benefit+=workload_weight_list[k]

        # print(step_id,"subset:",str(all_subsets[subset_id]),temp_benefit,true_weight)    

        benefit_list.append(temp_benefit)

      subset_to_use=benefit_list.index(min(benefit_list))  

      for k in range(0,len(workload)):     
        if set(all_subsets[subset_to_use])<=set(workload[k][0]):
          temp_list=[]
          for t in range(0,len(all_subsets[subset_to_use])):
            temp_index=workload[k][0].index(all_subsets[subset_to_use][t])
            temp_list.append(workload[k][1][temp_index])
            
          if str(temp_list) not in all_filter_dict[str(all_subsets[subset_to_use])][block_index_to_test[j]]:
            workload_weight_list[k]=workload_weight_list[k]*fpr

      # print("stepid",step_id," benefit: ",sum(workload_weight_list),true_weight)
      weight_per_step[step_id]+=(sum(workload_weight_list)*1.00/(true_weight*1.00*len(block_index_to_test)))

  total_space=0.0
  for i in tqdm(range(0,len(block_index_to_test))):
    total_space+=5*len(block_list[block_index_to_test[i]])*64

  space_per_step=[]
  for i in range(0,steps):
    space_per_step.append((i+1)*step_size*len(block_index_to_test))


  print("------------Space Usage for Smart Dim Bloom Filter---------")
  # print("True Positives:",true_positives*1.00/total_queries)
  for i in range(0,len(space_per_step)):
    print("Prop Space Used: ",space_per_step[i]*1.00/total_space," PROP BLOCKS Read WEIGHTED: ",weight_per_step[i])

  # print("True Positives:",true_positive_weighted*1.00/total_weight)
  # for i in range(0,len(final_prop_passed)):
  #   print("Prop Space Used: ",final_prop_space[i]," PROP BLOCKS Read WEIGHTED: ",final_prop_weighted[i])
  # print(final_prop_passed)
  # print(final_prop_space)   
  print("-------Space Usage for Smart Dim Bloom Filter--------")  

  return 1,1

# Import the model we are using
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix

def learned_bf(block_list,col_list,col_names):

  block_index_list=[]
  for i in range(0,len(block_list)):
    block_index_list.append(i)

  sampled_block_index_list=random.sample(block_index_list,10)

  workload=generate_workload(block_list,col_list,col_names,sampled_block_index_list)

  print(workload[1:10])

  X=[]
  Y=[]
  X_test=[]
  Y_test=[]

  for i in range(0,len(workload)):
    feature_vec=[]
    for j in range(0,7):
      if j in workload[i][0]:
        temp=workload[i][0].index(j)
        feature_vec.append(workload[i][1][temp])
      else:
        feature_vec.append(0-1000)

    block_index=0
    if block_index in workload[i][-1]:
      for rep in range(0,15):
        Y.append(0)
        X.append(feature_vec)
    else:
      temp=random.randint(0,100)
      if temp<1000:
        Y.append(1)
        X.append(feature_vec)
      else:
        Y_test.append(1)
        X_test.append(feature_vec)


 

  print("false values prop: ",sum(Y)*1.00/len(Y))    
  print("true values prop: ",1.00-(sum(Y)*1.00/len(Y)))    

  # Instantiate model with 1000 decision trees
  rf = RandomForestClassifier(n_estimators = 20,max_depth=4, random_state = 42)
  # Train the model on training data
  rf.fit(X, Y)

  predictions = rf.predict(X)

  print(confusion_matrix(Y,predictions))

  conf=confusion_matrix(Y,predictions)
  KL_div=0.0
  f=conf[0][0]*1.00/(conf[0][0]+conf[0][1])
  g=conf[1][0]*1.00/(conf[1][0]+conf[1][1])

  KL_div+=(f*(math.log(f*1.00/g,2)))
  f=1-f
  g=1-g
  KL_div+=(f*(math.log(f*1.00/g,2)))
  print("KL Div is: ",KL_div)


  s=[1,2,3,4,5,6]
  subset_2=findsubsets(s,2)
  subset_2=[[1,2,3,4],[1,2,3]]

  for i in range(0,len(subset_2)):
    X_test=[]
    Y_test=[]

    for j in range(0,len(block_list[0])):
      feature_vec=[0-1000]*7
      for k in range(0,len(subset_2[i])):
        feature_vec[subset_2[i][k]]=block_list[0][j][subset_2[i][k]]
      X_test.append(feature_vec)
      Y_test.append(0)

    predictions = rf.predict(X_test)
    print("columns done: ",col_names[subset_2[i][0]],col_names[subset_2[i][1]])
    print(confusion_matrix(Y_test,predictions))


  # Calculate the absolute errors
  # errors = abs(predictions - test_labels)
  # Print out the mean absolute error (mae)
  # print('Mean Absolute Error:', round(np.mean(errors), 2), 'degrees.')







def plot_effectiveness(block_list,col_list,col_names):

  block_index_list=[]
  for i in range(0,len(block_list)):
    block_index_list.append(i)

  sampled_block_index_list=random.sample(block_index_list,40)

  workload=generate_workload(block_list,col_list,col_names,sampled_block_index_list)

  print("Workload Generated")

  smart_dim_filter_fpr,single_dim_filter_space=smart_dim_filter_analysis(block_list,[2,3,4,5,6],col_names,workload,sampled_block_index_list)

  single_dim_filter_fpr,single_dim_filter_space=single_dim_filter_analysis(block_list,[0,1,2,3,4,5,6],col_names,workload,sampled_block_index_list)

  # smart_alloc_fpr,smart_alloc_space=smart_alloc_analysis(block_list,col_names,workload)


def benchmark_vortex(file_name,acceptance_list):

  col_name=[]
  if 'dmv' in file_name:
    df=pd.read_csv(file_name)
  else:

    if "synthetic" in file_name:
      column_card=10
      num_rows=100000
      num_cols=len(acceptance_list)
      df=generate_data(column_card,num_rows,num_cols)
      df.columns=acceptance_list
    else:  
      df=pd.read_csv(file_name,header=None,sep='|')  
      for i in range(0,len(df.columns)):
        col_name.append(str(i))

    # df.columns=col_name 

  col_list=list(df.columns)  
  
  if "synthetic" not in file_name:
    for i in range(0,len(col_list)):
      # if "VIN" not in col_list[i]:
      #   continue
      # if "Date" in col_list[i]:
      #   continue
      # continue
      if col_list[i] in acceptance_list:
        continue
      df=df.drop(col_list[i], 1)  

  if "dmv" in file_name:    
    df=df.sort_values('VIN')  
    df=df.sort_values('City')
  if "synthetic" in file_name:
    df=df.sort_values('1')


  columns_list=list(df.columns)
  
  print(col_list)
  print('df stuff',df.columns)
  print("___________________")
  print("Kapil")
  print("___________________")
  matrix=df.values.tolist()
  print(columns_list)

  card_list={}

  for i in tqdm(range(0,len(matrix[0]))):

    print("col type",df.columns[i],type(matrix[0][i]),matrix[0][i])
    if type(1.00)==type(matrix[0][i]):
      df[df.columns[i]] = df[df.columns[i]].fillna(0-1)

      for j in range(0,len(matrix)):
        # print("i:",i,"j:",j,matrix[j][i])  
        if type(matrix[j][i]) == str:
            continue
        if np.isnan(matrix[j][i]):
          matrix[j][i]=0-1

    if True or type("string")==type(matrix[0][i]):
    # if True:  
      map_dict={}
      # print(matrix[i])
      temp_col=[]
      for t in range(0,len(matrix)):
        map_dict[matrix[t][i]]=0
        temp_col.append(matrix[t][i])

      print("size of col is:",len(set(temp_col)))
      print("sorting by freq")  
      # new_list = sorted(temp_col, key = temp_col.count, reverse=True)
      # new_list= a.sort(key=Counter(a).get, reverse=True)
      new_list=list(set(temp_col))
      # new_list= list(chain.from_iterable(repeat(i, c) for i,c in Counter(temp_col).most_common()))
      print("sorting by freq done")
  
      count=0  
      map_dict[new_list[0]]=0
      for t in range(0,len(new_list)):
        # print(new_list[t])
        if t==0:
          continue
        if new_list[t]!=new_list[t-1]:
          count+=1
        map_dict[new_list[t]]=count  

      print("cardinality of col",i,"is:",count) 
      
      card_list[i]=count 
    #   card_order.append(count)
      df[df.columns[i]] = df[df.columns[i]].map(map_dict)  
      for j in range(0,len(matrix)):
        # print(matrix[j][j],map_dict[matrix[i][j]])
        matrix[j][i]=map_dict[matrix[j][i]]

  # for i in tqdm(range(0,len(matrix[0]))):
  #   for j in tqdm(range(0,len(matrix[0]))):
  #     if i<=j:
  #       continue
  #     check_map_val(matrix,df,i,j)
  #     check_map_val(matrix,df,j,i)

#   random.shuffle(matrix)

  block_list=chop_blocks(matrix,100000)  

#   validate(chop_blocks,columns_list)

  learned_bf(block_list,[1,2,3,4,5,6],columns_list)

  return

  plot_effectiveness(block_list,[2,3,4,5,6],columns_list)

  return

  for subset_size in tqdm(range(0,5)):
    s=[2,3,4,5,6]
  #   subset_1=findsubsets(s,1)
    # subset_1=findsubsets(s,subset_size+1)
    subset_1=findsubsets(s,subset_size+1)
    # subset_1=[[6]]

    qd_tree_benefit=[]
    multi_dim_filiter_benefit=[]

    for i in tqdm(range(0,len(subset_1))):
        temp_col_list=[]
        for j in range(0,len(subset_1[i])):
            temp_col_list.append(columns_list[subset_1[i][j]]) 
        a,b=print_stats(block_list,subset_1[i],temp_col_list) 
        multi_dim_filiter_benefit.append(a)
        qd_tree_benefit.append(b)

    print("---------FINAL CONCLUSION-----------")
    print("Subset Size:",subset_size+1)
    # print("Subset Size:",6)
    print("QD TREE BENEFIT:",gmean(qd_tree_benefit))
    print("MULTI DIM FILTER BENEFIT:",gmean(multi_dim_filiter_benefit))
    print("---------FINAL CONCLUSION-----------")    
    # break



#   print_stats(block_list,[1,2],[columns_list[1],columns_list[2]])

  return 

  




#DMV
file_name="blockpartitioning/build/cache/dmv/small_dmv.csv"
# file_name="dmv_tiny.csv"
# acceptance_list=["City","State"]
# acceptance_list=["State","Zip"]
# acceptance_list=["Record Type","Registration Class"]
acceptance_list=["City","Zip","County","State","Color","Model Year","VIN"]

# file_name="synthetic"
# acceptance_list=["0","1","2","3","4","5","6"]
# acceptance_list=["City","Zip"]
# acceptance_list=["City","Zip","County","State","Color"]
# query_cols=[["City","Color"],["City"],["Color"],["Body Type"],["State","Color"],["City","Color","State"],["Color","Body Type"]]
# target_fpr=0.01

# benchamrk_bloom_filter_real(file_name,acceptance_list)
benchmark_vortex(file_name,acceptance_list)