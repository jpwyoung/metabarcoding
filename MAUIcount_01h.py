#!/usr/bin/env python
import os
import subprocess
import copy
from Levenshtein import hamming
# Requires python-Levenshtein module. You can install it with pip install python-Levenshtein.


"""
This script reads amplicon sequences from a set of fastq files.
The first seqid_len bases of each read are a random tag (seqid or Unique Molecular Identifier).
The script keeps track of how many times each seqid is used with each unique sequence.
For the set of samples, outputs are files with a list of fasta sequences in descending rank order of abundance.
and corresponding tables with the counts of each sequence in each sample.
Three sets are produced using UMIs:
1. the most abundant sequences, unfiltered
2. filtered to remove sequences that occur too often, over all samples, as the second sequence
when a seqid is found with more than one sequence - these are usually chimeras or PCR mutants.
3. filtered on a per-sample basis, rather than on totals across all samples. If allele frequencies
vary greatly across samples, this would be preferable in principle, but can lead to sequences
being stochastically deleted from some samples but not others unless read counts are very high.

In addition, two output sets are produced using 'conventional' analysis without UMIs:
4. the most abundant sequences
5. clustering together sequences that differ by just one nucleotide from a more abundant one.

In all cases, the outputs are truncated to discard very rare sequences that would have
frequencies below add_limit in the overall set of samples. 

Output 2 ("accepted_sequences") is the optimal output of the MAUI-seq method; 
the other outputs are provided for comparison, validation and problem-solving.

A file MC_parameters.py specifies the files to be analysed and gene-specific parameters.
It should be in the same folder as this script (or another path that will be found).

If there is a file fastq_file_list.txt in the same folder as the data, only the files listed 
in this file will be included in the analysis. If this file is not present, it will be  
created with a list of all files that have the extension .fasta.

Written by Peter Young. Version 01h on 5 June 2019.
"""

#Parameters that control the stringency of the analysis
#These default values will be replaced by any listed in MC_parameters.py
read_diff = 2 
#Count seqid only if most abundant sequence has at least read_diff more reads than the next
reject_threshold = 1.0 
#Reject sequences that occur as second sequences with seqids at least reject_threshold
#times as often as they occur as the primary sequence
add_limit = 0.001
#sequences are included in rank order until the next would add a fraction less than add_limit
#i.e. this discards sequences with an overall relative abundance less than add_limit

#Get info on file locations and gene-specific parameters from MC_parameters.py
from MC_parameters import *


def find_match(line,dic):
    """
    Split a sequence line into seqid and sequence (removing the primers);
    then make a dictionary of all the different seqids with counts of each of their sequences
    Calls increment()
    """
    seqid = line[0:seqid_len]
    sequence = line[(seqid_len + f_primer_len):(len(line) - r_primer_len)]
    if seqid in dic:
        increment(dic[seqid],sequence,1)
    else:
        dic[seqid] = {sequence:1}

def increment(dic,key,count):
    """
    Increments the total for key in dic by count, or creates key:count if not existing
    """
    if key in dic:
        dic[key] += count
    else:
        dic[key] = count
    

sample_tuples = []
all_samples_table = {}
total_counts = {}
total_clean_by_sample_counts = {} #counts only those that pass per-sample criterion
clean_all_samples_table = {}
total_sec_seq_counts = {} 
total_cleaned_sec_counts = {}

all_samples_reads = {}
total_reads = {}
raw_read_count = 0
                    
#Get a list of the samples to process
#The sample ID will be the fasta filename up to the first "."

if not os.path.isfile(working_folder+"fastq_file_list.txt"):
    subprocess.call("ls " + working_folder + "*.fastq > " + working_folder +"fastq_file_list.txt", shell=True)

with open(working_folder + "fastq_file_list.txt") as file_list:
    for fastq_filename in file_list:

        if "/" in fastq_filename:
            fastq_filename = fastq_filename[fastq_filename.rindex("/")+1:]
        
        fastq_filepath = working_folder + fastq_filename.rstrip("\n")
        sample_ID = fastq_filename[0:fastq_filename.find(".")]
        sample_tuples.append((sample_ID,fastq_filepath))
            
for (sample_ID,fastq_filepath) in sample_tuples:            
            
    seqid_dict = {}
    read_dict = {} #keys will be sequences, items will be number of reads for each sequence in the sample

    #Read in a fastq file one sequence record at a time (4 lines) and process the DNA sequence (2nd line) with find.match
    #Create seqid_dict with structure {seqid:{sequence:count,...},...}
    with open(fastq_filepath) as fastq_file:
        ctr = 0
        record = []
        for next_line in fastq_file:
            record.append(next_line.rstrip("\n"))
            ctr += 1
            if ctr == 4:
                raw_read_count +=1
                if total_len -2 <= len(record[1]) <= total_len + 2: #Only process sequences that are expected length +/- 2 bases 
                    find_match(record[1],seqid_dict)
                    sequence = record[1][(seqid_len + f_primer_len):(len(record[1]) - r_primer_len)]
                    increment(read_dict, sequence, 1) 
                    increment(total_reads, sequence, 1) 
                record = []
                ctr = 0

    #Extract sequence count data from seqid_dict
    sample_counts = {}  #keys will be sequences, items will be number of seqids for each sequence in the sample
    sec_seq_counts = {} #keys will be sequences, items will be how often they are second sequence in a seqid           
    clean_sample_counts = {} #only includes sequences that are below the threshold for second sequence count
        
    for seqid, matches in sorted(seqid_dict.items(), key=lambda item: sum(item[1].values()), reverse=True):
   
        sorted_list = sorted(seqid_dict[seqid].items(), key=lambda item: item[1], reverse=True)
        sequence = sorted_list[0] #Choose the most abundant sequence for each seqid
        if len(sorted_list) > 1: #There is a second sequence with this seqid
            sec_seq = sorted_list[1]
        else:
            sec_seq = ("null",0)
            
        if sequence[1] - sec_seq[1] >=read_diff: #Only include sequences that have at least read_diff more reads than the sec_seq                  
            increment(sample_counts, sequence[0], 1)
            increment(total_counts, sequence[0], 1)
             
            #Count second sequences across all seqids
            if sec_seq[0] != "null":
                increment(sec_seq_counts, sec_seq[0], 1)
                increment(total_sec_seq_counts, sec_seq[0], 1)
                   
    for seq in sample_counts:
        seq_count = sample_counts[seq]
        if seq in sec_seq_counts: 
            sec_count = sec_seq_counts[seq]
        else: sec_count = 0
        
        if sec_count < seq_count * reject_threshold: 
            clean_sample_counts[seq] = seq_count
            increment(total_clean_by_sample_counts, seq, seq_count)
            increment(total_cleaned_sec_counts, seq, sec_count)

    all_samples_table[sample_ID] = sample_counts
    clean_all_samples_table[sample_ID] = clean_sample_counts    
    
    all_samples_reads[sample_ID] = read_dict #Record read data for this sample
    
#The next part goes through the total_reads in order of decreasing number of reads and clusters lower-ranking reads if they 
# differ by just one position from the focal sequence or are a subsequence of it.
# Read counts are amalgamated in all_samples_clusters for each sample and in cluster_dict overall.
# Clustering stops when the number of reads in a new cluster drops below add_limit.

remaining_seqs = copy.deepcopy(total_reads)
cluster_dict = {}
all_samples_clusters = {}
for (sample_ID,fastq_filepath) in sample_tuples:
    all_samples_clusters[sample_ID] = {}
amalgamated_seq_counts = {}
cumul_count = 0
end_flag = 0 #end_flag terminates the clustering if the last attempted cluster had fewer than cutoff_count sequences

while (len(remaining_seqs) > 0) and (end_flag == 0):
    max_seq = max(remaining_seqs, key=remaining_seqs.get)
    current_cluster = remaining_seqs[max_seq]
    remaining_seqs.pop(max_seq)
    seqs_in_cluster = 1
    for sample_ID in all_samples_reads:
        if max_seq in all_samples_reads[sample_ID]:
            all_samples_clusters[sample_ID][max_seq] = all_samples_reads[sample_ID][max_seq]

    for next_seq, matches in sorted(remaining_seqs.items(), key=lambda x:x[1], reverse=True):

        if len(next_seq) != len(max_seq):
            mismatches = 99
            #wrong length
        else:
            mismatches = hamming(max_seq, next_seq)

        if mismatches <= 1 or (next_seq in max_seq): #allow missing bases at ends
        
            #add_to_cluster
            current_cluster += remaining_seqs[next_seq]
            
            for sample_ID in all_samples_reads:
                if next_seq in all_samples_reads[sample_ID]:
                    increment(all_samples_clusters[sample_ID], max_seq, all_samples_reads[sample_ID][next_seq])

            remaining_seqs.pop(next_seq)
        
            seqs_in_cluster += 1

    if current_cluster >= cumul_count*add_limit:
        cluster_dict[max_seq] = current_cluster
        cumul_count += current_cluster
        amalgamated_seq_counts[max_seq] = seqs_in_cluster
    else:
        end_flag = 1

unclustered_seqs = len(remaining_seqs) + seqs_in_cluster
unclustered_reads = sum(remaining_seqs.values()) + current_cluster
  
#Write sets of accepted sequences in fasta format
#Headers have raw sequence rank, total counts, total secondary counts  

os.mkdir(working_folder + "MAUIcount_output")
output_folder = working_folder + "MAUIcount_output/"                    

fasfilename0 = output_folder + "all_primary_sequences.fas"
fasfile0 = open(fasfilename0, "w")
fasfilename1 = output_folder + "accepted_sequences.fas"
fasfile1 = open(fasfilename1, "w")
fasfilename2 = output_folder + "accepted_by_sample_sequences.fas"
fasfile2 = open(fasfilename2, "w")
fasfilename3 = output_folder + "read_sequences.fas"
fasfile3 = open(fasfilename3, "w")
fasfilename4 = output_folder + "cluster_sequences.fas"
fasfile4 = open(fasfilename4, "w")


rank = 0
ranked_sequence_list = []
sequence_list = [[],[],[],[],[]]
cumultotal = [0,0,0,0,0]
total_accepted_seqs = 0
total_accepted_counts = 0
reported_seqid_seqs = 0
reported_accepted_seqs = 0
reported_read_seqs = 0

for sequence, seq_count in sorted(total_counts.items(), key=lambda item: item[1], reverse=True):
    
    if sequence in total_sec_seq_counts:
        sec_seq_count = total_sec_seq_counts[sequence]
    else: sec_seq_count = 0
  
    rank += 1
    ranked_sequence_list.append(sequence)
    
    if seq_count > cumultotal[0]*add_limit:
        cumultotal[0] += seq_count
        sequence_list[0].append(sequence)
        reported_seqid_seqs +=1
    
        fasfile0.write(">seq_%d_%d_%d" % (rank, seq_count, sec_seq_count))
        fasfile0.write("\n")
        fasfile0.write(sequence)
        fasfile0.write("\n")
    
    if sec_seq_count < seq_count*reject_threshold:
        total_accepted_seqs += 1
        total_accepted_counts += seq_count 
        if seq_count > cumultotal[1]*add_limit:
            cumultotal[1] += seq_count
            sequence_list[1].append(sequence)
            reported_accepted_seqs +=1
        
            fasfile1.write(">seq_%d_%d_%d" % (rank, seq_count, sec_seq_count))
            fasfile1.write("\n")
            fasfile1.write(sequence)
            fasfile1.write("\n")

for sequence, seq_count in sorted(total_clean_by_sample_counts.items(), key=lambda item: item[1], reverse=True):
    if seq_count > cumultotal[2]*add_limit:   #omit sequences with few counts
        cumultotal[2] += seq_count
        
        if sequence in total_cleaned_sec_counts:
            sec_seq_count = total_cleaned_sec_counts[sequence]
        else: sec_seq_count = 0
      
        rank = ranked_sequence_list.index(sequence) + 1
        sequence_list[2].append(sequence)
        
        fasfile2.write(">seq_%d_%d_%d" % (rank, seq_count, sec_seq_count))
        fasfile2.write("\n")
        fasfile2.write(sequence)
        fasfile2.write("\n")

rank_r = 0
ranked_sequence_list_r = []

for sequence, seq_count in sorted(total_reads.items(), key=lambda item: item[1], reverse=True):

    rank_r +=1
    ranked_sequence_list_r.append(sequence)
    
    if seq_count > cumultotal[3]*add_limit:
        cumultotal[3] += seq_count
        sequence_list[3].append(sequence)
        reported_read_seqs +=1
        
        #Headers have raw sequence rank, total counts                          
        fasfile3.write(">seqr_%d_%d" % (rank_r, seq_count))
        fasfile3.write("\n")
        fasfile3.write(sequence)
        fasfile3.write("\n")
        
    if sequence in cluster_dict:
        #add_limit has already been applied to clusters, so don't need cumultotal
        sequence_list[4].append(sequence)

#Headers have raw sequence rank, number of sequences amalgamated, total counts                          
        fasfile4.write(">seqr_%d_%d_%d" % (rank_r, cluster_dict[sequence], amalgamated_seq_counts[sequence]))
        fasfile4.write("\n")
        fasfile4.write(sequence)
        fasfile4.write("\n")
    
    
                
fasfile0.close()
fasfile1.close() 
fasfile2.close()  
fasfile3.close() 
fasfile4.close()  
 

#Write tables of counts for each sequence in each sample for seqid-based analyses
        
tablefilename0 = output_folder + "all_primary_sequences.tab"
tablefile0 = open(tablefilename0, "w")
tablefilename1 = output_folder + "accepted_sequences.tab"
tablefile1 = open(tablefilename1, "w")
tablefilename2 = output_folder + "accepted_by_sample_sequences.tab"
tablefile2 = open(tablefilename2, "w")
tablefilename3 = output_folder + "read_sequences.tab"
tablefile3 = open(tablefilename3, "w")
tablefilename4 = output_folder + "cluster_sequences.tab"
tablefile4 = open(tablefilename4, "w")


#Write sequence ranks as column headers

for sequence in sequence_list[0]:
    rank = ranked_sequence_list.index(sequence) +1
    tablefile0.write("\tseq_%d" % (rank))
    
for sequence in sequence_list[1]:
    rank = ranked_sequence_list.index(sequence) +1
    tablefile1.write("\tseq_%d" % (rank))
    
for sequence in sequence_list[2]:
    rank = ranked_sequence_list.index(sequence) +1
    tablefile2.write("\tseq_%d" % (rank))

tablefile0.write("\n")
tablefile1.write("\n")
tablefile2.write("\n")

#Write a row of counts for each sample

for sample_ID in sorted(all_samples_table):

    tablefile0.write(sample_ID+"\t")
    tablefile1.write(sample_ID+"\t")
        
    for sequence in ranked_sequence_list:
        if sequence in all_samples_table[sample_ID]:
            seq_count = all_samples_table[sample_ID][sequence]
        else:
            seq_count = 0
        if sequence in sequence_list[0]:    
            tablefile0.write(str(seq_count)+"\t")
        if sequence in sequence_list[1]:
            tablefile1.write(str(seq_count)+"\t")
                    
    tablefile0.write("\n")
    tablefile1.write("\n")

for sample_ID in sorted(clean_all_samples_table):

    tablefile2.write(sample_ID+"\t")
        
    for sequence in sequence_list[2]:
        if sequence in clean_all_samples_table[sample_ID]:
            seq_count = clean_all_samples_table[sample_ID][sequence]
        else:
            seq_count = 0
        tablefile2.write(str(seq_count)+"\t")
                    
    tablefile2.write("\n")

#Write total counts for each sequence under corresponding column

tablefile0.write("total")
tablefile1.write("total")
tablefile2.write("total")

for sequence in sequence_list[0]:
    tablefile0.write("\t"+str(total_counts[sequence]))
    
for sequence in sequence_list[1]:
    tablefile1.write("\t"+str(total_counts[sequence]))
    
for sequence in sequence_list[2]:
    tablefile2.write("\t"+str(total_clean_by_sample_counts[sequence]))

tablefile0.write("\nseconds")
tablefile1.write("\nseconds")
tablefile2.write("\nseconds")

for sequence in sequence_list[0]:
    if sequence in total_sec_seq_counts:
        sec_seq_count = total_sec_seq_counts[sequence]
    else: sec_seq_count = 0    
    tablefile0.write("\t"+str(sec_seq_count))

for sequence in sequence_list[1]:
    if sequence in total_sec_seq_counts:
        sec_seq_count = total_sec_seq_counts[sequence]
    else: sec_seq_count = 0    
    tablefile1.write("\t"+str(sec_seq_count))
    
for sequence in sequence_list[2]:
    if sequence in total_cleaned_sec_counts:
        sec_seq_count = total_cleaned_sec_counts[sequence]
    else: sec_seq_count = 0    
    tablefile2.write("\t"+str(sec_seq_count))

tablefile0.close()
tablefile1.close()
tablefile2.close()

#Write a similar table based on reads not seqids

#Write sequence ranks as column headers
for sequence in sequence_list[3]:
    rank_r = ranked_sequence_list_r.index(sequence) +1
    tablefile3.write("\tseqr_%d" % (rank_r))  
tablefile3.write("\n")

#Write a row of counts for each sample
for sample_ID in sorted(all_samples_reads):
    tablefile3.write(sample_ID+"\t")
    for sequence in sequence_list[3]:
        if sequence in all_samples_reads[sample_ID]:
            seq_count = all_samples_reads[sample_ID][sequence]
        else:
            seq_count = 0
        tablefile3.write(str(seq_count)+"\t")                    
    tablefile3.write("\n")

#Write total counts for each sequence under corresponding column
tablefile3.write("total")
for sequence in sequence_list[3]:
    tablefile3.write("\t"+str(total_reads[sequence]))

tablefile3.close()

#Write a similar table based on clustered sequences (within 1 nt)

#Write sequence ranks as column headers
for sequence in sequence_list[4]:
    rank_r = ranked_sequence_list_r.index(sequence) +1
    tablefile4.write("\tseqr_%d" % (rank_r))  
tablefile4.write("\n")

#Write a row of counts for each sample
for sample_ID in sorted(all_samples_reads):
    tablefile4.write(sample_ID+"\t")
    for sequence in sequence_list[4]:
        if sequence in all_samples_clusters[sample_ID]:
            seq_count = all_samples_clusters[sample_ID][sequence]
        else:
            seq_count = 0
        tablefile4.write(str(seq_count)+"\t")                    
    tablefile4.write("\n")

#Write total counts for each sequence under corresponding column
tablefile4.write("total")
for sequence in sequence_list[4]:
    tablefile4.write("\t"+str(cluster_dict[sequence]))
    
tablefile4.close()

#Write a text file with a summary of parameters and various counts
with open(output_folder + "summary.txt", "w") as summary:
    summary.write(str(read_diff) + "\tread_diff\n")
    summary.write(str(reject_threshold) + "\treject_threshold\n")
    summary.write(str(add_limit) + "\tadd_limit\n\n")
    summary.write(str(raw_read_count) + "\tRaw read count\n")
    summary.write(str(sum(total_reads.values())) + "\tTotal reads used\n")
    summary.write(str(len(total_reads)) + "\tTotal unique sequences\n")
    summary.write("\nBefore truncation:\n")
    summary.write(str(sum(total_counts.values())) + "\tTotal UMI counts\n")
    summary.write(str(total_accepted_counts) + "\tTotal accepted counts\n")
    summary.write(str(len(total_counts)) + "\tTotal UMI sequences\n")
    summary.write(str(total_accepted_seqs) + "\tTotal accepted sequences\n")
    summary.write("\nAfter truncating to threshold:\n")
    summary.write(str(reported_seqid_seqs) + "\tReported UMI sequences\n")
    summary.write(str(reported_accepted_seqs) + "\tReported accepted sequences\n\n")
    summary.write(str(reported_read_seqs) + "\tReported read sequences\n")
    summary.write(str(len(cluster_dict)) + "\tReported number of clusters \n")
    summary.write(str(sum(cluster_dict.values())) + "\tReported number of reads in clusters \n")
    summary.write(str(unclustered_seqs) + "\tSequences remaining unclustered\n")
    summary.write(str(unclustered_reads) + "\tReads remaining unclustered\n")
    
