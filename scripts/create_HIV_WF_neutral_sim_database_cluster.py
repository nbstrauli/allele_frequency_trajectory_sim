#!/usr/bin/python
#$ -S /usr/bin/python
#$ -o out
#$ -e error
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00

import sys
import os
sys.path.insert(0, './')
import numpy
import math

geneticCode = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
               "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
               "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
               "GTT":"V","GTC":"V","GTA":"V","GTG":"V",

               "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
               "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
               "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
               "GCT":"A","GCC":"A","GCA":"A","GCG":"A",

               "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
               "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
               "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
               "GAT":"D","GAC":"D","GAA":"E","GAG":"E",

               "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
               "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
               "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
               "GGT":"G","GGC":"G","GGA":"G","GGG":"G", '---':'-'}

###############################
#####estimating_error_rate#####
###############################

def prime_error_rate_dic(aa_order):
    """Initializes the error_rate_dic using the correct amino acid order"""
    aa_error_rate_dic = {}
    for i in aa_order:
        #first element of definitions are the from mutation rate
        #and the second element is the to mutation rate
        aa_error_rate_dic[i] = [0.0, 0.0]
    return aa_error_rate_dic

def make_capital(codon):
    """Changes lower case codons to upper case codons"""
    dna_caps = {'c':'C', 'a':'A', 'g':'G', 't':'T'}
    codon = ''.join([dna_caps[i] for i in codon])
    return codon

def estimate_error(error_rate_data_dirpath, return_overall_error_rate=False):
    """
    This script takes in allele frequency data from sequenced wild type virus, such that each allele should theoretically be identical to the reference allele, and any deviations from this are the result of experimental error (i.e. sequencing error, pcr error, random mutations, devine intervention, etc.). See the README for proper formatting of this error rate file. The script uses this data to estimate the error rate for each possible amino acid in the data. It returns the estimated probability that a given amino acid will change to any other amino acid, and also the probability that any other amino acid will change to a given amino acid. It estimates these 'from' and 'to' rates for each amino acid.
    'error_rate_data_dirpath' - the path to the directory that contains the frequency information for each possible allele for the sequenced wildtype virus.
    """
    if error_rate_data_dirpath[-1] != '/':
        error_rate_data_dirpath += '/'

    #quickly get amino acid order from one of the input files                                              
    for i in os.listdir(error_rate_data_dirpath):
        if i[0] == '.' or i[:-4] == 'README':
            continue
        filein = open(error_rate_data_dirpath + i, "r")
        aa_order = filein.readline()[:-1].split('\t')[1:-1]
        filein.close()
        break

    #initialize the error rate dictionary using the correct
    #amino acid order
    aa_error_rate_dic = prime_error_rate_dic(aa_order)

    grand_total = 0
    error_counts = 0.0
    for i in os.listdir(error_rate_data_dirpath):
        if i[0] == '.' or i[:-4] == 'README':
            continue
        input_filepath = error_rate_data_dirpath + i
        filein = open(input_filepath, "r")
        filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            ref_codon = line[0].split(' ')[1]
            #make codon in capital letters if it is
            #not already
            if ref_codon[0] == 'a' or ref_codon[0]== 'c' or ref_codon[0]== 't' or ref_codon[0]== 'g':
                ref_codon = make_capital(ref_codon)
            #make reference allele an amino acid
            ref_aa = geneticCode[ref_codon]
            rates = [float(k) for k in line[1:-1]]
            total_reads = int(line[-1])
            row_sum = 0#sanity check (should sum to 1)
            for k in xrange(len(aa_order)):
                #use the data in each row of the matrix to update the 'from'
                #and 'to' mutation rates for each amino acid. These 'rates'
                #are actually counts right now, but they will be normalized
                #by the grand total later.
                if aa_order[k] == ref_aa:
                    from_rate = (1 - rates[k]) * total_reads
                    aa_error_rate_dic[ref_aa][0] += from_rate
                    error_counts += from_rate
                else:
                    to_rate = rates[k] * total_reads
                    aa_error_rate_dic[aa_order[k]][1] += to_rate
                row_sum += rates[k]
            if round(row_sum, 2) != 1.0:
                print 'uh oh, row should sum to one. row_sum:', row_sum#should be 1.0
            grand_total += total_reads
        filein.close()

    #now normalize all the elements of 'aa_error_rate_dic' by the
    #'grand_total' to get a rate
    for i in aa_error_rate_dic:
        aa_error_rate_dic[i][0] /= grand_total
        aa_error_rate_dic[i][1] /= grand_total
    #if the overall prob of error is desired
    #(regardless of amino acid identity) then
    #return this as well
    if return_overall_error_rate == True:
        overall_error_rate = error_counts / grand_total
        print aa_error_rate_dic
        print overall_error_rate
        return aa_error_rate_dic, overall_error_rate
    else:
        return aa_error_rate_dic

###############################
###end_estimating_error_rate###
###############################

###############################
####estimateing_pop_size#######
###############################

def get_pop_sizes(start_pop_size, max_pop_size, gens):
    """This function will calculate the pop size for all the generations up until 'gens'. The equation used to calculate the population size was found by using logistic regression to fit a logistic curve to the observed population size data. We used a logistic regression tool in the R programming language called 'glm.out' (tutorial found at: http://ww2.coastal.edu/kingw/statistics/R-tutorials/logistic.html).
    'start_pop_size' - starting population size at time 0
    'max_pop_size' - carrying capacity of the population
    'gens' - number of generations after generation 0 to estimate the population size.
    """
    pop_sizes = [start_pop_size]
    for i in xrange(gens):
        pop_size = round(1/(1+math.exp(10.74974-(2.646389*(i+1))))*max_pop_size, 0)
        pop_sizes.append(pop_size)
    return pop_sizes

###############################
###end_estimateing_pop_size####
###############################

###################################
###simulate_sequencing_and_error###
###################################

def sub_sample_add_seq_error_population(sub_samp_size, end_allele_freq, error_rate_from, error_rate_to):
    """This script adds extra sources of variance into the simulations that are the result of subsampling the true population when sequencing, and also the error in allele frequency measurement caused by sequencing error, amplifictation error, or mutation error."""
    if end_allele_freq == 0.0:
        sub_samp_allele_count = 0
        reads_lost = 0
        reads_gained = numpy.random.binomial(sub_samp_size-sub_samp_allele_count, error_rate_to)
    elif end_allele_freq == 1.0:
        sub_samp_allele_count = sub_samp_size
        reads_lost = numpy.random.binomial(sub_samp_allele_count, error_rate_from)
        reads_gained = 0
    else:
        sub_samp_allele_count = numpy.random.binomial(sub_samp_size, end_allele_freq)
        if sub_samp_allele_count == sub_samp_size:
            reads_lost = numpy.random.binomial(sub_samp_allele_count, error_rate_from)
            reads_gained = 0
        elif sub_samp_allele_count == 0:
            reads_lost = 0
            reads_gained = numpy.random.binomial(sub_samp_size-sub_samp_allele_count, error_rate_to)
        else:
            reads_lost = numpy.random.binomial(sub_samp_allele_count, error_rate_from)
            reads_gained = numpy.random.binomial(sub_samp_size-sub_samp_allele_count, error_rate_to)
    sub_samp_seq_error_allele_count = sub_samp_allele_count - reads_lost + reads_gained
    return sub_samp_seq_error_allele_count

###################################
#end_simulate_sequencing_and_error#
###################################

def run_sims(output_filepath, pop_sizes, trials, error_rate_from, error_rate_to, start_allele_freq, ending_num_reads, use_norm_approx, s):
    """
    This script actually runs the simulation. It is a relatively simple Wright-Fisher simulation. If the population is super large, then it is advised to use a normal distribution approximation to the binomial. This can cause problems if the allele freq gets close to 1 or 0, but if the pop size is REALLY big then the likelihood of an allele freq drifting to 1 or 0 becomes quite low anyway.
    'output_filepath' - The path to the output file, for which the ending allele frequency for each simulation will be written
    'pop_sizes' - a list of population sizes for each generation
    'trials' - number of simulations to run
    'error_rate_from' - rate at which the allele will mutate to any other allele
    'error_rate_to' - rate at which any other allele will mutate to the given allele
    'start_allele_freq' - Frequency at which the simulated allele starts
    'ending_num_reads' - sequencing depth (coverage) for the position at the end of the competition experiment.
    """
    fileout = open(output_filepath, "w")
    #write header for output file
    fileout.write('ending_allele_frequency_with_start_freq_of_' + str(start_allele_freq) + '\n')
    #for each simulation
    for i in xrange(trials):
        allele_freq = start_allele_freq
        if allele_freq != 0.0 and allele_freq != 1.0:
            for j in xrange(len(pop_sizes)-1): #for each gen
                if use_norm_approx == True:
                    #since pop is super big, use normal distribution to
                    #approximate the binomial
                    var = pop_sizes[j+1] * allele_freq * (1-allele_freq)
                    stdev = math.sqrt(var)
                    if stdev <=0:
                        print j, allele_freq
                        print stdev, var
                    mean = pop_sizes[j+1] * allele_freq * (1+s)
                    next_gen_allele_count = round(numpy.random.normal(mean, stdev), 0)
                elif use_norm_approx == False:
                    #use binomial if computation time is managable
                    allele_freq = allele_freq * (1+s)
                    next_gen_allele_count = float(numpy.random.binomial(pop_sizes[j+1], allele_freq))
                allele_freq = next_gen_allele_count / pop_sizes[j+1]
                if allele_freq <= 0.0:
                    allele_freq = 0.0
                    break
                elif allele_freq >=1.0:
                    allele_freq = 1.0
                    break
        #introduce noise from subsampling when sequencing
        #and also from mutation/sequencing error
        sub_samp_seq_error_allele_count = sub_sample_add_seq_error_population(ending_num_reads, allele_freq, error_rate_from, error_rate_to)
        allele_freq = sub_samp_seq_error_allele_count / ending_num_reads
        fileout.write(str(allele_freq) + '\n')
    fileout.close()
    return

def run_sims_foreach_allele(params_foreach_allele_dirpath, output_dirpath, start_popsize, max_popsize, gens, trials, seq_error_dirpath, use_norm_approx):
    """
    This module runs the Wright-Fisher simulations of the HIV population. It runs a set of simulations for each allele in the data. The simulations are meant to simulate the 'null' hypothesis, which is that each allele is neutral (i.e. has a selection coefficient equal to 0. However, in principal there is no reason why a non-zero selection coefficient could not be simulated with this framework. Most likely the user has many alleles, which require many simulations each. So, this module requires a computational cluster that uses the Sun Grid Engine (SGE) software to submit array jobs to a cluster. The user must also provide the '-t' arguement when running this script (see the README for proper usage).
    'params_foreach_allele_dirpath' - this is the path to the directory that contains the simulation parameters for each allele in the data.
    'output_dirpath' - This is the path to the directory for which the results of each set of simulations will be written
    'start_popsize' - This gives the starting population size of the simulated population
    'max_popsize' - This gives the carrying capacity (max size) of the population. This is used in the logistic function that is used to model the population size over time.
    'gens' - This gives the number of generations in the simulation
    'trials' - This gives the number of simulations that are run for each allele
    'seq_error_dirpath' - This is the path to the directory that contains the information on sequencing error.
    """

    ################ hard_coded_parameters ###################
    #this gives the selection coefficient. These are 'neutral'
    #simulations, so this is set at 0, but if one wants to 
    #simulate alleles under selection, then change this
    #accordingly
    s = 0.0
    ################ hard_coded_parameters ###################

    if seq_error_dirpath[-1] != '/':
        seq_error_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if params_foreach_allele_dirpath[-1] != '/':
        params_foreach_allele_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)

    #get error rate for each possible amino acid
    error_rate_dic = estimate_error(seq_error_dirpath)

    #now get the population size at each
    #generation
    pop_sizes = get_pop_sizes(start_popsize, max_popsize, gens)

    #choose the correct parameter file based upon
    #the value of the global variable 'SGE_TASK_ID'
    param_filepaths = []
    for i in os.listdir(params_foreach_allele_dirpath):
        if i == 'README' or i[0] == '.':
            continue
        param_filepaths.append(params_foreach_allele_dirpath + i)
    sge_task_id = int(os.environ['SGE_TASK_ID'])
    param_filepath = param_filepaths[sge_task_id-1]

    #read in parameters from the parameter
    #file
    filein = open(param_filepath, "r")
    filein.readline()#burn header
    line = filein.readline()[:-1].split('\t')
    filein.close()
    gene = line[0]
    position = line[1]
    derived_allele = line[3]
    rep = line[4]
    start_allele_freq = line[5]
    ending_num_reads = line[6]

    #check if sufficient data in param file
    if start_allele_freq == 'NA' or ending_num_reads == 'NA':
        output_filepath = '%s%s_%s_%s_%s' % (output_dirpath, gene, derived_allele, position, rep)
        fileout = open(output_filepath, "w")
        fileout.write("No simulations because start allele freq and/or ending number of reads was 'NA'.\nnull\n")
        fileout.close()
        return
    else:
        start_allele_freq = float(start_allele_freq)
        ending_num_reads = float(ending_num_reads)
        #check if the ending number of reads is zero
        if ending_num_reads == 0.0:
            output_filepath = '%s%s_%s_%s_%s' % (output_dirpath, gene, derived_allele, position, rep)
            fileout = open(output_filepath, "w")
            fileout.write("No simulations because the ending number of reads was 0 so the ending allele frequency is undefined.\nnull\n")
            fileout.close()
            return

    #this is the mutation rate that the given allele will mutate
    #to no longer be that allele
    error_rate_from = error_rate_dic[derived_allele][0]
    #this is the mutation rate where other alleles will mutate
    #to become the given allele
    error_rate_to = error_rate_dic[derived_allele][1]

    #define the output filepath based on allele identity
    output_filepath = '%s%s_%s_%s_%s' % (output_dirpath, gene, derived_allele, position, rep)

    #run simulation!
    run_sims(output_filepath, pop_sizes, trials, error_rate_from, error_rate_to, start_allele_freq, ending_num_reads, use_norm_approx, s)

    return

if __name__ == '__main__':
    if sys.argv[-1] == 'True':
        use_norm_approx = True
    elif sys.argv[-1] == 'False':
        use_norm_approx = False
    else:
        print '"use_norm_approx" variable must be either "True" or "False"'
    run_sims_foreach_allele(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7], use_norm_approx)
