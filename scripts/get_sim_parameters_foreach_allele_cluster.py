import sys
import os
import itertools
import math

def make_catalogue(allele_freq_matrix_dirpath, output_dirpath):
    """
    This script simply takes the matrices of allele frequencies and takes each entry in the matrix and turns those entries into a long catalogued list of alleles and their frequency information. Each allele has information such as position, ancestral and derived alleles, and frequency (for that time-point).
    'allele_freq_matrix_dirpath' = the directory that contains the tab delimited of the allele frequency matrix. Examples of these alleles can be found in ../input_data/allele_freq_data/ (see the README there for more info).
    'output_dirpath' - the path to the directory for which the output will be written
    """
    if allele_freq_matrix_dirpath[-1] != '/':
        allele_freq_matrix_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    for i in os.listdir(allele_freq_matrix_dirpath):
        if i[0] == '.' or i[:-4] == 'README':
            continue
        output_subdirpath = output_dirpath + i + '/'
        if not os.path.exists(output_subdirpath):
            os.makedirs(output_subdirpath)
        for j in os.listdir(allele_freq_matrix_dirpath + i):
            if j[0] == '.' or j[:-4] == 'README':
                continue
            print j
            input_filepath = allele_freq_matrix_dirpath + i + '/' + j
            filein = open(input_filepath, "r")
            output_filepath = output_subdirpath + j[:-4]
            fileout = open(output_filepath, "w")
            fileout.write('position\tancestral_allele\tderived_allele\ttotal_reads_for_population\tderived_allele_frequency\n')
            derived_alleles = filein.readline()[:-1].split('\t')[2:]
            if derived_alleles[-1] == '':
                del derived_alleles[-1]
            position = 0
            for k in filein:
                position += 1
                line = k[:-1].split('\t')
                total_reads = line[-1]
                ancestral_allele = line[1]
                if ancestral_allele == 'null':
                    pass
                else:
                    ancestral_allele = ancestral_allele.split(' ')[1]
                count = 0
                for l in line[2:-1]:
                    allele_freq = l
                    if allele_freq == 'NaN' or allele_freq == '' or allele_freq == '#VALUE!':
                        allele_freq = 'NA'
                    else:
                        allele_freq = float(allele_freq)
                    fileout.write(str(position) + '\t' + ancestral_allele + '\t' + derived_alleles[count] + '\t' + total_reads + '\t' + str(allele_freq) + '\n')
                    count += 1
            filein.close()
            fileout.close()
    return

def calc_fitness(pre_sel_filepath, post_sel_filepath, output_filepath):
    filein_pre = open(pre_sel_filepath, 'r')
    filein_pre.readline()
    filein_post = open(post_sel_filepath, 'r')
    filein_post.readline()
    fitness_catalogue = []
    print pre_sel_filepath
    print post_sel_filepath
    for i, j in itertools.izip(filein_pre, filein_post):
        line_pre = i[:-1].split('\t')
        line_post = j[:-1].split('\t')
        if [line_pre[0], line_pre[2]] != [line_post[0], line_post[2]]:
            print line_pre
            print line_post
            print 'uh oh'
            return
        #added the last three conditionals on 9.11.15
        if line_pre[4] == 'NA' or line_post[4] == 'NA' or float(line_pre[3]) == 0.0 or float(line_post[3]) == 0.0 or line_pre[1] == 'null' or line_post[1] == 'null':
            pre_freq = line_pre[4]
            post_freq = line_post[4]
            fitness = 'NA'
        else:
            pre_freq = float(line_pre[4])
            post_freq = float(line_post[4])
            pre_freq_for_fitness = pre_freq
            post_freq_for_fitness = post_freq
            if post_freq == 0.0:
                post_freq_for_fitness = 1 / float(line_pre[3])
            if pre_freq == 0.0:
                pre_freq_for_fitness = 1 / float(line_pre[3])
            fitness = math.log(post_freq_for_fitness / pre_freq_for_fitness, 10)
        #fitness is calculated by adding 1 to the numerator or denominator if
        #the value is zero, but the recorded freqs are actually zero if this 
        #is the case. May want to change this later, we don't know...
        fitness_catalogue.append([fitness, line_pre[0], line_pre[1], line_pre[2], pre_freq, post_freq])
    filein_pre.close()
    filein_post.close()
    fileout = open(output_filepath, 'w')
    fileout.write('position\tancestral_allele\tderived_allele\tpre_sel_freq\tpost_sel_freq\tfitness\n')
    fitness_catalogue = sorted(fitness_catalogue, reverse=True)
    #need to move the fitness values that are 'NA' from the front of the
    #list to the back
    index = 0
    l = fitness_catalogue[:]
    for i in l:
        if i[0] == 'NA':
            fitness_catalogue.append(i)
            index += 1
        else:
            break
    del fitness_catalogue[0:index]
    for i in fitness_catalogue:
        line = [str(j) for j in i]
        fileout.write('\t'.join(line[1:]) + '\t' + line[0] + '\n')
    fileout.close()
    return

def run_calc_fitness(allele_freq_catalogues_dirpath, output_dirpath):
    if allele_freq_catalogues_dirpath[-1] != '/':
        allele_freq_catalogues_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    for i in os.listdir(allele_freq_catalogues_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        output_genetype_dirpath = output_dirpath + i + '/'
        if not os.path.exists(output_genetype_dirpath):
            os.makedirs(output_genetype_dirpath)
        pre_sel_filepath = allele_freq_catalogues_dirpath + i + '/' + i + '_pre_selection'
        output_filepath = output_genetype_dirpath + 'replicate_1.txt'
        post_sel_filepath = allele_freq_catalogues_dirpath + i + '/' + i + '_post_selection_1'
        calc_fitness(pre_sel_filepath, post_sel_filepath, output_filepath)
        output_filepath = output_genetype_dirpath + 'replicate_2.txt'
        post_sel_filepath = allele_freq_catalogues_dirpath + i + '/' + i + '_post_selection_2'
        calc_fitness(pre_sel_filepath, post_sel_filepath, output_filepath)
    return

def get_post_sel_num_reads(allele_freq_matrix_dirpath):
    """This script walks through the 'allele_freq_matrix_dirpath' to get the number of reads collected after selection, for each of the competition experiments. It saves this information in a dictionary, where the identifying information for the allele is the index, and the number of reads is the definition. Returns this dicitonary."""
    post_sel_num_reads_dic = {}
    for i in os.listdir(allele_freq_matrix_dirpath):
        if i[0] == '.' or i[:-4] == 'README':
            continue
        for j in ['1', '2']:
            input_filepath = '%s%s/%s_post_selection_%s.txt' % (allele_freq_matrix_dirpath, i, i, j)
            filein = open(input_filepath, "r")
            filein.readline()
            for k in filein:
                line = k[:-1].split('\t')
                num_reads = line[23]
                position = line[1].split(' ')[0]
                post_sel_num_reads_dic[i + '_' + position + '_rep' + j] = num_reads
            filein.close()
    return post_sel_num_reads_dic

def get_parameters(allele_freq_matrix_dirpath, allele_freq_catalogue_dirpath, output_dirpath):
    """This script retrieves certain info for each of the alleles in the observed dataset for the purpose of getting the parameters for doing the simulations for each of the alleles. It gets the number of reads collected after selections for each of the competition experiments (i.e. subsample size) from 'allele_freq_matrix_dirpath'. And it gets the starting allele frequency prior to selection from 'allele_freq_catalogue_dirpath'. It gets this information and then writes it (one line for each allele/repetition) to 'output_filepath'. This script writes each allele's parameters as its own file (saved in 'output_dirpath'). It does this so that we can turn the simulation program into an array job, and each allele's simulation is one job on the computational cluster. The simulaitons will be done by 'create_HIV_WF_neutral_sim_database_cluster.py'."""
    if allele_freq_matrix_dirpath[-1] != '/':
        allele_freq_matrix_dirpath += '/'
    if allele_freq_catalogue_dirpath[-1] != '/':
        allele_freq_catalogue_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    #first get the number of reads for each position after
    #the selection experiment
    post_sel_num_reads_dic = get_post_sel_num_reads(allele_freq_matrix_dirpath)
    for i in os.listdir(allele_freq_catalogue_dirpath):
        if i[0] == '.' or i[:-4] == 'README':
            continue
        gene = i
        input_filepath = '%s%s/%s_pre_selection' % (allele_freq_catalogue_dirpath, i, i)
        filein = open(input_filepath, 'r')
        filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            position = line[0]
            ancestral_allele = line[1]
            derived_allele = line[2]
            start_freq = line[4]
            output_filepath_rep1 = '%s%s_%s_%s_rep1' % (output_dirpath, gene, position, derived_allele)
            fileout_rep1 = open(output_filepath_rep1, "w")
            fileout_rep1.write('gene\tposition\tancestral_allele\tderived_allele\trepatition\tstarting_freq\tending_num_of_reads\n')
            output_filepath_rep2 = '%s%s_%s_%s_rep2' % (output_dirpath, gene, position, derived_allele)
            fileout_rep2 = open(output_filepath_rep2, "w")
            fileout_rep2.write('gene\tposition\tancestral_allele\tderived_allele\trepatition\tstarting_freq\tending_num_of_reads\n')
            try:
                post_sel_num_reads = post_sel_num_reads_dic[i + '_' + position + '_rep1']
                fileout_rep1.write('%s\t%s\t%s\t%s\trep1\t%s\t%s\n' % (i, position, ancestral_allele, derived_allele, start_freq, post_sel_num_reads))
            except KeyError:
                print i + '_' + position + '_rep1'
                fileout_rep1.write('%s\t%s\t%s\t%s\trep1\t%s\tNA\n' % (i, position, ancestral_allele, derived_allele, start_freq))
            try:
                post_sel_num_reads = post_sel_num_reads_dic[i + '_' + position + '_rep2']
                fileout_rep2.write('%s\t%s\t%s\t%s\trep2\t%s\t%s\n' % (i, position, ancestral_allele, derived_allele, start_freq, post_sel_num_reads))
            except KeyError:
                print i + '_' + position + '_rep2'
                fileout_rep2.write('%s\t%s\t%s\t%s\trep2\t%s\tNA\n' % (i, position, ancestral_allele, derived_allele, start_freq))
            fileout_rep1.close()
            fileout_rep2.close()
        filein.close()
    return

def run(allele_freq_matrix_dirpath, allele_freq_catalogue_dirpath, allele_fitness_catalog_dirpath, output_dirpath):
    make_catalogue(allele_freq_matrix_dirpath, allele_freq_catalogue_dirpath)
    run_calc_fitness(allele_freq_catalogue_dirpath, allele_fitness_catalog_dirpath)
    get_parameters(allele_freq_matrix_dirpath, allele_freq_catalogue_dirpath, output_dirpath)
    return

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
