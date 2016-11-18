These are the software associated with the manuscript 'Evolution drives a division of labor in overlapped genes in HIV'. The code herein functions to run neutral simulations of alleles 'drifting' in frequency over time. Broadly speaking, there are three steps to running this code. The first is to make parameter files for each of the alleles for which simulations must be run. The second is to run the simulations. And the third is to compare the observed results to the simulations.

Software in this directory:

create_HIV_WF_neutral_sim_database_cluster.py

This runs the Wright-Fisher neutral simulations of alleles. It takes a directory of parameter files (where each file in the directory gives the simulations parameters for a single allele) and simulates this allele over a given number of generations. It then writes the ending allele frequency for each of the neutral simulations to and output file. It also incorporates experimental error and subsampling (from sequencing) into the simulations. User-be-ware that this software requires a computational cluster to run. The data associated with this manuscript had a tremendous amount of alleles that required many simulations each (with large population sizes to boot), so array jobs on a computational cluster was necessary. One must have the Sun Grid Engine (SGE) software for submitting jobs to a computational cluster. If one does not have this, then the software could easily be adapted to run on a single core, or any other environment for that matter.

	USAGE:
	
	qsub -t 1-[num_parameter_files] create_HIV_WF_neutral_sim_database_cluster.py [parameter_file_directory] [output_directory] [starting_population_size] [max_population_size] [num_generations] [num_simulations] [error_directory] [use_normal_approximation]
	
	num_parameter_files - This should be an integer. It gives the number of files in the 'parameter_file_directory', which corresponds to the number of alleles to run simulations for. This tells the qsub submission how large of an array job to request. Make sure this value is right!
	parameter_file_directory - This is the path to the directory that contains the parameter files for each of the alleles. These files give the unique parameters for each of the alleles. Examples of these files can be found in ../input_data/parameter_files . These files can be made (if your input allele frequency data is formatted correctly) by the script 'get_sim_parameters_foreach_allele_cluster.py'. 
	output_directory - This is the path to a directory that has been dedicated to writing the output of the simulations. Here, the ending allele frequency for each of the simulations for a given allele will be written to a single file. So, each output file in this directory corresponds to an allele in the data.
	starting_population_size - An integer that give the size of the population at the start of each simulation (assumed to be the same for each allele).
	max_population_size - An integer that gives the maximum that the population can reach
	num_generations - An integer that gives the number of generations to run each simulation
	error_directory - This is the path to the directory that contains the allele frequency information for the sequenced wildtype virus. This is used to calculate that rate at which each amino acid changes to anything else, which is used in the simulations. See ../input_data/wildtype_mutation_rates/ and the README file therein for examples of these files. There is one file for the Rev gene and one for the Tat gene.
	use_normal_approximation - This can be either 'True' or 'False'. If true then a normal approximation is used for the binomial distribution in the Wright-Fisher simulations. Use this if there is large population size, many generations, many alleles, or just large amounts of computation necessary in general (it greatly speeds things up). Otherwise use 'False' as this uses the correct (yet slower) binomial distribution for the sims.

get_sim_parameters_foreach_allele_cluster.py

This creates and formats the parameter files that are used by 'create_HIV_WF_neutral_sim_database_cluster.py' to run the simulations. It also creates a catalog of allele frequencies for each of the alleles in the data. Mainly, it creates a directory that has one file for each allele in the data, and each of these files has unique parameter values that will be used for the simulations. See ../input_data/parameter_files/ for examples of these files. A computational cluster is not necessary for this script.

	USAGE:
	
	python get_sim_parameters_foreach_allele_cluster.py [allele_frequency_info_directory] [allele_frequency_catalog_directory] [allele_fitness_catalog_directory] [parameter_file_directory]
	
	allele_frequency_info_directory - This is a path to the the directory that contains the frequency information for each possible allele in the data. It contains files that have this frequency information both before and after the competition experiments. These files are tab delimited. See ../input_data/allele_freq_data/ for examples of these files and more info.
	allele_frequency_catalog_directory - This is the path to the directory that will contain frequency information for each of the alleles in the data. These files are made by this script. The information contained in these files is that same as that in 'allele_frequency_info_directory', but they are formatted as a long list as opposed to a matrix, which can make reading in the data easier.
	allele_fitness_catalog_directory - This is the path to a directory that will contain rough fitness estimates for each allele in the data. The script uses the frequency before and after competition to calculate rough estimates of the fitness for each of the alleles. It writes this information as a tab delimited file where each line has the rough fitness information for each allele.
	parameter_file_directory - This is the path to the output directory that will contain the simulation parameters for each of the alleles in the data.
