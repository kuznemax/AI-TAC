import numpy as np
import json
import os
import sys
import _pickle as pickle
import logging
import preprocess_utils

# Configure the logging format
logging.basicConfig(
    level=logging.INFO,  # Set the desired logging level
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),  # Output log messages to console
        logging.FileHandler('log.txt')  # Output log messages to file
    ]
)

#parameters
num_chr = 20

#input files:
data_file = sys.argv[1] #bed file with peak locations
intensity_file = sys.argv[2] #text file with normalized peak heights
reference_genome = sys.argv[3] #directory with reference genome fasta files 
output_directory = sys.argv[4] #output directory

directory = os.path.dirname(output_directory)
if not os.path.exists(directory):
    os.makedirs(directory)



# read bed file with peak positions, and keep only entries with valid activity vectors
logging.info("read bed started")
positions = preprocess_utils.read_bed(data_file)
logging.info("read bed completed")

# read reference genome fasta file into dictionary
logging.info("read genome started")
if not os.path.exists('../data/chr_dict.pickle'):
    chr_dict = preprocess_utils.read_fasta(reference_genome, num_chr)
    pickle.dump(chr_dict, open('../data/chr_dict.pickle', "wb"))

chr_dict = pickle.load(open('../data/chr_dict.pickle', "rb"))
logging.info("read genome completed")

logging.info("get sequences started")
one_hot_seqs, peak_seqs, invalid_ids, peak_names = preprocess_utils.get_sequences(positions, chr_dict, num_chr)
logging.info("get sequences completed")
logging.info("one_hot_seqs size: " + one_hot_seqs.size)
logging.info("peak_seqs size: " + peak_seqs.size)
logging.info("invalid_ids size: " + invalid_ids.size)
logging.info("peak_names size: " + peak_names.size)

# remove invalid ids from intensities file so sequence/intensity files match
logging.info("format_intensities started")
cell_type_array, peak_names2 = preprocess_utils.format_intensities(intensity_file, invalid_ids)
logging.info("format_intensities completed")

cell_type_array = cell_type_array.astype(np.float32)


#take one_hot_sequences of only peaks that have associated intensity values in cell_type_array
peak_ids = np.intersect1d(peak_names, peak_names2)

idx = np.isin(peak_names, peak_ids)
peak_names = peak_names[idx]
one_hot_seqs = one_hot_seqs[idx, :, :]
peak_seqs = peak_seqs[idx]


idx2 = np.isin(peak_names2, peak_ids)
peak_names2 = peak_names2[idx2]
cell_type_array = cell_type_array[idx2, :]

if np.sum(peak_names != peak_names2) > 0:
    print("Order of peaks not matching for sequences/intensities!")



# write to file
np.save(output_directory + 'one_hot_seqs.npy', one_hot_seqs)
np.save(output_directory + 'peak_names.npy', peak_names)
np.save(output_directory + 'peak_seqs.npy', peak_seqs)

with open(output_directory + 'invalid_ids.txt', 'w') as f:
    f.write(json.dumps(invalid_ids))
f.close()

#write fasta file
with open(output_directory + 'sequences.fasta', 'w') as f:
    for i in range(peak_seqs.shape[0]):
        f.write('>' + peak_names[i] + '\n')
        f.write (peak_seqs[i] + '\n')
f.close()


np.save(output_directory + 'cell_type_array.npy', cell_type_array)