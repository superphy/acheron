import numpy as np
import pandas as pd
import gc
from concurrent.futures import ProcessPoolExecutor
import os, sys
from itertools import repeat
from sklearn.feature_selection import f_classif
from acheron.workflows import supervised_model

def make_row(filename, dataset, kmer_length):
    from Bio import Seq, SeqIO
    """
    Given a genome file, create and return a row of kmer counts
    to be inserted into the kover_11mer_matrix.mer matrix.
    """
    relevant_feats = np.load("data/{}/wgs/master_{}mers/all_kmers.npy".format(
    dataset, kmer_length))

    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    # Create a temp row to fill and return (later placed in the kmer_matrix)
    temp_row = [0]*len(relevant_feats)

    jf_path = "data/{}/wgs/{}mer_jellyfish_results/{}".format(
    dataset, kmer_length, filename)

    for record in SeqIO.parse(jf_path[:-3], "fasta"):
        # Retrieve the sequence as a string
        kmer_seq = record.seq
        kmer_seq = str(kmer_seq)
        #kmer_seq = kmer_seq._get_seq_str_and_check_alphabet(kmer_seq)

        kmer_count = int(record.id)
        if kmer_count>255:
            kmer_count = 255
        temp_row[cols_dict[kmer_seq]] = kmer_count

    return filename, temp_row


def make_matrix(dataset, kmer_length, matrix_dtype, num_threads, split_num):

    relevant_feats = np.load("data/{}/wgs/master_{}mers/all_kmers.npy".format(
    dataset, kmer_length))

    splits_array = np.load("data/{}/wgs/master_{}mers/matrix_splits.npy".format(
    dataset, kmer_length), allow_pickle=True)

    genomes = splits_array[int(split_num)-1]
    genomes = [i.split('/')[-1] for i in genomes if i != None]

    total = len(genomes)

    runs = [i.split('.')[0] for i in genomes]

    # declaring empty kmer matrix to fill
    kmer_matrix = np.zeros((len(genomes),len(relevant_feats)),dtype = matrix_dtype)

    # making dicts for faster indexing
    # note that rows dict is in filenames not genome/run names
    rows_dict = { genomes[i] : i for i in range(0, len(genomes))}
    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    # Use concurrent futures to get multiple rows at the same time
    # Then place completed rows into the matrix and update the row dictionary
    with ProcessPoolExecutor(max_workers=min(16, num_threads)) as ppe:
        for genome_name,temp_row in ppe.map(make_row, genomes, repeat(dataset),
        repeat(kmer_length)):
            for i, val in enumerate(temp_row):
                kmer_matrix[rows_dict[genome_name]][i] = val

    df = pd.DataFrame(data = kmer_matrix, columns = relevant_feats, index = runs)
    return df

def get_feat_splits(all_feats, RAM_denom):
    """
    Takes a list of features, splits it into a dictionary
    {0:[feat1,feat2,feat3], 1:[feat4,feat5......
    """
    assignments = []
    num_feats = len(all_feats)

    feats_per = num_feats//RAM_denom

    if num_feats % RAM_denom != 0:
        feats_per += 1

    for i in range(RAM_denom):
        assignments.append(all_feats[i*feats_per:(i+1)*feats_per])

    return assignments

def large_feat_select(file_paths, RAM_denom, num_feats, labels):
    """
    Takes a list of filepaths to slices of a feature matrix
    Does feature selection on a subset of data at a time as to not blow through the ram.

    If RAM_denom is set to 5, only 20% of each slice will be kept in ram.
    You are passing the denominator of a fraction, i.e. 1/2, 1/3, 1/4, 1/20
    i.e. if you pass 10, then it will reduce ram per slice to 1/10th

    Ram use will then be 0.2*(memory of 1 slice)*(number of slices loaded) + (memory of 1 slice)
    i.e., it will use more than 20% of the RAM as before, one part of the problem is reduced by 80%
    """

    f_values = []

    #  First load one slice, prime everything, then we will loop through the rest after
    slice1 = pd.read_pickle(file_paths[0])
    all_feats = slice1.columns
    del slice1

    feat_splits = np.array_split(all_feats, RAM_denom)

    for slice_num in range(RAM_denom):
        # as stated above, now we loop through the rest
        slices = []
        gc.collect()

        for i, file in enumerate(file_paths):
            if i == 0:
                slices = pd.read_pickle(file)[feat_splits[slice_num]]
            else:
                new_slice = pd.read_pickle(file)[feat_splits[slice_num]]
                slices = pd.concat([slices,new_slice])
                del new_slice
                gc.collect()

        f_vals, p_vals = f_classif(slices.values,labels)
        for i in range(len(f_vals)):
            f_values.append([f_vals[i],feat_splits[slice_num][i]])

    feats_to_keep = sorted(f_values, reverse=True)[:num_feats]
    feats_to_keep = [i[1] for i in feats_to_keep]

    slices = []
    for i, file in enumerate(file_paths):
        if i == 0:
            slices = pd.read_pickle(file)[feats_to_keep]
        else:
            new_slice = pd.read_pickle(file)[feats_to_keep]
            slices = pd.concat([slices,new_slice])
            del new_slice
            gc.collect()

    return slices


    # slices look like This
    """
                          AGCCTGGCTGGTCGGCTGTACAAAGACGAAT  ...  GTAGCATGAATGGGGGTAATCTGGAATGGAA
    SAMN02640777                                0  ...                                0
    SAMN02640827                                0  ...                                0
    SAMN02640878                                0  ...                                0
    SAMN02640928                                0  ...                                0
    SAMN02699192                                0  ...                                0
    ...                                       ...  ...                              ...
    SAMN05596337                                0  ...                                0
    SAMN05596788                                0  ...                                0
    SAMN05771752                                0  ...                                0
    SAMN05771803                                0  ...                                0
    SAMN05907766                                0  ...                                0

    [127 rows x 72504712 columns]
    """

def filter_feature(feature_data, labels, cv_num, num_of_cv, cv_df, hyp, attribute):
    """
    Takes full ~6300 genomes x 1 feature
    Returns only the kept sequences and their labels according to cv
    """
    training_samples = cv_df # needs to factor if hyp or not, check feat select to see if using mod+cv_num

    split_priority = [(i+int(cv_num)-1)%num_of_cv for i in range(num_of_cv)]
    """
    split_priority above defines which folds become validation or testing
    if fold 3 for a 5-fold cross-validation is given, split_priority will look like:
        [2, 3, 4, 0, 1]
    where fold 2 becomes the testing set, and 3 becomes validation (if requested, else becomes training)
    if fold 7 for 10-fold cross-validation is given, split_priority will look like:
        [6, 7, 8, 9, 0, 1, 2, 3, 4, 5]

    This means for cv we take 2,3,4,0 and for nested cv we take 2,3,4
    """
    keep_fold = split_priority[:3]

    mask = []
    pass #testing if other supervised split works first
    #return filtered_feats, filtered_labels

def build_variance_matrix(file_paths, RAM_denom, label_name, label_df, cv, hyp, type, trials, dataset_name):
    """
    Takes a list of filepaths to slices of a feature matrix
    Evaluates a subset of features at a time, under defined RAM contraints

    If RAM_denom is set to 5, only 20% of each slice will be kept in ram.
    You are passing the denominator of a fraction, i.e. 1/2, 1/3, 1/4, 1/20
    i.e. if you pass 10, then it will reduce ram per slice to 1/10th

    Ram use will then be 0.2*(memory of 1 slice)*(number of slices loaded) + (memory of 1 slice)
    i.e., it will use more than 20% of the RAM as before, one part of the problem is reduced by 80%

    For each attribute*trial*fold*slice, a df is stored that looks like:
        f_val,    p_val,    feature
    0   0.01      0.01      GTAGCATGAATGGGGGTAATCTGGAATGGAA
    1
    2
    """
    files = os.listdir("data/{}/wgs/raw/".format(dataset_name))
    files = [i.split('.')[0] for i in files]
    label_df = label_df[~label_df.index.duplicated(keep='first')]
    label_df = label_df.reindex(files) # label_df now only has the 6328 rows

    f_values = []

    #  First load one slice, prime everything, then we will loop through the rest after
    slice1 = pd.read_pickle(file_paths[0])
    all_feats = slice1.columns

    del slice1

    feat_splits = np.array_split(all_feats, RAM_denom)

    for slice_num in range(RAM_denom):
        # as stated above, now we loop through the rest
        slices = []
        gc.collect()

        for i, file in enumerate(file_paths):
            if i == 0:
                slices = pd.read_pickle(file)[feat_splits[slice_num]]
            else:
                new_slice = pd.read_pickle(file)[feat_splits[slice_num]]
                slices = pd.concat([slices,new_slice])
                del new_slice
                gc.collect()

        # At this point we have a full matrix for x amount of features. Called slices
        """
                f0, f1, f2, .... f128
        genome1
        genome2
        genome6328mask = pd.read_pickle("data/{}/features/masks/{}_{}.df".format(dataset_name,type,label_name))[attribute]
        """
        already_done = 0
        new = 0
        for trial in range(trials):
            split_df = pd.read_pickle("data/{}/splits/split{}_{}_{}_{}xCV.df".format(dataset_name,trial,type,label_name,cv))
            split_df = split_df.reindex(files) #split now only has 6328 samples
            for attribute in split_df.columns:
                labels = label_df[attribute]
                mask = pd.read_pickle("data/{}/features/masks/{}_{}.df".format(dataset_name,type,label_name))[attribute]
                attr_slices, labels = supervised_model.apply_mask(slices, labels, mask)
                for fold in range(cv):
                    out_path = "data/{}/variance/slice{}_of_{}_trial={}_type={}_label={}_attribute={}_fold={}_of_{}xCV.df".format(
                    dataset_name,slice_num+1,RAM_denom,trial,type,label_name,attribute,fold,cv)

                    if os.path.isfile(out_path):
                        already_done+=1
                    else:
                        x_train, y_train, x_test, y_test = supervised_model.split_data(attr_slices, labels, split_df, attribute, hyp, fold, cv, False)
                        f_vals, p_vals = f_classif(x_train, y_train)

                        for deletable in [x_train, y_train, x_test, y_test]:
                            del deletable
                        gc.collect()

                        df = pd.DataFrame(data = list(zip(f_vals, p_vals, attr_slices.columns)),columns=['f_val','p_val','feature'])
                        df.to_pickle(out_path)

                        new+=1
        print("{} variance files created, {} were already completed".format(new, already_done))

    return 0

def manual_variance_call(slice_num, trial, attribute):
    """
    This functions as a way to call the
    build_variance_matrix() above manually from the command line. Please dont use
    """
    #file_paths, RAM_denom, label_name, label_df, cv, hyp, type, trials, dataset_name
    #save_path = "slice{}_of_5_trial={}_type=31mer_label=AMR_MIC_attribute={}_fold={}_of_5xCV.df".format(slice,trial,attribute,i)
    if np.sum([os.path.isfile("slice{}_of_5_trial={}_type=31mer_label=AMR_MIC_attribute={}_fold={}_of_5xCV.df".format(slice_num+1,trial,attribute,i)) for i in range(5)]) ==5:
        print("already created these files")
        sys.exit()
    file_paths = ["data/salm_amr/wgs/master_31mers/sub_df_{}_of_50.df".format(i+1) for i in range(50)]
    RAM_denom = 5
    label_name = 'AMR_MIC'
    label_df = pd.read_pickle("data/salm_amr/labels/AMR_MIC.df")
    cv = 5
    hyp = False
    type = '31mer'
    trials = 1
    dataset_name = 'salm_amr'

    files = os.listdir("data/{}/wgs/raw/".format(dataset_name))
    files = [i.split('.')[0] for i in files]
    label_df = label_df[~label_df.index.duplicated(keep='first')]
    label_df = label_df.reindex(files) # label_df now only has the 6328 rows

    f_values = []

    #  First load one slice, prime everything, then we will loop through the rest after
    slice1 = pd.read_pickle(file_paths[0])
    all_feats = slice1.columns

    del slice1

    feat_splits = np.array_split(all_feats, RAM_denom)


    slices = []
    gc.collect()

    for i, file in enumerate(file_paths):
        if i == 0:
            slices = pd.read_pickle(file)[feat_splits[slice_num]]
        else:
            new_slice = pd.read_pickle(file)[feat_splits[slice_num]]
            slices = pd.concat([slices,new_slice])
            del new_slice
            gc.collect()

    split_df = pd.read_pickle("data/{}/splits/split{}_{}_{}_{}xCV.df".format(dataset_name,trial,type,label_name,cv))
    split_df = split_df.reindex(files) #split now only has 6328 samples

    labels = label_df[attribute]
    mask = pd.read_pickle("data/{}/features/masks/{}_{}.df".format(dataset_name,type,label_name))[attribute]
    attr_slices, labels = supervised_model.apply_mask(slices, labels, mask)
    for fold in range(cv):
        out_path = "data/{}/variance/slice{}_of_{}_trial={}_type={}_label={}_attribute={}_fold={}_of_{}xCV.df".format(
        dataset_name,slice_num,RAM_denom,trial,type,label_name,attribute,fold,cv)

        x_train, y_train, x_test, y_test = supervised_model.split_data(attr_slices, labels, split_df, attribute, hyp, fold, cv, False)
        f_vals, p_vals = f_classif(x_train, y_train)

        for deletable in [x_train, y_train, x_test, y_test]:
            del deletable
        gc.collect()

        df = pd.DataFrame(data = list(zip(f_vals, p_vals, attr_slices.columns)),columns=['f_val','p_val','feature'])
        df.to_pickle(out_path)


    return 0


if __name__ == "__main__":
    raise Exception("Please use the acheron commands to call this module")
    # if you know what you are doing, comment out the exception above
    # then run the command
    """
    do sbatch -c 1 --mem 700G --partition NMLResearch --job-name acheron --wrap="python acheron/workflows/over_11mer_matrix.py $slice $trial $drug"
    """
    manual_variance_call(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
