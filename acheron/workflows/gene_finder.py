import os,sys
import pandas as pd
import numpy as np
import skbio.io
import gffpandas.gffpandas as gffpd
from statistics import stdev

def find_locs(kmer, blast_df):
    """
    Finds the start and stop locations of this k-mer in each genome
    """
    locs = []

    # filter blast results to just our kmer of interest
    kmer_df = blast_df[blast_df['qseqid'] == kmer]

    if kmer_df.shape[0]==0:
        raise Exception("The k-mer {} was not found in the blast search".format(kmer))

    # for every kmer hit in the genome, append the location
    for i in range(kmer_df.shape[0]):
        send = kmer_df['send'].values[i]
        sstart = kmer_df['sstart'].values[i]
        genome_id = kmer_df['sseqid'].values[i].split('_')[0]
        contig_name = kmer_df['sseqid'].values[i]
        locs.append([send,sstart,genome_id,contig_name])

    # locs is 2d list, each row is (start, stop, genome_id, contig_name)
    return locs

def find_gene(start, stop, genome_id, contig_name, prokka_loc):
    """
    Finds the nearest gene upstream and downstream of the k-mer,
    reports their distance using the prokka annotations.
    """
    gene_up = ''
    dist_up = -1
    gene_down = ''
    dist_down = -1

    # prokka renames contigs but the numbers are consistent, so we need to pull the number
    if("NODE" in contig_name):
        orig_contig_name = contig_name.split('|')
        assert(len(orig_contig_name)==2)
        orig_contig_name = orig_contig_name[1]
        contig_num = orig_contig_name.split('_')[1]

    elif(contig_name.split('_')[0] == genome_id and len(contig_name.split('_'))==2):
        contig_num = contig_name.split('_')[1]

    # SRR5573065_SRR5573065.fasta|33_length=41292_depth=1.01x
    elif(contig_name.split('_')[0] == genome_id and len(contig_name.split('_')) in [4,5]):
        contig_num = contig_name.split('|')[1].split('_')[0]

    else:
        raise Exception("Unexpected contig name found: {}".format(contig_name))

    if(prokka_loc[-5:-1]=='ncbi'):
        gff_loc = "annotation/annotated_genomes/"
    else:
        gff_loc = "annotation/annotated_grdi_genomes/"

    # scan through contigs until the correct number is found, then keep the contig name
    with open("{0}{1}/{1}.gff".format(gff_loc,genome_id)) as file:
        for line in file:
            if("_{} ".format(contig_num) in line):
                prokka_contig = line.split(" ")[1]
                break

    if('prokka_contig' not in locals()):
        print("Contig number {2} and contig name {3} not located in {0}{1}/{1}.gff".format(gff_loc,genome_id, contig_num, contig_name))
        return [gene_up, dist_up, gene_down, dist_down]

    # columns are: ['seq_id','source','type','start','end','score','strand','phase','attributes']
    #with open(prokka_loc+genome_id+'.pkl', 'rb') as fh:
    #    gff_df = skbio.io.read(fh, format="blast+6",into=pd.DataFrame,default_columns=True)
    gff_df = pd.read_pickle(prokka_loc+genome_id+'.pkl')

    # keep only the contig the kmer was found in and only show coding sequences (Prodigal)
    contig_df = gff_df[gff_df['seq_id']==prokka_contig]
    contig_df = contig_df[contig_df['type']=='CDS']

    start = int(start)
    stop = int(stop)

    df_length = contig_df.values.shape[0]

    # find the nearest gene/genes
    # for every gene found by prokka, does it contain the kmer or is it near?
    for gene_num, gene_anno in enumerate(contig_df.values):
        gene_start = int(gene_anno[3])
        gene_stop = int(gene_anno[4])

        try:

            if(start > gene_stop):
                if(gene_num==df_length-1):
                    # we are after the last gene
                    gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                    dist_up = start - gene_stop
                    if(gene_dict['product']=='hypothetical protein'):
                        gene_up = "HypoProt:hypothetical protein"
                    else:
                        gene_up = gene_dict['gene']+':'+gene_dict['product']
                    break

                # we are not at the correct gene yet
                continue
            elif(stop < gene_start):
                # the kmer is before the current gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_down = gene_start-stop
                if(gene_dict['product']=='hypothetical protein'):
                    gene_down = "HypoProt:hypothetical protein"
                else:
                    try:
                        gene_down = gene_dict['gene']+':'+gene_dict['product']
                    except KeyError:
                        gene_down = 'none:'+ dict(i.split('=') for i in gene_anno[8].split(';'))['product']

                prev_gene_anno = contig_df.values[gene_num-1]

                gene_dict = dict(i.split('=') for i in prev_gene_anno[8].split(';'))
                dist_up = start - prev_gene_anno[4]
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            elif(start >= gene_start and stop <= gene_stop):
                # the kmer is inside of a gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_up = 0
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            elif(start <= gene_stop <= stop):
                # kmer hanging over right end of a gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_up = stop-gene_stop
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            elif(start <= gene_start <= stop):
                # kmer hanging over the left end of a gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_up = gene_start-start
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            else:
                raise Exception("Unexpected kmer start: {} stop: {} in relation to gene start: {} stop: {}".format(start, stop, gene_start, gene_stop))
        except KeyError:
            gene_up = 'none:'+ dict(i.split('=') for i in gene_anno[8].split(';'))['product']
            break

    return [gene_up, dist_up, gene_down, dist_down]

def find_hits(kmers,blast_search):
    """
    Input:
    list of k-mers to search for, preformed 1 at a time
    A PATH TO a blast search of those k-mers through your genomes

    Output:
    A pandas df where each row shows where a kmer was found in a genome,
    Which gene is upstream and downstream and how far away they are.
    Will have (# of kmers * number of genomes) rows
    """

    with open(blast_search) as fh:
        blast_df = skbio.io.read(fh, format="blast+6",into=pd.DataFrame,default_columns=True)

    # each row in gene_hits will be [kmer, gene_up, dist_up, gene_down, dist_down, start, stop, genome_id, contig_name]
    gene_hits = []
    for kmer in top_feats:
        # locs is 2d list, each row is (start, stop, genome_id, contig_name)
        locs = gene_finder.find_locs(kmer, blast_df)
        for loc in locs:
            # gene_info is 1D list: gene_up, dist_up, gene_down, dist_down
            gene_info = gene_finder.find_gene(*loc, prokka_loc) #TODO prokka loc as "annotation/gffpandas_ncbi/"
            gene_hits.append([kmer]+gene_info+loc)

    hits_df = pd.DataFrame(data = gene_hits,columns = ['kmer', 'gene_up', 'dist_up', 'gene_down', 'dist_down', 'start', 'stop', 'genome_id', 'contig_name'])

    return hits_df

def merge_df(df_path, drug, dataset):

    """
    Takes path to pandas df with the columns:
    [kmer, gene_up, dist_up, gene_down, dist_down, start, stop, genome_id, contig_name]
    and returns a df with the columns:
    [drug, dataset, kmer, gene_up, gene_up_count, avg_dist_up, gene_down, gene_down_count, avg_dist_down]
    """
    df = pd.read_pickle(df_path)
    hit_summary = []

    for kmer in set(df['kmer']):
        # filter for only a single kmer
        kmer_df = df[df['kmer']==kmer]

        for gene_direc, gene_dist in [['gene_up','dist_up'],['gene_down','dist_down']]:
            for gene in set(kmer_df[gene_direc]):
                if(len(gene)==0):
                    continue
                # filter for only a single gene
                gene_df = kmer_df[kmer_df[gene_direc]==gene]

                total_dist = 0

                for dist in gene_df[gene_dist]:
                        total_dist += abs(float(dist))

                count = gene_df.shape[0]
                average_dist = total_dist/count
                if(len(gene_df[gene_dist])<2):
                    std_dev = 0
                else:
                    std_dev  = stdev([abs(int(float(i))) for i in gene_df[gene_dist]])

                try:
                    gene, name = gene.split(':')
                except:
                    print("Gene: {}".format(gene))
                    print("{} {}".format(drug,dataset))
                    gene, carb, name = gene.split(':')

                hit_summary.append([dataset, drug, kmer, gene, name, count, average_dist, std_dev])

    return hit_summary


def hit_summary(dataset,out_path):
    drugs = ['AMP','AMC','AZM','CHL','CIP','CRO','FIS','FOX','GEN','NAL','SXT','TET','TIO','STR','KAN']

    all_hits = []
    for drug in drugs:
        if(dataset == 'grdi' and drug in ['FIS']):
            continue
        df_path = "annotations/{}_{}/hits_df.pkl".format(dataset, drug)

        drug_list = merge_df(df_path, drug, dataset)
        for hit in drug_list:
            all_hits.append(hit)

    data = np.asarray(all_hits)

    what_amg = np.zeros((data.shape[0]), dtype = object)

    amgs = pd.read_csv("data/gene_labels.tsv",sep='\t')

    for i, val in enumerate(what_amg):
        for j, amg_list in enumerate(amgs['gene_codes']):
            for amg in amg_list:
                if(data[i][3].split('_')[0] in amg):
                    what_amg[i] = amgs['AMR Gene Family'][j]

    all_df = pd.DataFrame(data=np.concatenate((data,np.asarray([what_amg]).T),axis=1),columns=['dataset','drug','kmer','gene','name','count','average_dist', 'std_dev' ,"AMG"])

    all_df.to_csv(output)

def score_summary(type, dataset, top_feats):
    """
    where top_feats is a 1d array
    """
    for drug in ['AMC','AMP','AZM','FOX','TIO','CRO','CHL','CIP','GEN','NAL','STR','FIS',
    'TET','SXT','KAN']:
        if dataset == 'grdi' and drug in ['FIS']:
            continue

        # TODO: This function put on hold, need f_classif scores, imp arr is kmer,XGBoost_score, f_classif
        # Will need to figure out how storing f_classif first
        imp_path = "data/multi-mer/feat_ranks/{}_1000000_{}_{}mer_feature_ranks.npy".format(dataset,drug,kmer_length)
        imp_arr = np.load(imp_path, allow_pickle=True)

        for kmer in top_feats:
            indx = np.where(imp_arr[0]==kmer)
            scores = imp_arr[:,indx]
            row.append([dataset,drug,kmer,scores[1],scores[2]])

    df = pd.DataFrame(data=row,columns=['Dataset','Antimicrobial',kmer_length+'mer','XGBoost Score','f_classif'])
    df.to_csv("data/multi-mer/feat_ranks/{}mer_score_summary.csv".format(kmer_length))
