run2bio = pd.read_excel("data/old_ncbi_salm/NCBI_dataset_metadata.xlsx")

def get_bio(run):
    return run2bio[run2bio['run'] == run]['biosample'].iloc[0]

for file in glob.glob("data/ncbi_salm/wgs/raw/SRR*"):
     try:
         os.system("cp {} data/ncbi_salm/wgs/raw/{}.fasta".format(file,get_bio(file.split('.')[0].split('/')[-1])))
     except:
         print(file)
