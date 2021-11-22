
import glob
import pandas as pd
import subprocess

def extract_value(sacct_str, sacct_type):
    sacct_str = sacct_str.replace(' ','')

    # NOTE: Resident Set Size is measure in kibibytes, not kilobytes,
    # despite using the kB suffix, rather than the kiB

    if sacct_type == "RAM":
        # check if last char is a SI suffix, else assume bytes
        # always return in GiB for consistency
        # (where 1024 MiB == 1 GiB) {'GiB': 1024**3,'MiB': 1024**2,'KiB': 1024}

        last_char = sacct_str[-1]
        if last_char.isalpha():
            ram_int = int(float(sacct_str[:-1]))
            if last_char == 'K':
                return ram_int/(1024**2)
            elif last_char == 'M':
                return ram_int/1024
            elif last_char == 'G':
                return ram_int
            else:
                raise exception("Suffix {} was unexpected in RAM calculation".format(last_char))

        else:
            return int(sacct_str)/(1024**3)


    if sacct_type == "TIME":
        times = sacct_str.split(":")
        if '-' in times[0]:
            days, hours = times[0].split('-')
            times[0] = int(days)*24+int(hours)
        times = [int(i) for i in times]
        total_mins = times[0]*60+times[1]+times[2]/60
        return total_mins


def add_resources_to_summaries():
    """
    This script scans for results that have slurm JOBID's and adds in RAM and
    time usage
    """
    # to get as specific as you want, change the following superstars to your search
    # for example, change model to XGB to only include those
    model = '*'
    train = '*'
    test = '*'
    validate = '*'
    feats = '*'
    type = '*'
    hyp = '*'
    cvfolds = '*'
    attribute = '*'
    trial = '*'

    search_space = "results/model={}_train={}_test={}_validate={}_feats={}_type={}_hyp={}_cvfolds={}_attribute={}_trial={}/summary.df".format(
        model, train, test, validate, feats, type, hyp, cvfolds, attribute, trial)


    file_list = glob.glob(search_space)
    print("Checking {} results files".format(len(file_list)))

    num_slurm = 0

    for file in file_list:
        try:
            summary_df = pd.read_pickle(file)
        except AttributeError:
            print("Ensure the version of Pandas used to pickle is the same as the version used to read the pickle")
            raise

        if 'SLURM_JOBID' in summary_df.columns and 'Time (m)' not in summary_df.columns:
            num_slurm += 1
            slurm_id = summary_df['SLURM_JOBID'][0]

            elapsed = subprocess.getoutput(["sacct -j {} --format='Elapsed' | tail -n 1".format(slurm_id)])
            max_ram = subprocess.getoutput(["sacct -j {} --format='MaxRSS' | tail -n 1".format(slurm_id)])

            # time is expressed in minutes, max ram in Gibibytes (GiB)
            try:
                elapsed = extract_value(elapsed, 'TIME')
            except:
                print("Time unextractable from slurm job with id {}".format(slurm_id))
                raise

            try:
                max_ram = extract_value(max_ram, 'RAM')
            except:
                print("RAM unextractable from slurm job with id {}".format(slurm_id))
                raise

            summary_df['Time (m)'] = elapsed
            summary_df['Max Ram (GiB)'] = max_ram

            summary_df.to_pickle(file)

        else:
            # This section intentionally left blank, if no slurm info is available,
            # we arent going to touch the summary file
            pass
    return num_slurm

if __name__ == "__main__":
    num = add_resources_to_summaries()
    print("RAM and Time use has been added to {} results files".format(num))
