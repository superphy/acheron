"""
The type of the feature the models were trained on were forgotten when included
in the results directories. These are all 11mers so we just need to add _type=11mer
to the correct spot

NB: dataset names and other variables will have _ in them, so cant split on that parameter
"""

from glob import glob
import os, sys

def main():
    dirs = glob("results/*")

    for dir_str in dirs:
        feats_loc = dir_str.find('feats')
        post_feats_underscore = dir_str.find('_',feats_loc,-1)

        # This takes the dir name from the start, until after the number of features
        # 'results/model=XGB_train=salm_amr_test=grdi_validate=none_feats=1000'
        first_half = dir_str[0:post_feats_underscore]

        # This is everything after the feats
        # 'hyp=False_cvfolds=5_attribute=AMP_trial=0'
        second_half = dir_str[post_feats_underscore+1:]

        if 'type'  in second_half:
            continue
        final_dir = first_half + "_type=11mer_" + second_half

        os.system("mv {} {}".format(dir_str,final_dir))

if __name__ == "__main__":
    main()
