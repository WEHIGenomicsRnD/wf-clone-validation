#!/usr/bin/env python

import argparse
import os
import pandas as pd
import glob
import time

from .util import get_named_logger, wf_parser
#parser = argparse.ArgumentParser()
#parser.add_argument('-f','--anal_folder', help='Path to analysis folder')
time.sleep(5)

def main(args):

   anal_path=args.anal_folder
   ## Parsing Histogram files to get the highest peak
   fname=anal_path.strip().split("/")[-2]

#   qc_rules="/vast/scratch/users/gupta.i/forked_wf/wf-clone-validation/qc_rules.txt"
   qc_rules=args.rule_file
   rules_df=pd.read_csv(qc_rules,header=None,delimiter="\t")

   outfile=f"sample_QC.txt"
   out=open(outfile,'w')
   out.write("SampleId\tSize(User derfined)\tSize (Raw fastq)\tSize (After Assmbly)\n")

   hist_dict={}
   hist_files=f"{anal_path}/*.length_hist"
   for file in args.hist_file:
      fname=os.path.splitext(os.path.basename(file))[0]
      df=pd.read_csv(file,header=None,delimiter="\t",skiprows=200)
      df2=df.loc[df[2].idxmax()][1]
      lower=round(df2 * .9,2)
      upper=round(df2 * 1.2,2)
      hist_dict[fname]=[lower,upper,df2]

   print(f"Hist {hist_dict}")
   ## Parsing assembly length file

   af_dict={}
   assm_file=args.stats

   asm_df=pd.read_csv(assm_file,header=0)
   for r in range(len(asm_df)):
      size=asm_df['Length'][r]
      reason = asm_df['Assembly completed / failed reason'][r]
      if 'nan' in str(size):
        size=0
      af_dict[asm_df['Sample'][r]]=[size,reason]

#   print(f"AF-- {af_dict}")
   ## Parsing user define length file

   user_dict={}
   user_file=f"{anal_path}/sample_sheet.csv"

   user_df=pd.read_csv(user_file,header=0)
   for r in range(len(user_df)):
       size= user_df['approx_size'][r]
       lower= size - 20
       upper= round(size * 1.2,2)
       user_dict[user_df['alias'][r]]=[lower,upper,size]

#   print(f"user: {user_dict}")

   ## Checking the rules for each sample

   res={}
   for k,v in user_dict.items():
       if k not in hist_dict:
           hist_dict[k]=[0,0,0]
       x=af_dict[k][0]
       reason=af_dict[k][1]
       test_lower=user_dict[k][0]
       test_upper=user_dict[k][1]
       exact = user_dict[k][2]
       peak_lower=hist_dict[k][0]
       peak_upper=hist_dict[k][1]
#       print(f"Rules {k} -- {test_lower} -- {test_upper}--{x} -- {peak_lower}-- {peak_upper} ")
       for r in range(len(rules_df)):
          rule_c1=rules_df[1][r]
          rule_c2=rules_df[2][r]
#          print(f"Rules -- {rule_c1}   --- {rule_c2}-- {test_lower} -- {test_upper} -- {x} -- {peak_lower}-- {peak_upper}")
          if eval(rule_c1) and eval(rule_c2):
             status=rules_df[3][r]
             if reason not in 'Completed successfully':
                 status=reason
             res[k]=[k,user_dict[k][2],hist_dict[k][2],af_dict[k][0],status]
             break
#          else:
#             print(f"No rule for sample {k}")

#print(res)
##Writing results to file
   for k,v in res.items():
      out.write("\t".join(str(rv) for rv in v) + '\n')

   out.close()


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("eval_qc")
#    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--anal_folder', help='Path to analysis folder')
    parser.add_argument('-q','--rule_file', help='QC rule file')
    parser.add_argument('-s','--stats', help='Stats file')
    parser.add_argument('-r','--hist_file', help='Histogram file',nargs='+')
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
