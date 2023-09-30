import os
import glob

import pandas as pd
import numpy as np

script_path = os.getcwd()
annotation_path = os.path.join(script_path, "annotations")
PWM_path = os.path.join(script_path, "PWMs") # proteases with 8 and more cleavages

seq_file = input("Enter the name of sequence (FASTA) file: ")

''' PWM score function '''
def get_PWM_score(sequence, PWM):
    
    processed_sequence = f"--{sequence}--"
    PWM_scores = []
    for i in range(2, len(sequence) + 1):
        seq_frame = processed_sequence[i - 2: i + 4]
        local_score = 0
        for aa, pos in zip(seq_frame, ["P3", "P2", "P1", "P1'", "P2'", "P3'"]):
            if aa == '-': continue
            local_score += PWM.loc[PWM["AA"] == aa, pos].values[0]
        PWM_scores.append(round(local_score, 4))  
    
    PWM_scores.append(np.nan)
    return PWM_scores

''' Sequence extracting '''
with open(seq_file, 'r') as file:
    data = file.read()
    seq = ''.join(data.strip('\n').split('\n')[1:])# FASTA-file: The first string - the name and additional information about protein sequence #

result = {'AA':list(seq), 'num_AA':list(range(1, len(seq) + 1))}

PWMs = glob.glob(os.path.join(PWM_path, "*_PWM.txt"))
for index, pwm_file in enumerate(PWMs):
    merops_id = os.path.basename(pwm_file).split('.txt')[0].split('_')[0]
    print(f"{index + 1} --- {merops_id}")
    
    pwm = pd.read_csv(pwm_file)
    PWM_scores = get_PWM_score(seq, pwm)
    result[merops_id] = PWM_scores

''' Matrix transforming '''
df = pd.DataFrame(result)
df["AA"] = df["AA"] + '.' + df["num_AA"].astype(str)
del df["num_AA"]
df = df.set_index("AA")
df = df.T.reset_index().rename(columns={"index":"MEROPS_code"})

''' Annotating with name from MEROPS data '''
annotation_df = pd.read_csv(os.path.join(annotation_path, "HomoSapiens_MEROPS_proteases.tsv"), sep='\t').rename(columns={"Code":"MEROPS_code", "Name":"MEROPS_name"})
df = df.merge(annotation_df[["MEROPS_code", "MEROPS_name"]], on="MEROPS_code", how="left").drop_duplicates().reset_index(drop=True)
df.loc[df["MEROPS_code"] == "S01.247", "MEROPS_name"] += "/TMPRSS2" # append TMPRSS2 name #
df = df[["MEROPS_code", "MEROPS_name"] + [i for i in df.columns if '.' in i]] 
df.to_csv(input("Enter the name of output file (without extension): ") + '_result.csv', index=False)