import os
import glob
import argparse

import pandas as pd
import numpy as np

''' Define required paths '''
script_path = os.getcwd()
annotation_path = os.path.join(script_path, "annotations")
PWM_path = os.path.join(script_path, "PWMs") # proteases with 8 and more cleavages

''' Define program parameters '''
parser = argparse.ArgumentParser(
                                prog = "PWM-scores Calculation",
                                description = "Calculation of PWM-scores for your protein sequence",
                                epilog = "For more information, see https://github.com/KazanovLab/ProteaseSpecificityModels"
                                )

parser.add_argument('path_to_FASTA_file', nargs='+', type=str, help="Enter the whole path to your FASTA file")
args = parser.parse_args()

''' Define PWM score calculation function '''
def get_PWM_score(sequence, PWM):

    frame = ["P3", "P2", "P1", "P1'", "P2'", "P3'"]
    processed_sequence = f"--{sequence}--"
    start = 2
    finish = 4
        
    PWM_scores = []
    for i in range(start, len(sequence) + start - 1):
        seq_frame = processed_sequence[i - start: i + finish]
        local_score = 0
        for aa, pos in zip(seq_frame, frame):
            if aa == '-': continue
            local_score += PWM.loc[PWM["AA"] == aa, pos].values[0]
        PWM_scores.append(round(local_score, 4))  
    
    PWM_scores.append(np.nan)
    return PWM_scores

''' Define list of FASTA files '''
FASTA_list = []
for file in args.path_to_FASTA_file:
    if '*' in file:
        FASTA_list.extend(glob.glob(file))
    else:
        FASTA_list.append(file)
FASTA_list = list(set(FASTA_list))

''' Calculation of PWM-scores '''
for FASTA_file in FASTA_list:
    print('#'*40 + f" {FASTA_file} " + '#'*40)

    ''' Extract sequence from FASTA file '''
    with open(FASTA_file, 'r') as file:
        data = file.read()
        seq = ''.join(data.strip('\n').split('\n')[1:])# FASTA-file: The first string - the name and additional information about protein sequence, the second and next strings - the sequence #
        print(f"\nThe sequence extracted:\n{seq}")
        print(f"\nThe length of sequence extracted:\n{len(seq)}")
    
    result = {'AA':list(seq), 'num_AA':list(range(1, len(seq) + 1))}
    PWMs = glob.glob(os.path.join(PWM_path, "*_PWM.txt"))
    print(f"\nCalculation of PWM scores:")
    
    ''' For each PWM '''
    for index, pwm_file in enumerate(PWMs):
        merops_id = os.path.basename(pwm_file).split('_')[0]
        print(f"{index + 1} --- {merops_id}")
        
        ''' Calculation PWM-scores '''
        pwm = pd.read_csv(pwm_file)
        PWM_scores = get_PWM_score(seq, pwm)
        result[merops_id] = PWM_scores
    print(f"Calculation of PWM scores is well done!")

    ''' Matrix transforming '''
    df = pd.DataFrame(result)
    df["AA"] = df["AA"] + '.' + df["num_AA"].astype(str)
    del df["num_AA"]
    df = df.set_index("AA")
    df = df.T.reset_index().rename(columns={"index":"MEROPS_code"})

    ''' Annotating with name from MEROPS data '''
    print('\nAnnotating with MEROPS, TISSUES and COMPARTMENTS data...')
    protease_classes = {"A":"Aspartic class", "C":"Cysteine class", "G":"Glutamic class", "M":"Metallo class", "N":"Asparagine class", "P":"Mixed class", "S":"Serine class", "T":"Threonine class", "U":"Unknown class", "X":"Compound Peptidase class"}
    MEROPS_df = pd.read_csv(os.path.join(annotation_path, "HomoSapiens_MEROPS_proteases.tsv"), sep='\t').rename(columns={"Code":"MEROPS_code", "Name":"MEROPS_name"})
    df = df.merge(MEROPS_df[["MEROPS_code", "MEROPS_name"]], on="MEROPS_code", how="left").drop_duplicates().reset_index(drop=True)
    df['Protease_class'] = df['MEROPS_code'].apply(lambda x: protease_classes[x[0]] if x[0] in protease_classes else '-')
    df.loc[df["MEROPS_code"] == "S01.247", "MEROPS_name"] += "/TMPRSS2" # append TMPRSS2 name #

    ''' Annotating with TISSUES and ProteinAtlas data '''
    TISSUES_data = pd.read_csv(os.path.join(annotation_path, "TISSUES_ProteinAtlas_MEROPS_HumanProteases.csv"))
    df = df.merge(TISSUES_data[["MEROPS_code", "Respiratory system", "Lung", "Gastrointestinal tract"]], on="MEROPS_code", how="left").drop_duplicates()

    ''' Annotating with COMPARTMENTS data '''
    COMPARTMENTS_data = pd.read_csv(os.path.join(annotation_path, "COMPARTMENTS_MEROPS_HumanProteases.csv"))
    COMPARTMENTS_values = ["Extracellular space", "Cytoplasm", "Plasma membrane", "Lysosome", "Endosome", "Golgi apparatus"]
    df = df.merge(COMPARTMENTS_data[["MEROPS_code"] + COMPARTMENTS_values], on="MEROPS_code", how="left").drop_duplicates()

    for col in COMPARTMENTS_values + ["Lung", "Respiratory system", "Gastrointestinal tract"]:
        df[col].fillna(0, inplace=True)
        df[col] = df[col].astype(int)
    df = df[["MEROPS_code", "Protease_class", "MEROPS_name"] + ["Lung", "Respiratory system", "Gastrointestinal tract"] + COMPARTMENTS_values + [i for i in df.columns if '.' in i]] 
    matrix = df[[i for i in df.columns if '.' in i] + ['MEROPS_code']].set_index('MEROPS_code')
    print('Annotating with MEROPS, TISSUES and COMPARTMENTS data is well done!')
    
    result_file = FASTA_file.split('.fasta')[0] + '.xlsx'
    matrix_file = FASTA_file.split('.fasta')[0] + '.matrix.xlsx'
    df.to_excel(result_file, index=False)
    matrix.to_excel(matrix_file)
    print(f"\nResults are written in '{result_file}'.")
    print(f"\nMatrix are written in '{matrix_file}'.")