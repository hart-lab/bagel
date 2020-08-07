#!/usr/bin/python

# a script for run BAGEL without copynumber artefact correction
# usage) python3 run_bagel.py -i [readcount table] -o [output file tag] -s [screen info]

import sys
import numpy as np
import pandas as pd
import subprocess

from os import path
import click



'''
### before run it ###
- for multi-targeting effect correction, download precalculated alignment file or build up using 'precalc_library_alignment_info.py'. 
'''

### modify if needed ###

# BAGEL path (Python 3)
BAGEL_PATH = "./BAGEL.py"

# Reference Essential file
REF_ESS = "CEGv2.txt"

# Reference Non-essential file
REF_NONESS = "NEGv1.txt"

# Screen quality check
QC_PATH = "qc_ess_dist.py"


## enable multi-targeting correction
MULTI_TARGETING_CORRECTION = True

# Precalculated genome alignment of sgRNAs (if multi targeting correcetion is enabled) via precalc_library_alignment_info.py
'''
Library
TKOv3 (Toronto, Reagent ID)
  -> TKOv3_align_summary_REAGENT_ID.txt
Avana (DepMap, sgRNA sequence)
  -> Avana_align_summary.txt
KYv1.0 (Project Score):
  -> KYv1_align_summary.txt
'''
GUIDE_ALIGNMENT_LIBRARY = "TKOv3_align_summary_REAGENT_ID.txt"




## Specific options

# Bootstrapping
BOOTSTRAPPING = False
BOOTSTRAPPING_ITERATION = 1000

# Normalize sgRNA counts per gene (default: 0 -> no normalization)
NORMALIZE_SGRNA_COUNT = 0  # NORMALIZE_SGRNA_COUNT = [Target number]

# Normalize the number of replicates (default: 0 -> no nomalization)
NORMALIZE_REPLICATES_N = 0  # NORMALIZE_REPLICATES_N = [Target number]







########################


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input-readcount', required=True, type=click.Path(exists=True))
@click.option('-s', '--screen-info', required=True, type=click.Path(exists=True))
@click.option('-o', '--output-prefix')
def run_bagel_script(
        input_readcount, screen_info,output_prefix
):
    print("Read count file  :", input_readcount)
    print("Screen info file :", screen_info)
    if output_prefix == None:
        output_prefix = input_readcount
    print("Output prefix :", output_prefix)

    readcount_data = pd.read_csv(input_readcount,index_col=0,header=0,sep="\t")
    readcount_data.index.rename('sgRNA',inplace=True)  # for CRISPRcleanR
    readcount_data.rename(columns={"Gene":'gene',"GENE":'gene'},inplace=True)  # for CRISPRcleanR
    readcount_data.head()

    screeninfodata = pd.read_table(screen_info,header=None,index_col=None,comment='#')
    screeninfodata.head(3)
    
    ## simple file check

    # Unique ids
    if readcount_data.index.is_unique == False:
        print("\n\n\nError! There are duplicated GUIDE IDs")
        sys.exit(1)

    # The number of of clumns
    if len(screeninfodata.columns) != 3:
        print(""" Screen information error!
        #All delimiter is tab
        #*: replicates/controls should be identical as headers of the readcount file.
        #[Screen Name/ID]   [replicates*]   [controls*]
        BATCH1_HAP1_T3  T3_A,T3_B,T3_C  T0_-
        BATCH1_HAP1_T18 T18_A,T18_B,T18_C   T0_-
        BATCH1_HAP1_T18_STARVED T18_A_Starved,T18_B_Starved,T18_C_Starved   T0_-""")
        sys.exit(1)

    # Bulid table to store file names

    filenames_rep = pd.DataFrame(columns = ['screenid', 'readcount','fc_prefix','fc','controls'])
    filenames_screen = pd.DataFrame(columns = ['combinedFC','BF','PR','replicates'])


    for i in screeninfodata.index:
        fc_prefix_screen = output_prefix + "." + screeninfodata[0][i] 
        fcfile = fc_prefix_screen + ".combined_fc"
        bffile = fc_prefix_screen + ".bf"
        prfile = fc_prefix_screen + ".pr"
        screen_id = screeninfodata[0][i]
        controls = screeninfodata[2][i]

        for rep in screeninfodata[1][i].split(","):
            fc_prefix_rep = fc_prefix_screen + "." + rep.replace(" ","_").replace("/","_")
            fcfile_rep = fc_prefix_rep + ".foldchange"
            
            rcfile = fc_prefix_rep + ".readcount"
            
            
            filenames_rep.loc[rep] = [
                                      screen_id,
                                      rcfile,
                                      fc_prefix_rep,
                                      fcfile_rep,
                                      controls 
                                     ]
            
        filenames_screen.loc[screen_id] = [
                                           fcfile,
                                           bffile,
                                           prfile,
                                           screeninfodata[1][i] 
                                          ]


    filenames_rep.to_csv(output_prefix + ".table_replicates.csv")
    filenames_screen.to_csv(output_prefix + ".table_screens.csv")

    screen2rep = filenames_rep.groupby(filenames_rep['screenid']).groups

    for rep in filenames_rep.index:
        fc_prefix = filenames_rep.loc[rep,'fc_prefix']
        ctrl = filenames_rep.loc[rep,'controls']
        replicates = rep
        readcount_screen = filenames_rep.loc[rep,'readcount']
        readcount_data[[readcount_data.dtypes.index[0]] + ctrl.split(",") + replicates.split(",")].to_csv(readcount_screen,sep="\t")
        

        subprocess.run(["echo", 'python' ,BAGEL_PATH, 'fc', '-i', readcount_screen, '-o', fc_prefix, '-c', ctrl])
        subprocess.run(['python', BAGEL_PATH, 'fc', '-i', readcount_screen, '-o', fc_prefix, '-c', ctrl])
        


    # run BAGEL

    for screen_id in filenames_screen.index:
        # load CRISPRcleanR applied foldchange data
        fcdata = pd.DataFrame(index=readcount_data.index,columns=['GENE'])
        firstcheck=True
        reps = list()
        for rep in screen2rep[screen_id]:
            if path.exists(filenames_rep.loc[rep,'fc']): # if exists
                fcfile_rep = pd.read_csv(filenames_rep.loc[rep,'fc'],sep="\t",index_col=0,header=0)
                if firstcheck:
                    fcdata['GENE'] =fcfile_rep['GENE']
                    firstcheck=False

                fcdata[rep] = fcfile_rep[rep].copy()
                reps.append(rep)
            else:
                print(f"\n{filenames_rep.loc[rep,'fc']} doesn't exist. Skipped!\n")

        # write gathered foldchanges into a file.  # drop NA values
        fcdata.dropna(how='any').to_csv(filenames_screen.loc[screen_id,'combinedFC'],sep="\t")

        # run qc for each screens
        if path.exists(QC_PATH):
            qc_command = ['python', QC_PATH,
                          '-i',filenames_screen.loc[screen_id,'combinedFC'],
                          '-e',REF_ESS,
                          '-n',REF_NONESS
                         ]
            try:
                subprocess.run(['echo'] + qc_command)
                subprocess.run(qc_command)
            except:
                print("\n## QC run Failed! Check error messages\n")
                sys.exit(1)
        else:
            print(f"{QC_PATH} doesn't exist. QC skipped!")

        # run BAGEL bf
        bagel_command_bf = ['python', BAGEL_PATH, 'bf', 
                            '-i', filenames_screen.loc[screen_id,'combinedFC'], 
                            '-o', filenames_screen.loc[screen_id,'BF'], 
                            '-e', REF_ESS,
                            '-n', REF_NONESS,
                            '-c', ",".join(reps)
                           ] # add parameters

        # option for multi targeting correction
        if MULTI_TARGETING_CORRECTION==True:
            bagel_command_bf.extend(['-m',  # option for multi-targeting
                                     '--align-info',GUIDE_ALIGNMENT_LIBRARY  # path for guide 
                                    ])

        # option for bootstrapping
        if BOOTSTRAPPING == True:
            bagel_command_bf.extend(['-b'])
            if int(BOOTSTRAPPING_ITERATION) != 1000:
                bagel_command_bf.extend(['-NB',str(int(BOOTSTRAPPING_ITERATION))])

        # option for Normalization
        if int(NORMALIZE_SGRNA_COUNT) > 0:
            bagel_command_bf.extend(['--equalise-sgrna-no',str(int(NORMALIZE_SGRNA_COUNT))])
        if int(NORMALIZE_REPLICATES_N) > 0:
            bagel_command_bf.extend(['--equalise-rep-no',str(int(NORMALIZE_REPLICATES_N))])


        try:
            subprocess.run(['echo'] + bagel_command_bf)
            subprocess.run(bagel_command_bf)
        except:
            print("\n## BAGEL run Failed! Check error messages\n")
            #sys.exit(1)

        # run BAGEL pr
        bagel_command_pr = ['python', BAGEL_PATH, 'pr', 
                            '-i', filenames_screen.loc[screen_id,'BF'], 
                            '-o', filenames_screen.loc[screen_id,'PR'], 
                            '-e', REF_ESS,
                            '-n', REF_NONESS,
                           ]
        if int(NORMALIZE_SGRNA_COUNT) > 0:
            bagel_command_pr.extend(["-k",'NormBF'])

        try:
            subprocess.run(['echo'] + bagel_command_pr)
            subprocess.run(bagel_command_pr)
        except:
            print("\n## BAGEL run Failed! Check error messages\n")
            #sys.exit(1)


    print("\n\n\nJob done")

if __name__ == '__main__':


    #####################
    # initial file check#
    #####################

    if path.exists(BAGEL_PATH) == False:
        print("BAGEL doen't exist. Please check ", BAGEL_PATH)
        sys.exit(1)

    if path.exists(REF_ESS) == False:
        print("The reference essential gene set doen't exist. Please check ", REF_ESS)
        sys.exit(1)

    if path.exists(REF_NONESS) == False:
        print("The reference non-essential gene set doen't exist. Please check ", REF_NONESS)
        sys.exit(1)


    if MULTI_TARGETING_CORRECTION == True:
        if path.exists(GUIDE_ALIGNMENT_LIBRARY) == False:
            print("The library for multi-targeting correction doesn't exist. Please check ", GUIDE_ALIGNMENT_LIBRARY)
            sys.exit(1)


    run_bagel_script()







