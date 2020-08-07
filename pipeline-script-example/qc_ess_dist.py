#!/usr/bin/python

# a script for run BAGEL and CRISPRcleanR
# usage) python3 run_bagel_crisprcleanr.py -i [readcount table] -o [output file tag] -c [control columns] -x [experiment columns]

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import numpy as np
import pandas as pd
import seaborn as sns
import math
sns.set_context('talk')
from os import path
import click


def cohensd_func(a,b):
    s=math.sqrt(((len(a)-1)*np.var(a) + (len(b)-1)*np.var(b)) /
           (len(a)+len(b)-2)
          )
    return (np.mean(b) - np.mean(a)) / s

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input-foldchange', required=True, type=click.Path(exists=True))
@click.option('-e', '--essential-genes', required=True, type=click.Path(exists=True))
@click.option('-n', '--non-essential-genes', required=True, type=click.Path(exists=True))
@click.option('-o', '--output-prefix')
def qc_ess_distribution(
        input_foldchange, essential_genes, non_essential_genes, output_prefix
):
    print("Fold-change file  :", input_foldchange)
    print("Reference essential genes  :", essential_genes)
    print("Reference non-essential genes  :", non_essential_genes)
    if output_prefix == None:
        output_prefix = input_foldchange
    print("Output prefix :", output_prefix)

    print("\n")

    fcdata = pd.read_csv(input_foldchange,index_col=0,header=0,sep="\t")
    ess = pd.read_csv(essential_genes, index_col=0,header=0,sep="\t").index  # load ref_core-essential
    non = pd.read_csv(non_essential_genes, index_col=0,header=0,sep="\t").index  # load ref_non-essential
    

    for i in range(1, len(fcdata.dtypes)):

        sample = fcdata.columns.values[i]
       
        fc_ess = fcdata[sample][fcdata[fcdata.columns[0]].isin(ess)]
        fc_noness = fcdata[sample][fcdata[fcdata.columns[0]].isin(non)]
        fig,ax = plt.subplots(1,figsize=(7,5))


        ax = sns.kdeplot(fc_ess,color='red',ax=ax)
        ax = sns.kdeplot(fc_noness,color='blue',ax=ax)
        ax.legend(['Ess','Non-ess'],loc='best')
        cohensd = cohensd_func(fc_ess.values,fc_noness.values)

        print (f"Quality Score ({sample}) = {cohensd:.4g}")  # -2.0 < QS < 2.0, QS > 1.0 is a good sample.
        
        
        ax.set_xlabel('Fold Change')
        ax.set_ylabel('Density')
        ax.set_title(f"{sample}\nQuality Score = {cohensd:.4g}")

        plt.subplots_adjust(top=0.85,bottom=0.2)
        fig.savefig(f"{output_prefix}.qc.{sample}.pdf",format='pdf')
        print(f"Figure saved -> {output_prefix}.qc.{sample}.pdf")
    
    
if __name__ == '__main__':
    qc_ess_distribution()
