# BAGEL2

BAGEL2 software

The Bayesian Analysis of Gene Essentiality 2 (BAGEL2) is a robust gene essentiality identifier from CRISPR-cas9 genome wide pooled library screen data.

## Download BAGEL2

```
git clone https://github.com/hart-lab/bagel
```

## Test BAGEL2-CRISPRcleanR Pipeline

```
cd pipeline-script-example
```

```
python run_bagel_crisprcleanr.py -i HAP1-TKOv3-EXAMPLE.txt -s HAP1-TKOv3-EXAMPLE-SCREENINFO.txt -o HAP1-TKOv3-bagel
```

## References
1. Kim, E., Hart, T. Improved analysis of CRISPR fitness screens and reduced off-target effects with the BAGEL2 gene essentiality classifier. *Genome Med* 13, 2 (2021). https://doi.org/10.1186/s13073-020-00809-3
2. Hart, T., Brown, K. R., Sircoulomb, F., Rottapel, R., & Moffat, J. Measuring error rates in genomic perturbation screens: gold standards for human functional genomics. *Molecular systems biology*, 10(7), 733. (2014).

## License

[MIT License](LICENSE)





