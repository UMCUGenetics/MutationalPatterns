#!/bin/bash

# Via bioconda channel
# conda install -c bioconda r-mutationalpatterns
# Author: chapmanb
# https://github.com/bioconda/bioconda-recipes/blob/master/recipes/r-mutationalpatterns/build.sh

# R 3.3.0 not strictly required but earlier versions not tested
# https://github.com/CuppenResearch/MutationalPatterns/issues/3#issuecomment-244312799
# Swap out strict requirements to allow builds and tests
sed -i.bak 's/R (>= 3.3.0)/R (>= 3.2.2)/' DESCRIPTION
sed -i.bak 's/GenomicRanges (>= 1.24.0)/GenomicRanges (>= 1.22.4)/' DESCRIPTION
sed -i.bak 's/IRanges (>= 2.6.0)/IRanges (>= 2.4.8)/' DESCRIPTION
sed -i.bak 's/Biostrings (>= 2.40.0)/Biostrings (>= 2.38.4)/' DESCRIPTION
sed -i.bak 's/GenomeInfoDb (>= 1.8.1)/GenomeInfoDb (>= 1.6.3)/' DESCRIPTION
sed -i.bak 's/SummarizedExperiment (>= 1.2.2)/SummarizedExperiment (>= 1.0.2)/' DESCRIPTION
sed -i.bak 's/VariantAnnotation (>= 1.18.1)/VariantAnnotation (>= 1.16.4)/' DESCRIPTION
sed -i.bak 's/gridExtra (>= 2.2.1)/gridExtra (>= 2.0.0)/' DESCRIPTION
sed -i.bak 's/ggplot2 (>= 2.1.0)/ggplot2 (>= 1.0.1)/' DESCRIPTION
sed -i.bak 's/BiocGenerics (>= 0.18.0)/BiocGenerics (>= 0.16.1)/' DESCRIPTION
#sed -i.bak 's///' DESCRIPTION

echo -e "\n##### DIFF #####\nDESCRIPTION file edited with following differences\n###############\n"
diff DESCRIPTION.bak DESCRIPTION

echo -e "\n##### INFO #####\nBuilding package under $(R --version)\n###############\n"
sleep 2
R CMD INSTALL --build .
EXITSTAT="$?"
echo -e "\n##### STATUS #####\nCompleted building package with exit code ${EXITSTAT}\n###############\n"

## END ##
