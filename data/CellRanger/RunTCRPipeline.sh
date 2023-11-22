# Run the TCR pipeilne for each sample

# Loop over each folder
for folder in *_*
do

echo $folder

# Change into that folder
cd $folder

# Run the pipeline
python3 ../CreateTCRData.py \
clonotypes.csv \
filtered_contig_annotations.csv \
filtered_feature_bc_matrix/barcodes.tsv.gz \
${folder}_TCR.csv

# Change out of that folder
cd ..

done
