# Original: https://divingintogeneticsandgenomics.rbind.io/post/split-a-10xscatac-bam-file-by-cluster/

# time bamtools filter -tag CB:Z:AAACTGCAGAGCAGCT-1 -in atac_v1_pbmc_5k_possorted_bam.bam -out AAACTGCAGAGCAGCT-1.bam

# time samtools view atac_v1_pbmc_5k_possorted_bam.bam | awk -v tag="CB:Z:AAACTGCAGAGCAGCT-1" 'index($0,tag)>0' >> AAACTGCAGAGCAGCT-1.sam


import pybam
import csv

cluster_dict = {}
with open('clusters.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    #skip header
    header = next(csv_reader)
    for row in csv_reader:
        cluster_dict[row[0]] = row[1]

clusters = set(x for x in cluster_dict.values())


# open the number of bam files as the same number of clusters, and map the out file handler to the cluster id

header = pybam.read(input_bam).file_header
fouts_dict = {}
for cluster in clusters:
    fout = open("cluster" + cluster + ".sam", "w")
    fout.write(header)
    fouts_dict[cluster] = fout

for read in pybam.read(input_bam):
        ## not always the same position in the list for the CB tag
        ## there could be no CB tag for a certian read as well
        ## it will return empty list
        CB_list = [ x for x in read.sam_tags_list if x[0] == "CB"]
        if CB_list:
            cell_barcode = CB_list[0][2]
            cluster_id = cluster_dict.get(cell_barcode)
            if cluster_id:
                fouts_dict[cluster_id].write(read.sam + '\n')
        
## do not forget to close the files
for fout in fouts_dict.values():
    fout.close()