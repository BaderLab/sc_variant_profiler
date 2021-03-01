import pysam
import os
unsplit_file = "CB_containing_reads.bamsorted.bam"  ## reading a file with 400milion reads (with 6milion name)
barcodesfile = open("barcodes.tsv", "r+")          ## reading a file with my interested names (11,000)
barcodes = [line.strip() for line in barcodesfile]
out_dir = "Cells_bamFiles"
os.mkdir(out_dir)
#variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile(unsplit_file, "rb")
for read in samfile.fetch(until_eof=True):          ## going over every line of 400milion reads
    # barcode itr for current read
    CB_itr = read.get_tag('CB')
    if len(barcodes) == 0:
        pass
    if CB_itr not in barcodes:                     ## if name of the read is not in my interested names continiue                  
        continue
    else:                                          ##esle make a file and store the reads for evey 11,000 names
        #if change in barcode or first line; open new file
        if (CB_itr != CB_hold or itr == 0):
            # close previous split file, only if not first read in file
            if (itr != 0):
                split_file.close()
            CB_hold = CB_itr
            itr += 1
            print(itr)
            split_file = pysam.AlignmentFile(out_dir + "/CB_{}.bam".format(itr), "wb", template=samfile)

        # write read with same barcode to file
        split_file.write(read)
split_file.close()
samfile.close()
