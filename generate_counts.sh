
cd /share/workshop/$USER/rnaseq_example
mkdir 04-Counts
mkdir 04-Counts/tmp
for sample in 03-STAR/*_ReadsPerGene.out.tab;
do 
    echo ${sample}
    bname=$(basename $sample)
    fbname=${bname%_ReadsPer*}
    cat ${sample} | tail -n +5 | cut -f4 > 04-Counts/tmp/${fbname}.count
done


for sample in 04-Counts/tmp/*_SE.count;
do 
    file1=$sample
    file2=$(echo ${file1} | sed -r 's/_SE./_PE./g')
    output_file=$(echo ${file1} | sed -r 's/_SE./_total./g')

    # Use paste to combine lines from both files and awk to sum them
    paste "$file1" "$file2" | awk '{print $1 + $2}' > "$output_file"
done

tail -n +5 03-STAR/WM_Mixed-10-ng-1_htstream_SE_ReadsPerGene.out.tab | cut -f1 > 04-Counts/tmp/geneids.txt

paste 04-Counts/tmp/geneids.txt 04-Counts/tmp/*_total.count > 04-Counts/tmp/tmp.out


# Extract sample names and combine files
for file in 04-Counts/tmp/*_htstream_total.count; do
    # Extract the sample name (remove '_htstream_total.count' from the filename)
    bname=$(basename $(basename $file))
    sample_name="${bname%%_htstream_total.count}"
    echo -n "$sample_name " # Output the sample name followed by a space
done > 04-Counts/tmp/sample_names.txt

samples='ls 04-Counts/tmp/*_total.count'
cat <(cat 04-Counts/tmp/sample_names.txt | paste -s) 04-Counts/tmp/tmp.out | less > 04-Counts/rnaseq_workshop_counts.txt
rm -rf 04-Counts/tmp
head 04-Counts/rnaseq_workshop_counts.txt
