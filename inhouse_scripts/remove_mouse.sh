# 20210403 created by pjh
# removes mouse aligned reads

input_bam=$1
output_bam=$2

samtools view -h $input_bam $(seq 1 22) X Y MT |
    awk '
    BEGIN{FS="\t" ; OFS="\t"}
    $0 ~ /^@/{print}
    $0 !~ /^@/{
        if ( $7 !~ /^Mus/ ) {print}
    }' |
    samtools view -hb -o $output_bam

samtools index $output_bam