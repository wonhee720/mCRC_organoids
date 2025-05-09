#!/bin/bash
# v2
	#201106
	#210303
		# changed main program to v2.2.3
		# changed splitter version to v2 (since some split results in bcf.gz format are not indexable)
	#210304
		# changed splitter version to v3 (since some split results in bcf.gz format are not indexable)
	#210306
		# edited log printing system
		# split input files and split output files are stored in separate directories
		# uses mktemp for creating temporary directories
		# temporary directories are created in the same directory as the final output file
	#210307
		# gnu parallel is used in local scheduler mode
# v3
	#210307
		# split unit modification is removed (too small split unit results in low efficiency)
		# parallelization algorithm is modified
# v4
	#210309
		# parallelization algorithm is modified
	#210312
		# does not remove error output files
		# uses main program v2.2.4
	#210315
		# Presence of index file for each input bam is checked at the beginning

# v5
	#210322
		# added --marginOverlap option
	#210402
		# main program version is v2.2.5
	#210428
		# made a logging function
	#210526
		# main program version is v2.2.6
		# conda executables are used

# v6
	#210602
		#init
		#uses main program v2.2.7
		# -w option is not used any more
		# -f option can accept an arbitrary fasta file path
	
# v7
	#210615
		#init
		#uses main program v2.2.8
		#temporary directory is removed by default

# v8
	#210708
		#init
		#adds contig header to input copy if there is not one

# v9
	# 210804
		# init
		# input vcf is divided into  min(max_split_number, min(nt*nt_factor, lineno)) files, instead of default width of 100 lines
		# job array is not used
		# nt_factor is 1

# v10
	# 211115
		# init
		# modified for use with slurm
		# --pbs-intv option removed
		# -N default value changed

# v10-lwh
	# 231206 : Changed after server breakdown by LWH


trap 'echo error $LINENO' ERR
set -eu


printErr()
{
	echo "readinfoAnnot_wrapper : $*" >> /dev/stderr
}


usage()
{
	(
		echo "Usage: $(basename $0) \\"
		echo $'\t''-i : input variant file. May be any of vcf, vcf.gz, bcf, or bcf.gz.'
		echo $'\t''-b : Comma-separated list of bam files to get information from. Must be in the same order with -s option value.'
		echo $'\t''-o : Path of result output file.'
		echo $'\t''-s : Comma-separated list of sample name. A new sample column will be created if not already present.'
		echo $'\t''-f : Reference fasta file. May be "hg19" or "hg38", which will lead to using a preset fasta file.'
		echo $'\t\t''"19" or "38" are treated as "hg19" or "hg38".'
		echo $'\t\t''Trinucleotide context is extracted from this file.'
		echo $'\t\t''If the input vcf file does not contain contig headers, the fai file is used to create contig header lines.'
		echo

		echo $'\t''[optional] -O : Format of output file. One of v,z,u,b(identical to bcftools). Default: z'
		echo $'\t''[optional] --sched : parallelization option. One of "local" or "slurm". Default: local'

		echo $'\t''[optional] -p : Number of subprocesses for parallelization. Default: 1'
		echo $'\t\t''(Input vcf split number: min( max_split_number(default 500), min(p*factor(default 2), lineno) ))'

		echo $'\t''[optional] -N : Job name prefix when slurm is used. Default: readinfoAnnot_<input vcf basename>'
		#echo $'\t''[optional] --pbs-intv : PBS job monitoring interval. Default: 600 seconds'
		echo $'\t''[optional] --dont-rm-tmpdir : If set, temporary directory is not removed. (Removed by default)'

		echo $'\t''[optional] --satelliteVer : Reference genome version for satellite information. Valid values are "hg19" or "hg38".'
		echo $'\t\t''If --satelliteVer option is not set and -f value is "19" or "hg19", --satelliteVer value becomes "hg19". Similar for -f "38" or "hg38".'
		echo $'\t\t''If -f is neither "19", "hg19", "38", "hg38" and --satelliteVer is not set, satellite region skipping is not performed.'
		echo $'\t''[optional] --level : level of annotation. Possible values are 2 (full annotation) or 1 (read count and vaf only). Default: 2'
		echo $'\t''[optional] --marginOverlap : Must be T or F. If T, marginOverlap function is applied. Default: F'
		echo $'\t''[optional] --limBQ : (SNV only) reads with BQ of target base below this value will be treated as lowqual reads. Default: 20'
		echo $'\t''[optional] --limMQ : reads with MQ below this value will be treated as lowqual reads. Default: 20'
		echo $'\t''[optional] --MMfrac : reads with mismatch (excluding the variant itself) more than'
		echo $'\t\t''<read length including soft-clipped bases>*<MMfrac> will be treated as lowqual reads. Default: 0.1'
		echo $'\t''[optional] --flankLen : <flankLen> bases on either left or right of the variant position must be cigar M to be considered as not an irrelevant read. Default: 1'

		echo $'\t''[optional] --limBQ : reads with BQ of target base below this value will be treated as non-supporting reads. Default: 20'
		echo $'\t''[optional] --limMQ : reads with MQ below this value will be treated as non-supporting reads. Default: 20'
		echo $'\t''[optional] --MMfrac : reads with mismatch (excluding the variant itself) more than'
		echo $'\t\t''<read length including soft-clipped bases>*<MMfrac> will be treated as non-supporting reads. Default: 0.1'
		echo $'\t''[optional] --flankLen : <flankLen> bases on either left or right of the variant position must be cigar M to be considered as ref/alt supporting. Default: 1'
	) > /dev/stderr
	exit 1
}


# parse arguments
parse_arguments()
{
	while [[ $# -gt 0 ]] ; do
		case "$1" in
			# wrapper options - required
			-i)
				shift ; vcf_original="$1" ; shift ;;
			-b)
				shift ; bam_list_str="$1" ; shift ;;
			-o)
				shift ; outfile="$1" ; shift ;;
			-s)
				shift ; sample_list_str="$1" ; shift ;;

			# main script options - required
			-f)
				shift ; main_opts_f="$1" ; shift ;;

			# wrapper options - not required
			-O)
				shift ; outfile_format="$1" ; shift ;;
			--sched)
				shift ; sched="$1" ; shift ;;
			-p)
				shift ; nt="$1" ; shift ;;
			-N)
				shift ; jobname_pf="$1" ; shift ;;
			#--pbs-intv)
				#shift ; pbs_intv="$1" ; shift ;;
			--dont-rm-tmpdir)
				rmtmp=0 ; shift ;;

			# main script options - not required
			--satelliteVer)
				shift ; main_opts_satelliteVer="$1" ; shift ;;
			--level)
				shift ; main_opts_level="$1" ; shift ;;
			--marginOverlap)
				shift ; main_opts_marginOverlap="$1" ; shift ;;
			--limBQ)
				shift ; main_opts_limBQ="$1" ; shift ;;
			--limMQ)
				shift ; main_opts_limMQ="$1" ; shift ;;
			--MMfrac)
				shift ; main_opts_MMfrac="$1" ; shift ;;
			--flankLen)
				shift ; main_opts_flankLen="$1" ; shift ;;

			# help and unknown
			-h|--help)
				usage ;;
			*)
				echo -e "Unknown option: ${1}\n" > /dev/stderr ; usage ;;
		esac
	done
}


sanity_check()
{
	# check required parameters
	required_params=(
		vcf_original
		outfile
		bam_list_str
		sample_list_str
		main_opts_f
	)
	for param in ${required_params[@]} ; do
		if [[ -z ${!param:-} ]] ; then usage ; fi
	done


	# -p option
	if [[ -n ${nt:-} ]] ; then
		if [[ $nt -lt 1 || $nt -gt 400 ]] ; then
			printErr "-p must 1 <= p <= 400"
			exit 1
		fi
	fi

	# --satelliteVer option
	if [[ -n ${main_opts_satelliteVer:-} ]] ; then
		case $main_opts_satelliteVer in
			hg19|hg38) : ;;
			*) printErr '--satelliteVer option must be "hg19" or "hg38"' ; exit 1 ;;
		esac
	fi

	# --level option
	if [[ -n ${main_opts_level:-} ]] ; then
		case $main_opts_level in
			1|2) : ;;
			*) printErr '--level option must be "1" or "2"' ; exit 1 ;;
		esac
	fi

	# --sched option
	if [[ -n ${sched:-} ]] ; then
		case $sched in
			local|slurm) : ;;
			*) printErr '--sched option must be "local" or "slurm".' ; exit 1 ;;
		esac
	fi

	# check output directory permission
	outdir_perm=$(stat -c %a $(dirname $outfile))
	if ! [[ $(( ( 8#${outdir_perm} & 8#600 ) == 8#600 )) = 1 ]] ; then
		printErr "You do not have read and write permission for the directory of the specified output file."
		exit 1
	fi

	# check bam file and sampleID numbers
	len1=$( echo $bam_list_str | tr ',' ' ' | wc -w )
	len2=$( echo $sample_list_str | tr ',' ' ' | wc -w )
	if ! [[ $len1 = $len2 ]] ; then
		printErr "Numbers of input bam files and sample names are different."
		exit 1
	fi

	# check if bam index files exist
	bam_list=( $( echo $bam_list_str | tr ',' ' ' ) )
	for bam in ${bam_list[@]} ; do
		if ! [[ -e ${bam}.bai || -e ${bam/%bam/bai} ]] ; then
			printErr "file name: $bam"
			printErr "Index file of this bam was not found."
			exit 1
		fi
	done
}


set_default_arguments()
{
	# default argument values
	if [[ -z ${outfile_format:-} ]] ; then outfile_format=z ; fi
	if [[ -z ${sched:-} ]] ; then sched=local ; fi
	if [[ -z ${nt:-} ]] ; then nt=1 ; fi
	if [[ -z ${jobname_pf:-} ]] ; then jobname_pf="readinfoAnnot_$(basename $vcf_original)" ; fi
	#if [[ -z ${pbs_intv:-} ]] ; then pbs_intv=600 ; fi
	if [[ -z ${rmtmp:-} ]] ; then rmtmp=1 ; fi
}


edit_raw_arguments()
{
	# absolute path
	if ! [[ $vcf_original =~ ^/ ]] ; then vcf_original=${PWD}/${vcf_original} ; fi
	if ! [[ $outfile =~ ^/ ]] ; then outfile=${PWD}/${outfile} ; fi

	# get bam file list and sampleID list
	bam_list=( $( echo $bam_list_str | tr ',' ' ' ) )
	sample_list=( $( echo $sample_list_str | tr ',' ' ' ) )

	# set --satelliteVer value according to -f value
	if [[ -z ${main_opts_satelliteVer:-} ]] ; then
		if [[ $main_opts_f = hg19 || $main_opts_f = 19 ]] ; then
			main_opts_satelliteVer=hg19
		elif [[ $main_opts_f = hg38 || $main_opts_f = 38 ]] ; then
			main_opts_satelliteVer=hg38
		fi
	fi
}


make_main_opts_string()
{
	main_opts_string=$(
		for varname in ${!main_opts_*} ; do
			key=${varname#main_opts_}
			val=${!varname}
			if [[ $key =~ ^.$ ]] ; then
				printf %s "-${key} ${val} "
			else
				printf %s "--${key} ${val} "
			fi
		done
	)
}


set_constants()
{
#	main_program=/home/users/pjh/scripts/annotation/short_variants/readinfoAnnot.v2.2.8.py
	main_program=/home/users/wonhee720/tools/scripts/pjh_scripts/short_variants/readinfoAnnot.v2.2.8-lwh.py # LWH
#	splitter=/home/users/pjh/scripts/vcf_editing/divide_vcf.v4.py
	splitter=/home/users/wonhee720/tools/scripts/pjh_scripts/vcf_editing/divide_vcf.v2-lwh.py # LWH
	bcftools=/home/users/pjh/conda_bin/bcftools
	parallel=/home/users/pjh/conda_bin/parallel
	bash=/usr/bin/bash
	#PBS_helper=/home/users/pjh/scripts/others/PBS_submit_helper.v4.sh
	sbatch_helper=/home/users/pjh/scripts/slurm_utils/sbatch-helper

	submit_interval=1
	tmpfile_format=z
	job_ncore=1

	nt_factor=1
	max_split_number=500
}




# main functions
make_tmpdir()
{
	printErr "Making temporary directories"

	tmpdir_base=$( mktemp -d -p $(dirname $outfile) readinfoAnnotTmpdir.$(basename $vcf_original).XXXXXXXX )

	tmpdir_scripts=${tmpdir_base}/scripts ; mkdir $tmpdir_scripts
	tmpdir_log=${tmpdir_base}/logs ; mkdir $tmpdir_log

	# tmpdir_step_list has length (number of bam files + 1)
	declare -ga tmpdir_step_list
	for idx in $(seq -w 0 $len1) ; do
		tmpdir_step=${tmpdir_base}/step${idx}
		mkdir $tmpdir_step
		tmpdir_step_list+=( $tmpdir_step )
	done

	tmpdir_input_list=( ${tmpdir_step_list[@]:0:$((${#tmpdir_step_list[@]}-1))} )
	tmpdir_output_list=( ${tmpdir_step_list[@]:1} )
}


make_input_copy()
{
	# check if the original file has contig headers
	if [[ $( $bcftools view -h $vcf_original | grep -c '^##contig' ) = 0 ]] ; then
		input_has_contig_hdr=0
	else
		input_has_contig_hdr=1
	fi

	# get reference fai file path
	if [[ $main_opts_f = hg19 || $main_opts_f = 19 ]] ; then
		fai=/home/users/pjh/References/reference_genome/GRCh37/hg1kv37/human_g1k_v37.fasta.fai
	elif [[ $main_opts_f = hg38 || $main_opts_f = 38 ]] ; then
		fai=/home/users/pjh/References/reference_genome/GRCh38/GCA_for_alignment_pipelines/full_plus_hs38d1_210607/GRCh38_full_decoy_HLA.fasta.fai
	else
		fai=${main_opts_f}.fai
	fi

	# create input copy
	input_copy=${tmpdir_base}/input_copy.vcf.gz
	if [[ $input_has_contig_hdr = 0 ]] ; then
		$bcftools reheader -f $fai $vcf_original | $bcftools view -O z -o $input_copy
	else
		$bcftools view -O z -o $input_copy $vcf_original
	fi
}


get_split_number()
{
	# split number = min(max_split_number, min(nt*nt_factor, lineno))

	split_number=$(( $nt * $nt_factor ))

	lineno=$($bcftools view -H $vcf_original | wc -l)
	if [[ $lineno -lt $split_number ]] ; then
		split_number=$lineno
	fi

	# Set maximum split number as 500
	if [[ $split_number -gt $max_split_number ]] ; then
		split_number=$max_split_number
	fi
}


split_input()
{
	printErr "Splitting input file"
#	$splitter -N $split_number -O $tmpfile_format -p input $input_copy ${tmpdir_input_list[0]}  # split files are in vcf.gz format and not indexed
	$splitter -N $split_number -O $tmpfile_format $input_copy ${tmpdir_input_list[0]}  # LWH

}


make_job_scripts()
{
	printErr "Making job scripts"

	split_input_list=( $( find ${tmpdir_input_list[0]} -type f ! -name '*.csi' | sort ) )

	for idx_pad in $( seq -w 0 $(( ${#split_input_list[@]} - 1 )) ) ; do
		idx=$( echo $idx_pad | awk '{print gensub(/^0*([^0].*)$/, "\\1", 1, $0)}' )
		jobname=${jobname_pf}_${idx_pad}
		split_input=${split_input_list[$idx]}
		bname=$(basename $split_input)
		tmpfile_input_list=( $( for tmpdir in ${tmpdir_input_list[@]} ; do echo ${tmpdir}/${bname} ; done ) ) 
		tmpfile_output_list=( $( for tmpdir in ${tmpdir_output_list[@]} ; do echo ${tmpdir}/${bname} ; done ) )

		script=${tmpdir_scripts}/script.${idx_pad}.sh

		# Added option to exclude specific node
		cat > $script <<-EOF
		#!/bin/bash
		#SBATCH -N 1 -n 1 -c $job_ncore --exclude=bnode7
		#SBATCH -J $jobname
		#SBATCH -o /dev/null
		set -eu
		
		tmpfile_input_list=( ${tmpfile_input_list[@]} )
		tmpfile_output_list=( ${tmpfile_output_list[@]} )
		sample_list=( ${sample_list[@]} )
		bam_list=( ${bam_list[@]} )
		log_prefix=${tmpdir_log}/${bname}
		main_program="${main_program}"
		main_opts_string="${main_opts_string}"
		tmpfile_format="${tmpfile_format}"
		
		for idx_pad in \$( seq -w 0 \$(( \${#tmpfile_input_list[@]} - 1 )) ) ; do
		    idx=\$( echo \$idx_pad | awk '{print gensub(/^0*([^0].*)$/, "\\\\1", 1, \$0)}' )

		    tmpfile_input=\${tmpfile_input_list[\$idx]}
		    bam_input=\${bam_list[\$idx]}
		    tmpfile_output=\${tmpfile_output_list[\$idx]}
		    sampleID=\${sample_list[\$idx]}
		    log=\${log_prefix}.step\${idx_pad}.log

		    echo \$(date) -- Beginning annotation for \$bam_input >> \$log
		    \$main_program \\
		        \$main_opts_string \\
		        -O \$tmpfile_format \\
		        -o \$tmpfile_output \\
		        -s \$sampleID \\
		        -i \$tmpfile_input \\
		        -b \$bam_input &>> \$log
		done
		EOF
	done
}


run()
{
	printErr "Running annotation jobs"

	case $sched in
		local )
			find ${tmpdir_scripts} -type f | sort | $parallel -j $nt $bash
			;;
		slurm)
			$sbatch_helper --dir $tmpdir_scripts --submit-interval $submit_interval
			;;
	esac
}


concat()
{
	printErr "Concatenating split result files"

	$bcftools concat -O $outfile_format -o $outfile $( find ${tmpdir_output_list[$(( ${#tmpdir_output_list[@]} - 1 ))]} -type f ! -name "*.csi" | sort )
	if [[ $outfile_format =~ ^(z|b)$ ]] ; then $bcftools index $outfile ; fi
}


rm_tmpdirs()
{
	if [[ $rmtmp = 1  ]] ; then
		rm -r $tmpdir_base
	fi
}


# main
parse_arguments "$@"
sanity_check
set_default_arguments
edit_raw_arguments
make_main_opts_string
set_constants

make_tmpdir
make_input_copy
get_split_number
split_input
make_job_scripts
run
concat
rm_tmpdirs

printErr "All finished"
