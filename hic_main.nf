#!/usr/bin/env nextflow

/*
 * Hi-C Workflow Summary:
 *
 *		Pre-processing Steps
 *		Two-Steps Read Alignment
 *		Combining reads
 *		Filtering
 *		Build Matrix
 *		Normalise Matrix
 *		MultiQC Report
 */

/*
 * Input parameters
 *
 *
 * --reads			Path to input sequence data in fastq format
 * --genome			Path to reference genome
 * --outdir			Path to directory for outputs
 * --restricion_site		Specify restriction enzyme cutting motif
 * --ligation_site		Specify ligation motif for trimming step
 * --bwt_index			Path to bowtie2 index including index base name
 * --bin_size			Specify the resolution for generating contact maps
 * --rm_singleton		Remove unpaired ligtion products
 * --rm_multi			Remove reads which were multi mapped
 * --min_mapq			Set minimum mapping quality
 */

/*
 * Set default input parameters
 *
 *
 * params.reads = "/home/conor/Documents/Git_Repositories/Hi-C/reads/*{1,2}.fastq"
 * params.genome = "/home/conor/Documents/Git_Repositories/Hi-C/reference/GRCh38.primary_assembly.genome.fa"
 * params.outdir = "/home/conor/Documents/Git_Repositories/Hi-C/results"
 * params.restriction_site = "hindiii"
 * params.bwt_index = "/home/conor/Documents/Git_Repositories/Hi-C/bwt_index/genome.index"
 * params.ligation_site = "AAGCTAGCTT"
 * params.bin_size = '1000000,500000'
 * params.rm_singleton = true
 * params.rm_multi = true
 * params.min_mapq = 10
 */
 
/*
 * Set Up Channnels
 */

	// Set up input reads channel
	raw_reads = Channel.create()
	raw_reads_2 = Channel.create()
	
	Channel
		.fromFilePairs( params.reads )
		.separate( raw_reads, raw_reads_2 ) { a -> [tuple(a[0], a[1][0]), tuple(a[0], a[1][1])] }
		
		raw_reads = raw_reads.concat( raw_reads_2 ).dump(tag: "data")

	// Set up bowtie2 index channel
	lastPath = params.bwt_index.lastIndexOf(File.separator)
	bwt_dir = params.bwt_index.substring(0,lastPath+1)
	bwt_base = params.bwt_index.substring(lastPath)

	Channel.fromPath( bwt_dir )
		.into { bwt_index_first_step; bwt_index_second_step }
		
	// Set up mapping resolution for building contact map channel
	map_res = Channel.from( params.bin_size.tokenize(',') )

/************************************************************************** 
 * Pre-Processing Steps 
 */ 

// Extract chromosome sizes from reference genome 

process get_chrom_size { 
	publishDir "${params.outdir}/reference" 

	input: 
	path fasta from params.genome // reference genome in fa format 

	output: 
	file "*.size" into chromosome_size, chromosome_size_cool 

	script: 
	""" 
	samtools faidx ${fasta} 
	cut -f1,2 ${fasta}.fai > chrom.size 
	""" 
}

// Extract restriction fragment coordinates from reference into BED file for interaction mapping 

process get_restriction_fragments { 
	tag "$fasta ${params.restriction_site}" 
	publishDir "${params.outdir}/reference" 

	input: 
	path fasta from params.genome
	
	output: 
	file "*.bed" into res_frag_file 
	
	script: 
	""" 
	/HiC-Pro-devel/bin/utils/digest_genome.py -r ${params.restriction_site} -o restriction_fragments.bed ${fasta} 
	""" 
} 

/**************************************************************************** 
 * 
 * Main Workflow
 * 
 */ 

/* 
 * Step 1: Two-Steps Read Mapping 
 */ 

/* 
 * First end-to-end alignment on all reads 
 * This will output unmapped reads to separate file as these are assumed to span ligation sites 
 * These reads will be rescued by trimming the 3' end of the reads to the ligation site 
 * 
 * Note: the bowtie2 index file can be generated using the build_index.nf script in the current folder 
 * This was kept separate in the interest of reducing the computation time of the pipeline 
 */ 

/* 
 * Command Flags: 
 * 
 *		--rg-id				set the read group ID 
 *		--rg				add sample name to read group header 
 *		--very-sensitive	preset sensitivity/accuracy/speed parameters 
 *		-L					length of the seed substrings to align 
 *		--score-min		set threshold for alignment score 
 *		--end-to-end		perform end-to-end alignment (alignment involving all read characters) 
 *		--reorder			output file in order of reads in original file 
 *		--un				output unpaired reads which fail to align 
 *		-U					output files containing unpaired reads 
 */ 

 process end_to_end_alignment_1 { 
	tag "$prefix" 
	publishDir "${params.outdir}/mapping", mode: 'copy'
	
	input: 
	set val(sample), file(reads) from raw_reads 
	file index from bwt_index_first_step.collect() 

	output: 
	set val(prefix), file("${prefix}_unmap.fastq") into unmapped_end_to_end // unmapped reads spanning ligation site 
	set val(prefix), file("${prefix}.bam") into end_to_end_bam // first set of mapped reads 

	script: 
	// extract sample name and read direction (paired) by parsing out extension 
	prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/ 
	
	""" 
	bowtie2 --rg-id BMG \
			--rg SM:${prefix} \
			-x ${index}/${bwt_base} \
			--very-sensitive \
			-L 30 \
			--score-min L,-0.6,-0.2 \
			--end-to-end \
			--reorder \
			--un ${prefix}_unmap.fastq \
			-U ${reads} | samtools view -F 4 -bS - > ${prefix}.bam 
	""" 
}

// Trim unmapped reads at ligation site 

process trim_unmapped_reads { 
	tag "$prefix" 
	publishDir "${params.outdir}/mapping", mode: 'copy' 

	input: 
	set val(prefix), file(reads) from unmapped_end_to_end 

	output: 
	set val(prefix), file("${prefix}_trimmed.fastq") into trimmed_reads // trimmed unmapped reads 

	script: 
	// Use cutsite trimming tool to trim reads at specified ligation site 
	""" 
	/HiC-Pro-devel/scripts/cutsite_trimming --fastq $reads --cutsite ${params.ligation_site} --out ${prefix}_trimmed.fastq 
	""" 
} 

// Align trimmed unmapped reads 

process end_to_end_alignment_2 { 
	tag "$prefix" 
	publishDir "${params.outdir}/mapping", mode: 'copy' 

	input: 
	set val(prefix), file(trimmed_reads) from trimmed_reads 
	file index from bwt_index_second_step.collect() 

	output: 
	set val(prefix), file("${prefix}_trimmed.bam") into trimmed_bam // aligned reads which were initially unmapped 

	script: 
	prefix = trimmed_reads.toString() - ~/(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/ 
	""" 
	bowtie2 --rg-id BMG \
			--rg SM:${prefix} \
			-x ${index}/${bwt_base} \
			--very-sensitive \
			-L 20 \
			--score-min L,-0.6,-0.2 \
			--end-to-end \
			--reorder \
			-U ${trimmed_reads} | samtools view -bS - > ${prefix}_trimmed.bam 
	""" 
} 

// Merge aligned reads from both alignment steps into single paired-end BAM 
// (i.e. incorporate the rescued reads which span ligation site into BAM) 

process merge_alignment_steps { 
	tag "$sample = $bam1 + $bam2" 
	publishDir "${params.outdir}/mapping", mode: 'copy' 
	
	input: 
	set val(prefix), file(bam1), file(bam2) from end_to_end_bam.join( trimmed_bam ) 

	output: 
	set val(sample), file("${prefix}_bwt2merged.bam") into bwt2_merged_bam // merged aligned reads 
	set val(oname), file("${prefix}.mapstat") into all_mapstat // mapping statistics 

	script: 
	// Extract sample names 
	sample = prefix.toString() - ~/(_R1$|_R2$|_val_1$|_val_2$|_1$|_2$)/ 
	tag = prefix.toString() =~/_R1$|_val_1$|_1$/ ? "R1" : "R2" 
	oname = prefix.toString() - ~/(\.[0-9]+)$/ 
	""" 
	samtools merge -@ ${task.cpus} -f ${prefix}_bwt2merged.bam ${bam1} ${bam2}
	samtools sort -@ ${task.cpus} -m 800M -n -T /tmp/ -o ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam 

	mv ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam 

	## Create mapping statistics files 
	echo "## ${prefix}" > ${prefix}.mapstat 
	echo -n "total_${tag}\t" >> ${prefix}.mapstat 
	samtools view -c ${prefix}_bwt2merged.bam >> ${prefix}.mapstat 
	echo -n "mapped_${tag}\t" >> ${prefix}.mapstat 
	samtools view -c -F 4 ${prefix}_bwt2merged.bam >> ${prefix}.mapstat 
	echo -n "global_${tag}\t" >> ${prefix}.mapstat 
	samtools view -c -F 4 ${bam1} >> ${prefix}.mapstat 
	echo -n "local_${tag}\t" >> ${prefix}.mapstat 
	samtools view -c -F 4 ${bam2} >> ${prefix}.mapstat 
	""" 
} 

/* 
 * Step 2: Combining reads 
 */ 

// Combine alignment files so that crosslinked reads can be paired for contact mapping 

process combine_alignment_files{ 
	tag "$sample = $r1_prefix + $r2_prefix" 
	publishDir "${params.outdir}/mapping", mode: 'copy'
	
	input: 
	set val(sample), file(aligned_bam) from bwt2_merged_bam.groupTuple() // Merged BAM files for each read 

	output: 
	set val(sample), file("${sample}_bwt2pairs.bam") into paired_bam // merged two read files into single BAM 
	set val(oname), file("*.pairstat") into all_pairstat // pairing statistics 

	script: 
	r1_bam = aligned_bam[0] 
	r1_prefix = r1_bam.toString() - ~/_bwt2merged.bam$/ 
	r2_bam = aligned_bam[1] 
	r2_prefix = r2_bam.toString() - ~/_bwt2merged.bam$/
	oname = sample.toString() - ~/(\.[0-9]+)$/

	def opts = "-t" 
	opts = params.rm_singleton ? "${opts}" : "--single ${opts}" // Specify remove singleton reads 
	opts = params.rm_multi ? "${opts}" : "--multi ${opts}" // Specify remove multimapped reads 
	if ("$params.min_mapq".isInteger())	opts="${opts} -q ${params.min_mapq}" // Specify minimum mapping quality 
	""" 
	/HiC-Pro-devel/scripts/mergeSAM.py -f ${r1_bam} -r ${r2_bam} -o ${sample}_bwt2pairs.bam ${opts} 
	""" 
} 

/* 
 * Step 3: Filtering 
 */

// Remove invalid reads such as unligated and self ligated products 

process get_valid_pairs {
	tag "$sample"
	publishDir "${params.outdir}/hic_results/data", mode: 'copy'

	input:
	set val(sample), file(pe_bam) from paired_bam
	file frag_file from res_frag_file.collect()

	output:
	set val(sample), file("*.validPairs") into valid_pairs // all valid ligation products
	set val(sample), file("*.DEPairs") into de_pairs // list of dangling ends
	set val(sample), file("*.SCPairs") into sc_pairs // list of self circles
	set val(sample), file("*.REPairs") into re_pairs // list of religation products
	set val(sample), file("*.FiltPairs") into filt_pairs // list of filtered pairs
	set val(sample), file("*RSstat") into all_rsstat // summary stats of number of read pairs in each of above

	script: 
	"""
	/HiC-Pro-devel/scripts/mapped_2hic_fragments.py -f ${frag_file} -r ${pe_bam} --all
	"""
}

/* 
 * Step 4: Build Matrix 
 */ 

// Remove PCR duplicates 

process remove_PCR_duplicates { 
	tag "$sample" 
	publishDir "${params.outdir}/hic_results/data", mode: 'copy' 

	input: 
	set val(sample), file(vpairs) from valid_pairs.groupTuple()

	output: 
	set val(sample), file("*.allValidPairs") into all_valid_pairs 
	file("stats/") into all_mergestat // all filtering statistics 

	script: 
	""" 
	mkdir -p stats/${sample} 
	cat ${vpairs} > ${sample}.allValidPairs 
	echo -n "valid_interaction\t" > stats/${sample}/${sample}_allValidPairs.mergestat 
	cat ${vpairs} | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat 
	echo -n "valid_interaction_rmdup\t" >> stats/${sample}/${sample}_allValidPairs.mergestat 
	cat ${sample}.allValidPairs | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat 

	## Count short range (<20000) vs long range contacts 
	awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}}\$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${sample}.allValidPairs >> stats/${sample}/${sample}_allValidPairs.mergestat
	""" 
} 

// Merge statistic files on filtering and duplicate removal 
process merge_sample { 
	tag "$ext" 
	publishDir "${params.outdir}/hic_results/stats/${sample}", mode: 'copy' 
	
	input: 
	set val(prefix), file(fstat) from all_mapstat.groupTuple().concat(all_pairstat.groupTuple(), all_rsstat.groupTuple())
	
	output: 
	file("mstats/") into all_mstats // summary of all filtering statistics 
	
	script:
	sample = prefix.toString() - ~/(_R1$|_R2$|_val_1$|_val_2$|_1$|_2$)/
	if ( (fstat =~ /.mapstat/) ){ ext = "mmapstat" }
	if ( (fstat =~ /.pairstat/) ){ ext = "mpairstat" }
	if ( (fstat =~ /.RSstat/) ){ ext = "mRSstat" }
	"""
	mkdir -p mstats/${sample}
	/merge_statfiles_new.py -f ${fstat} > mstats/${sample}/${prefix}.${ext}
	""" 
} 

// Split genome into equal sized bins of specified size and generate raw contact matrix 
// Note: the contact matrix is stored in the form: bin_i / bin_j / counts_ij 
// (such that the matrix is made up of three columns) 

process build_contact_matrix {
	tag "$sample - $mres"
	publishDir "${params.outdir}/hic_results/matrix/raw", mode: 'copy'

	input:
	set val(sample), file(vpairs), val(mres) from all_valid_pairs.combine(map_res) // mapping resolution
	file chrsize from chromosome_size.collect() // chromosome size

	output:
	file("*.matrix") into raw_maps // raw contact matrix
	file "*.bed" // BED file with genomic coordinates of bins

	script:
	"""
	/HiC-Pro-devel/scripts/build_matrix --matrix-format upper --binsize ${mres} --chrsizes ${chrsize} --ifile ${vpairs} --oprefix ${sample}_${mres}
	"""
}

/* 
 * Step 5: Normalise Matrix (using ICE) 
 */ 

// Run iterative correction normalisation via ICE to normalise for unknown biases 

/* 
 * Command flags:
 *
 *		--filter_low_counts_perc	percentage of reads to filter out 
 *		--max_iter					maximumum number of matrix-balancing iterations
 *		--eps						relative increment before declaring convergenece 
 *		--remove-all-zeros-loci		remove non-interacting loci 
 *		--output-bias				output the bias vector 
 */ 

process ice_normalise {
	tag "$rmaps" 
	publishDir "${params.outdir}/hic_results/matrix/iced", mode: 'copy' 

	input: 
	file(rmaps) from raw_maps 

	output: 
	file("*iced.matrix") into iced_maps // genome wide contact map after ICE nomalisation 

	script: 
	// extract sample name 
	prefix = rmaps.toString() - ~/(\.matrix)?$/
	""" 
	/HiC-Pro-devel/scripts/src/ice_mod/iced/scripts/ice \
		--results_filename ${prefix}_iced.matrix \
		   --filter_high_counts_perc 0 \
		   --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 --verbose 1 ${rmaps}
	"""
}

/* 
 * Step 6: MultiQC Report 
 */ 

// Run multiqc to visualise mapping and contact statistics 

process multiqc { 
	publishDir "${params.outdir}/MultiQC" 

	input: 
	file ("input_*/*") from all_mstats.concat(all_mergestat).collect() 

	output: 
	file "*multiqc_report.html" into multiqc_report
	file "*_data"

	script: 
	""" 
	multiqc ./ -f -o "${params.outdir}/MultiQC" 
	"""
}
