#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Usage
def help_msg() {
	help = """[nf-mag-depths]: Calculate depth of coverage from a metagenomic set of bins.
			|
			|Usage:
			|  --input  comma-separated csv file with header: {sample,read1,read2}
			|           [default: none]
			|
			|Optional arguments:
			|  --outdir results directory
			|           [default: ${params.outdir}]
			|
			|Print help:
			|  nextflow run nf-mag-depths/main.nf -c nf-mag-depths/mag_depths.config --help
			|
			""".stripMargin()
	println help
	exit 0
}

process BWAMEM2 {
	tag "${sample}"
	publishDir "${params.outdir}/bwa-mem2", mode: "copy"
	
	input:
	tuple val(sample), file(fasta), file(read1), file(read2)
	
	output:
	tuple val(sample), path("${sample}"), emit: bam
	
	script:
	def prefix = task.ext.prefix ?: "${sample}"
	//prefix = prefix.replaceAll(/.fna/, "")
	
	if (params.debug == true) {
		//test mode
		"""
		mkdir ${prefix}
		touch ${prefix}.${params.ext}.0123
		touch ${prefix}.${params.ext}.amb
		touch ${prefix}.${params.ext}.ann
		touch ${prefix}.${params.ext}.bwt.2bit.64
		touch ${prefix}.${params.ext}.pac
		echo bwa-mem2 index ${fasta}
		
		touch ${prefix}.log
		touch ${prefix}.bam
		touch ${prefix}_sorted.bam
		#echo bwa-mem2 mem -t $task.cpus ${fasta} ${read1} ${read2} 2> ${prefix}.log | samtools view --threads $task.cpus -bS -o ${prefix}.bam
		#echo samtools sort --threads $task.cpus ${prefix}.bam -o ${prefix}_sorted.bam
		
		mv ${prefix}_sorted.bam ${prefix}.log ${prefix}
		"""
	}
	else {
		//run
		"""
		mkdir ${prefix}
		
		#index
		bwa-mem2 index ${fasta}
		#mapping
		#INDEX=`find -L ${prefix} -name "*.amb" | sed 's/\\.amb\$//'`
		bwa-mem2 mem -t $task.cpus ${fasta} ${read1} ${read2} 2> ${prefix}.log | samtools view --threads $task.cpus -bS -o ${prefix}.bam
		samtools sort --threads $task.cpus ${prefix}.bam -o ${prefix}_sorted.bam
		
		mv ${prefix}_sorted.bam ${prefix}.log ${prefix}
		"""
	}
}

process CALC_DEPTHS {
	tag "${params.fasta}"
	publishDir "${params.outdir}/depths", mode: "copy"
	
	container "quay.io/biocontainers/metabat2:2.15--h986a166_1"
	if (workflow.containerEngine == "singularity") {
		container "https://depot.galaxyproject.org/singularity/metabat2:2.15--h986a166_1"
	}
	
	input:
	file(fasta)
	file(bams)
	
	output:
	path("samples_depth.txt"), emit: depth
	
	script:
	if (params.debug == true) {
		//test mode
		"""
		touch samples_depth.txt
		echo jgi_summarize_bam_contig_depths --outputDepth samples_depth.txt --referenceFasta ${fasta} ${bams}
		"""
	}
	else {
		//run
		"""
		export OMP_NUM_THREADS=$task.cpus
		
		jgi_summarize_bam_contig_depths --outputDepth samples_depth.txt --referenceFasta ${fasta} ${bams}
		"""
	}
}

// Main
workflow {
	// Check arguments
	if (params.help || params.input == null) {
		help_msg()
	}
	
	// Channel for -merged- fasta
	input_ch = Channel.fromPath(params.input)
		.splitCsv(header: true)
		.map { row -> tuple(row.sample, file(params.fasta), file(row.read1), file(row.read2)) }
		//.view()
	
	// Index & mapping fasta
	BWAMEM2( input_ch )
	
	// Collect bam files of all samples
	bams = BWAMEM2.out.bam
		.map { sample, dir -> file(dir.toString() + '/' + sample + '_sorted.bam') } //.map { sample, dir -> dir }
		.collect()
		//.view()
	
	// Generate coverage depth across all samples
	CALC_DEPTHS( file(params.fasta), bams )
	//output a txt: parser required?
}
//nextflow -bg run nf-mag-depths/main.nf -c nf-mag-depths/mag_depths.config --input input.csv -profile singularity {-resume}
