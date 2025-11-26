rule microbemod_call_methylation:
    """
    Call methylation in nanopore samples.
    """
    input: 
        sorted_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam",
        reference = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'reference_path'], 
    output: 
        scratch_dict['microbemod_call_methylation'] / "{sample}_motifs.tsv", 
    conda:  
        config["MicrobeMod"]["conda env name"]  
    params: 
        microbemod_path = config["MicrobeMod"]["installation path"], 
        methylation_confidence = config["MicrobeMod"]["methylation_confidence_threshold"], 
        percent_methylation_cutoff = config["MicrobeMod"]["percent_methylation_cutoff"], 
        percent_cutoff_streme = config["MicrobeMod"]["percent_cutoff_streme"], 
    shell:
        """
        {params.microbemod_path}/bin/MicrobeMod call_methylation \
            --methylation_confidence_threshold {params.methylation_confidence} \
            --percent_methylation_cutoff {params.percent_methylation_cutoff} \
            --percent_cutoff_streme {params.percent_cutoff_streme} \
            --bam_file {input.sorted_bam} \
            --reference_fasta {input.reference} \
            --threads {resources.cpus_per_task} \
            --output_prefix {wildcards.sample} \
            --output_directory $(dirname {output})

        # call_methylation doesn't produce output if no methylated sites identified. 
        touch {output}
        """

rule microbemod_annotate_rm:
    """
    Annotate reference genome restriction modifications. 
    """
    input: lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'reference_path'], 
    output: scratch_dict['microbemod_annotate_rm'] / "{sample}.rm.genes.tsv", 
    conda: config["MicrobeMod"]["conda env name"]
    params: microbemod_path = config["MicrobeMod"]["installation path"], 
    shell:
        """
        {params.microbemod_path}/bin/MicrobeMod annotate_rm \
            --fasta {input} \
            --output_prefix {wildcards.sample} \
            --output_directory $(dirname {output}) \
            --threads {resources.cpus_per_task}

        touch {output}
        """

rule parse_microbemod:
    """
    Parse and combine outputs from MicrobeMod annotate_rm & call_methylation
    """
    input: 
        call_methylation_output = expand(scratch_dict['microbemod_call_methylation'] / "{sample}_motifs.tsv", sample=SAMPLES), 
        annotate_rm_output = expand(scratch_dict['microbemod_annotate_rm'] / "{sample}.rm.genes.tsv", sample=SAMPLES), 
    output: 
        parsed_call_methylation = results_dict['parsed_call_methylation'],
        parsed_annotate_rm =  results_dict['parsed_annotate_rm'],
    conda: 
        "../envs/data_parse.yaml"
    script: 
        "../scripts/parse_microbemod.py"