rule read_mapping: 
    """
    Read mapping using dorado 
    """
    input: 
        read_path = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'read_path'],
        reference = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'reference_path'], 
    output: 
        temp(scratch_dict["read_mapping"] / "{sample}.bam"), 
    params: 
        dorado_path = config["dorado"]["installation path"], 
    shell:
        """
        {params.dorado_path}/bin/dorado aligner \
            {input.reference} {input.read_path} \
            --threads {resources.cpus_per_task} \
            --verbose \
            > {output}
        """

rule samtools_sort_index:
    """
    Description: sort index bam. 
    """
    input: 
        scratch_dict["read_mapping"] / "{sample}.bam", 
    output: 
        sorted_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
        indexed_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam.bai", 
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        # sort bam 
        samtools sort \
            --threads {resources.cpus_per_task} \
            --output-fmt BAM \
            -o {output.sorted_bam} \
            {input}

        # index the sorted_bam
        samtools index \
            --threads {resources.cpus_per_task} \
            --bai \
            --output {output.indexed_bam} \
            {output.sorted_bam}
        """