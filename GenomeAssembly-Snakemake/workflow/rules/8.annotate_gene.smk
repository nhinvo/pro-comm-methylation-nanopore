rule annotate_gene:
    input: 
        scratch_dict["metabat_binning"] / "bins" / "{sample}",
    output:
        scratch_dict['gene_annotation'] / "{sample}" / "{sample}_done.txt",
    conda: 
        "../envs/prokka.yaml"
    shell:
        """
        # run prokka on each bin in the directory 
        for file in {input}/*; do 
            echo Running prokka on: $file

            filename=$(basename "$file" .fa)

            prokka \
                --cpus {resources.cpus_per_task} \
                --force --prefix $filename \
                --outdir $(dirname {output}) \
                $file

            echo Completed prokka on $filename. 
        done

        touch {output}
        """

rule annotate_gene_eggnog:
    input:
        scratch_dict['gene_annotation'] / "{sample}", 
    output:
        scratch_dict['annotate_gene_eggnog'] / "{sample}" / "{sample}_done.txt",
    conda:
        "../envs/eggnog.yaml"
    params:
        eggnogg_db = config['database']['eggnogg database'], 
        temp_dir = scratch_dict['annotate_gene_eggnog'], 
    shell:
        """
        # run eggnog on each bin in the sample directory 
        for file in {input}/*.faa; do 
            echo Running eggnog on .faa: $file

            filename=$(basename "$file" .faa)

            emapper.py \
                --cpu {resources.cpus_per_task} --override \
                -i $file --itype proteins \
                --decorate_gff yes --report_no_hits \
                --temp_dir {params.temp_dir} \
                --data_dir {params.eggnogg_db} \
                --output $filename \
                --output_dir $(dirname {output})

            echo Completed eggnog on $filename. 
        done

        touch {output}
        """

rule finalize_annot:
    input: expand(scratch_dict['gene_annotation'] / "{sample}" / "{sample}_done.txt", sample=SAMPLES),
    output: results_dict['gff_path_table'],
    conda: "../envs/data.yaml"
    script: "../scripts/finalize_annot.py"