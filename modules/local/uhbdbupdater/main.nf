process UHBDBUPDATER {
    tag "${meta.genus}"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/781ddeab29aecce0f483b4e655218405f1335194a72a0ac638253b9b178bb1b6/data"
    // Singularity: https://wave.seqera.io/view/builds/bd-b9a71dceac28ddf5_1?_gl=1*1muplwv*_gcl_au*MTI1MzgxOTA5MC4xNzY4MjM1MzM1
    storeDir "${params.output_dir}/${meta.genus}"

    input:
    tuple val(meta), path(local_txt), path(meta_tsv), path(objects_txt), path(urls_txt)
    path(uhbdb_dir)

    output:
    tuple val(meta), path("${meta.genus}_k31_s1000.skm")    , emit: skm
    tuple val(meta), path("${meta.genus}_k31_s1000.skd")    , emit: skd
    tuple val(meta), path("${meta.genus}.agc")              , emit: agc
    tuple val(meta), path("${meta.genus}.metadata.tsv.gz")  , emit: tsv
    tuple val(meta), path(".command.log")                   , emit: log
    tuple val(meta), path(".command.sh")                    , emit: script

    script:
    """
    ### Download fastas (aria2c)
    # remove quotes from urls file if present
    sed -i 's/"//g' ${urls_txt}

    # download with aria2c
    aria2c \\
        --input=${urls_txt} \\
        --dir=url_fastas \\
        --max-concurrent-downloads=8 \\
        --max-tries=5 \\
        --retry-wait=60 \\
        --continue


    ### Rename local fastas to match IDs in TSV
    # create directory for renamed local fastas
    mkdir -p local_fastas

    # create symbolic links for local fastas with new names based on TSV
    while IFS=' ' read -r target link; do
    if [[ -n "\$target" && -n "\$link" ]]; then
        # Create the symbolic link
        ln -s "\$target" "\$link"
    fi
    done < "${local_txt}"


    ### Create sketchlib sketch (sketchlib)
    # create input file of fastas
    find ./ -path '*_fastas/*' > fasta_files.txt

    # create tsv with file base name for 1st column then path for second column
    awk -F/ '{print \$NF"\\t"\$0}' fasta_files.txt > fasta_files.tsv
    sed -i 's/\\.f.*\\t/\\t/g' fasta_files.tsv

    # create sketchlib sketch
    sketchlib sketch \\
        -f fasta_files.tsv \\
        -o ${meta.genus}_k31_s1000 \\
        -k 31 \\
        --sketch-size 1000 \\
        --threads ${task.cpus} \\
        --verbose

    ### IF UHBDB FILES DO NOT EXIST, CREATE NEW ARCHIVES ###
    if [[ ! -f ${uhbdb_dir}/${meta.genus}/${meta.genus}_k31_s1000.skm ]]; then
        ### SELF: all-v-all dist (sketchlib)
        # run sketchlib dist
        sketchlib dist \\
            ${meta.genus}_k31_s1000.skm \\
            -o ${meta.genus}_self_dist.tsv \\
            -k 31 \\
            --ani \\
            --threads ${task.cpus} \\
            --verbose \\
            --knn 100


        ### SELF: Extract unique sequences (clusty)
        # add header for clusty
        sed -i '1i query\\tref\\tani' ${meta.genus}_self_dist.tsv

        # run clusty greedy clustering and retain highest n50
        clusty \\
            ${meta.genus}_self_dist.tsv \\
            ${meta.genus}_self_unique_cdhit.tsv \\
            --objects-file ${objects_txt} \\
            --similarity \\
            --min ani ${params.cluster_ani} \\
            --out-representatives


        # Create output directory for final outputs
        mkdir -p final_files

        ### FINAL: Move sketchlib sketch to final outputs directory (sketchlib)
        mv ${meta.genus}_k31_s1000.skm final_files/${meta.genus}_k31_s1000.skm
        mv ${meta.genus}_k31_s1000.skd final_files/${meta.genus}_k31_s1000.skd
        
        ### FINAL: Create file of unique fasta files (polars)
        unique_fastas.py \\
            --input fasta_files.tsv \\
            --clusters ${meta.genus}_self_unique_cdhit.tsv \\
            --output ${meta.genus}_self_unique_fastas


        ### FINAL: Create AGC archive of unique fastas (agc)
        ref_fasta=\$(head -n 1 ${meta.genus}_self_unique_fastas.ref.txt)
        
        # create AGC archive with unique fastas
        agc create \\
            \$ref_fasta \\
            -i ${meta.genus}_self_unique_fastas.unique.txt \\
            -a true \\
            -b 500 \\
            -s 1500 \\
            -o final_files/${meta.genus}.agc \\
            -t ${task.cpus} \\
            -f 0.01 \\
            -v 1


        ### FINAL: Create metadata file for unique fastas (polars)
        metadata.py \\
            --input ${meta_tsv} \\
            --clusters ${meta.genus}_self_unique_cdhit.tsv \\
            --output final_files/${meta.genus}.metadata.tsv

        # Compress
        gzip final_files/${meta.genus}.metadata.tsv

        ### Cleanup
        rm -rf url_fastas/ local_fastas/ fasta_files.txt fasta_files.tsv \\
            ${meta.genus}_k31_s1000* ${meta.genus}_self_dist.tsv ${meta.genus}_self_unique_cdhit.tsv \\
            ${meta.genus}_rm_list.txt ${meta.genus}_self_unique_fastas*

        mv final_files/* ./

    else
        # Create output directory for final outputs
        mkdir -p final_files

        ### NEW2OLD: Compare new unique to existing unique (sketchlib)
        sketchlib merge \\
            ${meta.genus}_k31_s1000.skd \\
            ${uhbdb_dir}/${meta.genus}/${meta.genus}_k31_s1000.skd \\
            -o final_files/${meta.genus}_k31_s1000 \\
            --verbose

        ### NEW2OLD: Create objects file with new + old unique sequences (polars)
        combined_objects.py \\
            --new_metadata ${meta_tsv} \\
            --old_metadata ${uhbdb_dir}/${meta.genus}/${meta.genus}.metadata.tsv.gz \\
            --output ${meta.genus}_combined_objects.txt

        tail -n +2 ${meta.genus}_combined_objects.txt > ${meta.genus}_combined_objects_ids.txt

        sketchlib dist \\
            final_files/${meta.genus}_k31_s1000 \\
            -o ${meta.genus}_combined_dist.tsv \\
            -k 31 \\
            --ani \\
            --threads ${task.cpus} \\
            --knn 100 \\
            --subset ${meta.genus}_combined_objects_ids.txt \\
            --verbose

        ### NEW2OLD: Extract unique sequences (clusty)
        # add header for clusty
        sed -i '1i query\\tref\\tani' ${meta.genus}_combined_dist.tsv

        # run clusty greedy clustering and retain highest n50
        clusty \\
            ${meta.genus}_combined_dist.tsv \\
            ${meta.genus}_combined_cdhit.tsv \\
            --objects-file ${meta.genus}_combined_objects.txt \\
            --similarity \\
            --min ani ${params.cluster_ani} \\
            --out-representatives

        
        ### NEW2OLD: Remove redundant sketches from new + old sketches (sketchlib)
        # identify duplicate sketches (polars)
        update.py \\
            --old_metadata ${uhbdb_dir}/${meta.genus}/${meta.genus}.metadata.tsv.gz \\
            --new_metadata ${meta_tsv} \\
            --new_fasta_tsv fasta_files.tsv \\
            --combined_clusters ${meta.genus}_combined_cdhit.tsv \\
            --output_metadata ${meta.genus}.metadata.tsv \\
            --output_new_unique_fastas ${meta.genus}_new_unique_fastas.txt

        ### Compress
        gzip ${meta.genus}.metadata.tsv

        ### FINAL: Append new unique sequences to existing AGC archive (agc)
        agc append \\
            ${uhbdb_dir}/${meta.genus}/${meta.genus}.agc \\
            -i ${meta.genus}_new_unique_fastas.txt \\
            -a true \\
            -o final_files/${meta.genus}.agc \\
            -t ${task.cpus} \\
            -f 0.01 \\
            -v 1

        ### Cleanup
        rm -rf url_fastas/ local_fastas/ fasta_files.txt fasta_files.tsv \\
            ${meta.genus}_k31_s1000* ${meta.genus}_combined* ${meta.genus}_new2old_dist.tsv \\
            ${meta.genus}_combined_cdhit.tsv ${meta.genus}_combined_objects.txt \\
            ${meta.genus}_new_unique_fastas.txt

        mv final_files/* ./
    fi
    """
}