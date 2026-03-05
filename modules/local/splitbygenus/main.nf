process SPLITBYGENUS {
    label 'process_high'
    tag "split_by_genus"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a6b4319d6bba30fc01cff76aa7f26a0f43592697b8819bc7dcd0b50e61ac184e/data"
    storeDir "${params.tmp_dir}/splitbygenus"

    input:
    path(tsv)
    val(large_threshold)

    output:
    path("large_genera/*")  , emit: large   , optional: true
    path("small_genera/*")  , emit: small   , optional: true
    path(".command.log")    , emit: log
    path(".command.sh")     , emit: script

    script:
    """
    # Make output directories
    mkdir -p large_genera
    mkdir -p small_genera

    ### Split by input genus
    split_by_genus.py \\
        --input ${tsv} \\
        --large_threshold ${large_threshold}
    """
}
