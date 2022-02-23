# Paired end
NXF_VER=20.01.0 \
    nextflow \
    run \
    -c nextflow.config \
    -profile testing \
    midas_workflow.nf \
    --manifest test_data/manifest.paired.csv \
    --db db/ \
    --output_folder test_output/paired/ \
    --species_cov 0.01 \
    -with-docker ubuntu:18.04 \
    -w work/ \
    -process.executor local \
    -with-report \
    -with-trace \
    -resume

aws batch submit-job \
    --job-name nf-midas-0223-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/midas_nextflow,\
"--db_midas","s3://dev-scratch/ReferenceDBs/Midas/v1.2/midas_db_v1.2.tar.gz",\
"--db_knead","/mnt/efs/databases/Biobakery/kneaddata"
"--output_folder","s3://genomics-workflow-core/Results/midas/hCom-B6Gen",\
"--manifest","s3://genomics-workflow-core/Results/midas/hCom-B6Gen/00_seedfile/hCom-B6Gen.seedfile.csv"

