# Paired end
aws batch submit-job \
    --job-name nf-midas-0223-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-midas,\
"--db_midas","s3://dev-scratch/ReferenceDBs/Midas/v1.2/midas_db_v1.2.tar.gz",\
"--db_knead","s3://dev-scratch/ReferenceDBs/Biobakery/kneaddata/human_genome/bowtie2",\
"--output_folder","s3://genomics-workflow-core/Results/midas/hCom-B6Gen",\
"--manifest","s3://genomics-workflow-core/Results/midas/hCom-B6Gen/00_seedfile/hCom-B6Gen.seedfile.csv"
