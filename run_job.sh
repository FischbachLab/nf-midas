# TEST
aws batch submit-job \
    --job-name nf-midas-0225-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-midas,\
"--manifest","s3://genomics-workflow-core/Results/midas/00_TEST/00_seedfile/test.seedfile.csv"

# Paired end
aws batch submit-job \
    --job-name nf-midas-0224-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-midas,\
"--project","hCom-B6Gen",\
"--manifest","s3://genomics-workflow-core/Results/midas/hCom-B6Gen/00_seedfile/hCom-B6Gen.seedfile.csv"