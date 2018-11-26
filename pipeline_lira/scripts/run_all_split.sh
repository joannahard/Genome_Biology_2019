export PATH=/media/box2/Experiments/Joanna/bin/miniconda3/envs/LiRA/bin:/media/box2/Experiments/Joanna/bin/miniconda3/envs/snake2/bin:$PATH
export PATH=$PATH:/media/box2/Experiments/Joanna/LiRA/LiRA
export LIRA_DIR=/media/box2/Experiments/Joanna/LiRA/LiRA



#nohup snakemake -p -j 22 lira_output/110/split_all_chrom.chk > nohup_out/split_110.txt &
nohup snakemake -p -j 4 lira_output/111/split_all_chrom.chk > nohup_out/split_111.txt &
nohup snakemake -p -j 4 lira_output/bulk/split_all_chrom.chk > nohup_out/split_bulk.txt &
nohup snakemake -p -j 4 lira_output/112/split_all_chrom.chk > nohup_out/split_112.txt &
nohup snakemake -p -j 4 lira_output/113/split_all_chrom.chk > nohup_out/split_113.txt &
nohup snakemake -p -j 4 lira_output/114/split_all_chrom.chk > nohup_out/split_114.txt &
nohup snakemake -p -j 4 lira_output/118/split_all_chrom.chk > nohup_out/split_118.txt &
nohup snakemake -p -j 4 lira_output/119/split_all_chrom.chk > nohup_out/split_119.txt &
nohup snakemake -p -j 4 lira_output/120/split_all_chrom.chk > nohup_out/split_120.txt &
nohup snakemake -p -j 4 lira_output/121/split_all_chrom.chk > nohup_out/split_121.txt &
nohup snakemake -p -j 4 lira_output/134/split_all_chrom.chk > nohup_out/split_134.txt &

nohup snakemake -p -j 4 lira_output/135/split_all_chrom.chk > nohup_out/split_135.txt &
nohup snakemake -p -j 4 lira_output/136/split_all_chrom.chk > nohup_out/split_136.txt &
nohup snakemake -p -j 4 lira_output/137/split_all_chrom.chk > nohup_out/split_137.txt &
nohup snakemake -p -j 4 lira_output/138/split_all_chrom.chk > nohup_out/split_138.txt &
nohup snakemake -p -j 4 lira_output/142/split_all_chrom.chk > nohup_out/split_142.txt &
nohup snakemake -p -j 4 lira_output/143/split_all_chrom.chk > nohup_out/split_143.txt &
nohup snakemake -p -j 4 lira_output/144/split_all_chrom.chk > nohup_out/split_144.txt &
nohup snakemake -p -j 4 lira_output/145/split_all_chrom.chk > nohup_out/split_145.txt &

