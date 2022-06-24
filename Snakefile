
rule all:
    input:
        "plots/pbmc/bulk_scatter.png",
        "plots/pbmc/fragments.png",
        "plots/bmmc/fragments_bmmc.png",
        "plots/bmmc/hmap_ac.png",
        "plots/bmmc/trajectory.png",
        "plots/bmmc/rna_repressed.pdf",
        "plots/hek_k562/scatterplots_bulk.png"
        "plots/hek_k562/correlation.png",
        "plots/hek_k562/multi_cuttag_frip.png",
        "plots/hek_k562/multi_cuttag_scatter.png",
        "plots/pbmc/encode_cor_all.png",
        "plots/pbmc/encode_cor_spearman.png",
        "plots/pbmc/replicate_correlation.png"

### Downloads ###
rule get_chr_size:
    output: "data/hg38.chrom.sizes"
    shell:
        """
        cd data
        wget https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
        """

rule get_hg38:
    output:
        "genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.0123",
        "genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    threads: 1
    message: "Download hg38 genome"
    shell:
        """
        cd genome/hg38
        aws s3 sync s3://stuart-genomes/hg38_analysis/bwa-mem2/ .
        aws s3 cp s3://stuart-genomes/hg38_analysis/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz .
        """

rule get_ct_pro:
    input: "datasets/ct_pro.txt"
    output:
        "data/ct_pro/H3K27ac_fragments.tsv.gz",
        "data/ct_pro/H3K27me3_fragments.tsv.gz"
    message: "Downloading CUT&Tag-pro data"
    threads: 1
    shell:
        """
        wget -i {input} -P data/ct_pro
        """

rule get_mesc_peaks:
    input: "datasets/mesc.txt"
    output: "data/mesc/ENCFF008XKX.bed.gz", "data/mesc/ENCFF360VIS.bed.gz"
    threads: 1
    message: "Downloading mESC peaks"
    shell:
        """
        wget -i {input} -P data/mesc
        """
        
rule get_henikoff:
    input: "datasets/henikoff.txt"
    output:
        "data/henikoff/GSM5034342_K27me3_R1_PBMC.fragments.HG38.tsv.gz",
        "data/henikoff/GSM5034343_K27me3_R2_PBMC.fragments.HG38.tsv.gz",
        "data/henikoff/GSM5034344_K27ac_PBMC.fragments.HG38.tsv.gz"
    message: "Downloading PBMC CUT&Tag data"
    threads: 1
    shell:
        """
        wget -i {input} -P data/henikoff
        wget https://github.com/Henikoff/scCUT-Tag/files/8647406/meta.csv -P data/henikoff
        wget https://github.com/Henikoff/scCUT-Tag/files/8653804/K27Ac_pbmc_meta.csv -P data/henikoff

        cd data/henikoff
        tabix -p bed GSM5034342_K27me3_R1_PBMC.fragments.HG38.tsv.gz
        tabix -p bed GSM5034343_K27me3_R2_PBMC.fragments.HG38.tsv.gz
        tabix -p bed GSM5034344_K27ac_PBMC.fragments.HG38.tsv.gz
        """
        
rule update_ct_pro:
    input:
        "data/ct_pro/H3K27ac_fragments.tsv.gz",
        "data/ct_pro/H3K27me3_fragments.tsv.gz"
    output:
        "data/ct_pro/H3K27ac_updated.rds",
        "data/ct_pro/H3K27me3_updated.rds"
    message: "Fix CUT&Tag-pro paths"
    threads: 1
    shell:
        """
        Rscript code/CT_pro/update_objects.R
        """

rule get_multi_ct:
    output: "data/multict/fragments/H3K27ac-H3K27ac.tsv.gz"
    message: "Downloading multiCUT&Tag data"
    threads: 1
    shell:
        """
        cd data/multict/fragments
        aws s3 sync s3://multi-cuttag/fragments/ . --request-payer
        """

rule get_k562_atac:
    output: "data/K562_ATAC/ENCFF006OFA.bigBed"
    message: "Downloading K562 ENCODE data"
    shell:
        """
        cd data/K562_ATAC
        wget https://www.encodeproject.org/files/ENCFF006OFA/@@download/ENCFF006OFA.bigBed
        wget https://www.encodeproject.org/files/ENCFF600FDO/@@download/ENCFF600FDO.bigWig
        wget https://www.encodeproject.org/files/ENCFF558BLC/@@download/ENCFF558BLC.bed.gz
        """

rule download_bmmc_atac:
    input: "datasets/bmmc_atac.txt"
    output: touch("data/bmmc_atac/download.done")
    threads: 1
    message: "Download BMMC scATAC-seq hg38 fragment files"
    shell:
        """
        wget -i {input} -P data/bmmc_atac
        aws s3 sync s3://mpal-hg38/public/ ./data/bmmc_atac/ --request-payer
        """

rule download_bmmc_rna:
    input: "datasets/bmmc_rna.txt"
    output: "objects/fullref.Rds"
    shell:
        """
        aws s3 sync s3://bmmc-reference/public objects/ --request-payer
        """

rule download_encode:
    input: "datasets/encode.txt"
    output:
        "data/encode/ENCFF842JLZ.bigWig",
        "data/encode/ENCFF138FVQ.bed.gz",
        "data/encode/ENCFF291LVP.bed.gz",
        "data/encode/ENCFF832RWT.bed.gz",
    message: "Downloading ENCODE PBMC data"
    shell:
        """
        wget -i {input} -P data/encode
        """

rule collect_encode_peaks:
    input: "data/encode/ENCFF138FVQ.bed.gz"
    output:
        "data/encode/h3k27me3.bed",
        "data/encode/h3k27ac.bed",
        "data/encode/all.bed",
        "data/encode/all_bulk.bed"
    threads: 1
    shell:
        """
        Rscript code/encode/combine_peaks.R
        """

rule get_scatac_pbmc:
    output: "data/pbmc_atac/10k_pbmc_ATACv2_nextgem_Chromium_X_fragments.tsv.gz"
    message: "Downloading PBMC scATAC-seq"
    threads: 1
    shell:
        """
        cd data/pbmc_atac
        wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_X/10k_pbmc_ATACv2_nextgem_Chromium_X_fragments.tsv.gz
        wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_X/10k_pbmc_ATACv2_nextgem_Chromium_X_fragments.tsv.gz.tbi
        wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_X/10k_pbmc_ATACv2_nextgem_Chromium_X_singlecell.csv
        wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_X/10k_pbmc_ATACv2_nextgem_Chromium_X_filtered_peak_bc_matrix.h5
        """

### Map ###

rule map_cellculture_sc:
    input:
        r1="data/HEK_K562_sc/dmx/scCmix_{mark}_R1.barcoded.fastq.gz",
        r2="data/HEK_K562_sc/dmx/scCmix_{mark}_R3.barcoded.fastq.gz",
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output:
        frags="data/HEK_K562_sc/sinto/{mark}.bed.gz",
        bam="data/HEK_K562_sc/mapped/{mark}.bam"
    threads: 24
    shell:
        """
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            {input.r1} \
            {input.r2} \
            | samtools sort -@ {threads} -O bam - \
            > data/HEK_K562_sc/mapped/{wildcards.mark}.bam
        
        samtools index data/HEK_K562_sc/mapped/{wildcards.mark}.bam
        
        sinto fragments \
          -b data/HEK_K562_sc/mapped/{wildcards.mark}.bam \
          -p {threads} \
          -f data/HEK_K562_sc/sinto/{wildcards.mark}.tmp \
          --barcode_regex "[^:]*"
        
        sort -k1,1 -k2,2n data/HEK_K562_sc/sinto/{wildcards.mark}.tmp > data/HEK_K562_sc/sinto/{wildcards.mark}.bed
        bgzip -@ {threads} data/HEK_K562_sc/sinto/{wildcards.mark}.bed
        tabix -p bed {output.frags}
        rm data/HEK_K562_sc/sinto/{wildcards.mark}.tmp
        """

rule map_cellculture_bulk:
    input:
        r1="data/HEK_K562_bulk/{cell}-{plex}-{mark}_1.fastq.gz",
        r2="data/HEK_K562_bulk/{cell}-{plex}-{mark}_2.fastq.gz",
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output: "data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bam"
    threads: 24
    message: "Mapping {input.r1}"
    shell:
        """
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            {input.r1} \
            {input.r2} \
            | samtools sort -@ {threads} -O bam - \
            > data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.bam
        samtools index data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.bam
        """
        
rule create_bulk_cellculture_bigwig:
    input: "data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bam"
    output: "data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bw"
    threads: 6
    shell:
        """
        bamCoverage -b {input} -o {output} -p {threads} --normalizeUsing BPM
        """

rule create_bulk_cellculture_fragments:
    input:"data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bam"
    output: "data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bed.gz"
    threads: 6
    shell:
        """
        sinto fragments -p {threads} \
          -b {input} \
          --barcode_regex "[^:]*" \
          -f data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.frag
        
        sort -k1,1 -k2,2n data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.frag \
          > data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.bed
        bgzip data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.bed
        tabix -p bed {output}
        rm data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.frag
        """

rule process_cellculture_sc:
    input:
        ac="data/HEK_K562_sc/sinto/K27ac.bed.gz",
        me="data/HEK_K562_sc/sinto/K27me.bed.gz",
        pol2="data/HEK_K562_sc/sinto/Pol2.bed.gz"
    output: 
        obj="objects/hek_k562.rds",
        frags="data/HEK_K562_sc/split/K27ac_hek.bed.gz"
    message: "Processing HEK/K562 scNTT-seq data"
    threads: 6
    shell:
        """
        Rscript code/HEK_K562_sc/process.R {threads} {output.obj} data/HEK_K562_sc/sinto
        """

rule map_pbmc_bulk:
    input:
        me3_r1_mono="data/pbmc_bulk/PBMCs-MmK27me3_S1_ME_L001_R1_001.fastq.gz",
        me3_r2_mono="data/pbmc_bulk/PBMCs-MmK27me3_S1_ME_L001_R2_001.fastq.gz",
        ac_r1_mono="data/pbmc_bulk/PBMCs-OcK27Ac_S2_ME_L001_R1_001.fastq.gz",
        ac_r2_mono="data/pbmc_bulk/PBMCs-OcK27Ac_S2_ME_L001_R2_001.fastq.gz",
        me3_r1_plex="data/pbmc_bulk/PBMCs-plex-MmK27me3_S3_ME_L001_R1_001.fastq.gz",
        me3_r2_plex="data/pbmc_bulk/PBMCs-plex-MmK27me3_S3_ME_L001_R2_001.fastq.gz",
        ac_r1_plex="data/pbmc_bulk/PBMCs-plex-OcK27Ac_S4_ME_L001_R1_001.fastq.gz",
        ac_r2_plex="data/pbmc_bulk/PBMCs-plex-OcK27Ac_S4_ME_L001_R2_001.fastq.gz",
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output:
        me3_mono="data/pbmc_bulk/mapped/mono/me3.bam",
        me3_plex="data/pbmc_bulk/mapped/plex/me3.bam",
        ac_mono="data/pbmc_bulk/mapped/mono/ac.bam",
        ac_plex="data/pbmc_bulk/mapped/plex/ac.bam"
    threads: 24
    message: "Mapping PBMC bulk NTT-seq"
    shell:
        """
        # me3
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            {input.me3_r1_mono} \
            {input.me3_r2_mono} \
            | samtools sort -@ {threads} -O bam - \
            > data/pbmc_bulk/mapped/mono/me3.bam
        samtools index data/pbmc_bulk/mapped/mono/me3.bam

        bwa-mem2 mem {input.genome} \
            -t {threads} \
            {input.me3_r1_plex} \
            {input.me3_r2_plex} \
            | samtools sort -@ {threads} -O bam - \
            > data/pbmc_bulk/mapped/plex/me3.bam
        samtools index data/pbmc_bulk/mapped/plex/me3.bam
        
        # ac
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            {input.ac_r1_mono} \
            {input.ac_r2_mono} \
            | samtools sort -@ {threads} -O bam - \
            > data/pbmc_bulk/mapped/mono/ac.bam
        samtools index data/pbmc_bulk/mapped/mono/ac.bam
        
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            {input.ac_r1_plex} \
            {input.ac_r2_plex} \
            | samtools sort -@ {threads} -O bam - \
            > data/pbmc_bulk/mapped/plex/ac.bam
        samtools index data/pbmc_bulk/mapped/plex/ac.bam
        """

rule create_bulk_pbmc_bigwig:
    input: "data/pbmc_bulk/mapped/{plex}/{mark}.bam"
    output: "data/pbmc_bulk/mapped/{plex}/{mark}.bw"
    threads: 6
    shell:
        """
        bamCoverage -b {input} -o {output} -p {threads} --normalizeUsing BPM
        """

rule create_bulk_pbmc_fragments:
    input:"data/pbmc_bulk/mapped/{plex}/{mark}.bam"
    output: "data/pbmc_bulk/mapped/{plex}/{mark}.bed.gz"
    threads: 6
    shell:
        """
        sinto fragments -p {threads} \
          -b {input} \
          --barcode_regex "[^:]*" \
          -f data/pbmc_bulk/mapped/{wildcards.plex}/{wildcards.mark}.frag
        
        sort -k1,1 -k2,2n data/pbmc_bulk/mapped/{wildcards.plex}/{wildcards.mark}.frag \
          > data/pbmc_bulk/mapped/{wildcards.plex}/{wildcards.mark}.bed
        bgzip data/pbmc_bulk/mapped/{wildcards.plex}/{wildcards.mark}.bed
        tabix -p bed {output}
        rm data/pbmc_bulk/mapped/{wildcards.plex}/{wildcards.mark}.frag
        """

rule ct_pro_bigwig:
    input:
        me3="data/ct_pro/H3K27me3_fragments.tsv.gz",
        ac="data/ct_pro/H3K27ac_fragments.tsv.gz",
        chrom="data/hg38.chrom.sizes"
    output:
        me3_bg="data/ct_pro/H3K27me3.bg",
        me3_bw="data/ct_pro/H3K27me3.bw",
        ac_bg="data/ct_pro/H3K27ac.bg",
        ac_bw="data/ct_pro/H3K27ac.bw"
    threads: 6
    shell:
        """
        bedtools genomecov -i {input.me3} -g {input.chrom} -bg > data/ct_pro/H3K27me3.tmp
        sort -k1,1 -k2,2n data/ct_pro/H3K27me3.tmp > {output.me3_bg}
        bedGraphToBigWig {output.me3_bg} {input.chrom} {output.me3_bw}
        
        bedtools genomecov -i {input.ac} -g {input.chrom} -bg > data/ct_pro/H3K27ac.tmp
        sort -k1,1 -k2,2n data/ct_pro/H3K27ac.tmp > {output.ac_bg}
        bedGraphToBigWig {output.ac_bg} {input.chrom} {output.ac_bw}
        """

rule encode_k562_peaks:
    input: "datasets/k562.txt"
    output:
        "data/k562_peaks/ENCFF031FSF.bed.gz",
        "data/k562_peaks/ENCFF038DDS.bed.gz",
        "data/k562_peaks/ENCFF266OPF.bed.gz",
        "data/k562_peaks/ENCFF053XYZ.bed.gz"
    threads: 1
    shell:
        """
        wget -i {input} -P data/k562_peaks
        """

rule collect_k562_rna_peaks:
    input: 
        "data/k562_peaks/ENCFF266OPF.bed.gz",
        "data/k562_peaks/ENCFF053XYZ.bed.gz"
    output: "data/k562_peaks/rna.bed"
    threads: 1
    shell:
        """
        Rscript code/HEK_K562_bulk/collect_peaks.R
        """

rule demux_pbmc_ntt:
    input:
        read1="data/{exp}/ntt/Undetermined_S0_R1_001.fastq.gz",
        readi5="data/{exp}/ntt/Undetermined_S0_R2_001.fastq.gz",
        read2="data/{exp}/ntt/Undetermined_S0_R3_001.fastq.gz",
        bc="data/{exp}/barcodes.fa"
    output:
        "data/{exp}/ntt/out/H3K27me3.R1.fastq",
        "data/{exp}/ntt/out/H3K27ac.R1.fastq"
    params:
        outdir="data/{exp}/ntt/out"
    message: "Demultiplex NTT reads for {wildcards.exp}"
    threads: 1
    shell:
        """
        python code/demux.py \
           --read1 {input.read1} \
           --read_i5 {input.readi5} \
           --read2 {input.read2} \
           --tn5 {input.bc} \
           --output {params.outdir}
        """

rule map_pbmc_ntt:
    input:
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        reads="data/{exp}/ntt/out/{mk}.R1.fastq"
    output: "data/{exp}/ntt/{mk}.tsv.gz"
    message: "Mapping PBMC scNTT-seq for {wildcards.exp} {wildcards.mk}"
    threads: 24
    shell:
        """
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            data/{wildcards.exp}/ntt/out/{wildcards.mk}.R1.fastq \
            data/{wildcards.exp}/ntt/out/{wildcards.mk}.R2.fastq \
            | samtools sort -@ {threads} -O bam - \
            > data/{wildcards.exp}/ntt/out/{wildcards.mk}.bam
            
        samtools index data/{wildcards.exp}/ntt/out/{wildcards.mk}.bam
        
        sinto fragments -p {threads} \
          -b data/{wildcards.exp}/ntt/out/{wildcards.mk}.bam \
          --barcode_regex "[^:]*" \
          -f data/{wildcards.exp}/ntt/{wildcards.mk}.frag

        sort -k1,1 -k2,2n data/{wildcards.exp}/ntt/{wildcards.mk}.frag > data/{wildcards.exp}/ntt/{wildcards.mk}.tsv

        rm data/{wildcards.exp}/ntt/{wildcards.mk}.frag
        bgzip -@ {threads} data/{wildcards.exp}/ntt/{wildcards.mk}.tsv
        tabix -p bed {output}
        """

rule index_adt:
    input: "data/totalseq_a.tsv"
    output: directory("data/adt_index")
    message: "Create ADT index"
    threads: 1
    shell:
        """
        salmon index -t {input} -i {output} --features -k7
        """

rule join_adt_reads_pbmc:
    input:
        r1="data/pbmc_protein/adt/TAG_S2_R1_001.fastq.gz",
        r3="data/pbmc_protein/adt/TAG_S2_R3_001.fastq.gz"
    output: "data/pbmc_protein/adt/read1.fastq.gz"
    threads: 4
    message: "Combining ADT reads for quantification"
    shell:
        """
        paste <(gunzip -c {input.r1}) <(gunzip -c {input.r3}) \
          | paste - - - - \
          | awk -F'\\t' -v one="data/pbmc_protein/adt/read1.fastq" '{{ OFS="\\n"; print $1,$3$4,$5,$7$8 >> one;}}'
        
        pigz -p {threads} data/pbmc_protein/adt/read1.fastq
        """

rule map_pbmc_adt:
    input:
        index=directory("data/adt_index"),
        reads="data/pbmc_protein/adt/read1.fastq.gz"
    output: "data/pbmc_protein/adt/outs/alevin/quants_mat.gz"
    message: "Quantify PBMC ADTs for donor1"
    threads: 6
    shell:
        """
        salmon alevin -l ISR -i {input.index} \
            -o data/pbmc_protein/adt/outs \
            -p {threads} \
            --naiveEqclass \
            --keepCBFraction 0.8 \
            -1 data/pbmc_protein/adt/TAG_S2_R2_001.fastq.gz \
            -2 {input.reads} \
            --bc-geometry 1[1-16] \
            --umi-geometry 2[1-10] \
            --read-geometry 2[71-85] \
            --tgMap data/adt_map_totalseq_a.tsv
        """

rule barcode_bmmc_ntt:
    input: "data/bmmc_dual/DNA/scBM_k27ac_S2_R2_001.fastq.gz"
    message: "Attaching BMMC cell barcodes"
    output: directory("data/bmmc_dual/DNA/barcoded")
    threads: 10
    shell:
        """
        cd data/bmmc_dual/DNA
        cutadapt -j {threads} -u 14 -o K27ac_R2.fastq scBM_k27ac_S2_R2_001.fastq.gz
        cutadapt -j {threads} -u 14 -o K27me_R2.fastq scBM_k27me_S1_R2_001.fastq.gz
        
        gzip -d scBM_k27ac_S2_R1_001.fastq.gz scBM_k27ac_S2_R3_001.fastq.gz
        gzip -d scBM_k27me_S1_R1_001.fastq.gz scBM_k27me_S1_R3_001.fastq.gz        
        
        sinto barcode \
          --barcode_fastq K27ac_R2.fastq \
          --read1 scBM_k27ac_S2_R1_001.fastq \
          --read2 scBM_k27ac_S2_R3_001.fastq \
          -b 16
          
        sinto barcode \
          --barcode_fastq K27me_R2.fastq \
          --read1 scBM_k27me_S1_R1_001.fastq \
          --read2 scBM_k27me_S1_R3_001.fastq \
          -b 16
          
        pigz -p {threads} *.fastq
        mkdir barcoded
        mv scBM_k27me_S1_R1_001.barcoded.fastq.gz scBM_k27me_S1_R3_001.barcoded.fastq.gz ./barcoded/
        mv scBM_k27ac_S2_R1_001.barcoded.fastq.gz scBM_k27ac_S2_R3_001.barcoded.fastq.gz ./barcoded/
        """

rule map_bmmc_ntt:
    input:
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        reads="data/bmmc_dual/DNA/barcoded"
    output:
        ac="data/bmmc_dual/DNA/H3K27ac.bed.gz",
        me="data/bmmc_dual/DNA/H3K27me3.bed.gz"
    message: "Mapping BMMC scNTT-seq"
    threads: 24
    shell:
        """
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            data/bmmc_dual/DNA/barcoded/scBM_k27me_S1_R1_001.barcoded.fastq.gz \
            data/bmmc_dual/DNA/barcoded/scBM_k27me_S1_R3_001.barcoded.fastq.gz \
            | samtools sort -@ {threads} -O bam - \
            > data/bmmc_dual/DNA/H3K27me3.bam
        
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            data/bmmc_dual/DNA/barcoded/scBM_k27ac_S2_R1_001.barcoded.fastq.gz \
            data/bmmc_dual/DNA/barcoded/scBM_k27ac_S2_R3_001.barcoded.fastq.gz \
            | samtools sort -@ {threads} -O bam - \
            > data/bmmc_dual/DNA/H3K27ac.bam
            
        samtools index data/bmmc_dual/DNA/H3K27ac.bam
        samtools index data/bmmc_dual/DNA/H3K27me3.bam
        
        sinto fragments -p {threads} \
          -b data/bmmc_dual/DNA/H3K27ac.bam \
          --barcode_regex "[^:]*" \
          -f data/bmmc_dual/DNA/H3K27ac.frag
          
        sinto fragments -p {threads} \
          -b data/bmmc_dual/DNA/H3K27me3.bam \
          --barcode_regex "[^:]*" \
          -f data/bmmc_dual/DNA/H3K27me3.frag

        sort -k1,1 -k2,2n data/bmmc_dual/DNA/H3K27ac.frag > data/bmmc_dual/DNA/H3K27ac.bed
        sort -k1,1 -k2,2n data/bmmc_dual/DNA/H3K27me3.frag > data/bmmc_dual/DNA/H3K27me3.bed

        rm data/bmmc_dual/DNA/H3K27ac.frag data/bmmc_dual/DNA/H3K27me3.frag
        bgzip -@ {threads} data/bmmc_dual/DNA/H3K27ac.bed
        bgzip -@ {threads} data/bmmc_dual/DNA/H3K27me3.bed
        tabix -p bed {output.ac}
        tabix -p bed {output.me}
        """

rule process_bmmc_atac:
    input: "data/bmmc_atac/download.done"
    output: "objects/bmmc_atac.rds"
    threads: 8
    message: "Generate BMMC scATAC Seurat object"
    shell:
        """
        Rscript code/bmmc_atac/process_bmmc_atac.r
        """
        
rule process_henikoff:
    input:
        "data/henikoff/GSM5034342_K27me3_R1_PBMC.fragments.HG38.tsv.gz",
        "data/henikoff/GSM5034343_K27me3_R2_PBMC.fragments.HG38.tsv.gz",
        "data/henikoff/GSM5034344_K27ac_PBMC.fragments.HG38.tsv.gz",
        "objects/pbmc_protein.rds"
    output: "objects/henikoff/ac_filt.rds", "objects/henikoff/me3_filt.rds"
    threads: 1
    message: "Processing Henikoff PBMC datasets"
    shell:
        """
        Rscript code/henikoff/process.R
        """

### Analysis ###

# cellculture

rule bulk_cellculture_covplot:
    input:
        "data/HEK_K562_bulk/mapped/K562-mono-K27ac.bw",
        "data/HEK_K562_bulk/mapped/K562-mono-K27me.bw",
        "data/HEK_K562_bulk/mapped/K562-mono-Pol2.bw",
        "data/HEK_K562_bulk/mapped/K562-plex-K27ac.bw",
        "data/HEK_K562_bulk/mapped/K562-plex-K27me.bw",
        "data/HEK_K562_bulk/mapped/K562-plex-Pol2.bw"
    output: "plots/hek_k562_bulk/bulk_covplot_k562.png"
    shell:
        """
        Rscript code/HEK_K562_bulk/analysis.R
        """

rule bulk_cellculture_heatmap:
    input:
        mono_ac="data/HEK_K562_bulk/mapped/K562-mono-K27ac.bw",
        mono_me3="data/HEK_K562_bulk/mapped/K562-mono-K27me.bw",
        mono_rna="data/HEK_K562_bulk/mapped/K562-mono-Pol2.bw",
        plex_ac="data/HEK_K562_bulk/mapped/K562-plex-K27ac.bw",
        plex_me3="data/HEK_K562_bulk/mapped/K562-plex-K27me.bw",
        plex_rna="data/HEK_K562_bulk/mapped/K562-plex-Pol2.bw",
        me3_peaks="data/k562_peaks/ENCFF031FSF.bed.gz",
        ac_peaks="data/k562_peaks/ENCFF038DDS.bed.gz",
        rna_peaks="data/k562_peaks/rna.bed"
    output:
        me3="plots/hek_k562_bulk/heatmap_k562_me3.png",
        ac="plots/hek_k562_bulk/heatmap_k562_ac.png",
        rna="plots/hek_k562_bulk/heatmap_k562_rna.png"
    threads: 24
    shell:
        """
        gzip -dc {input.me3_peaks} > data/k562_peaks/ENCFF031FSF.bed
        gzip -dc {input.ac_peaks} > data/k562_peaks/ENCFF038DDS.bed

        computeMatrix reference-point \
          --referencePoint 'center' \
          --missingDataAsZero \
          -S {input.plex_ac} {input.plex_me3} {input.plex_rna} \
          {input.mono_ac} {input.mono_me3} {input.mono_rna} \
          -R data/k562_peaks/ENCFF049HUP.bed \
          --upstream 5000 --downstream 5000 \
          -p {threads} \
          -o data/HEK_K562_bulk/h3k27me3.mat.gz

        computeMatrix reference-point \
          --referencePoint 'center' \
          --missingDataAsZero \
          -S {input.plex_ac} {input.plex_me3} {input.plex_rna} \
          {input.mono_ac} {input.mono_me3} {input.mono_rna} \
          -R data/k562_peaks/ENCFF038DDS.bed \
          --upstream 5000 --downstream 5000 \
          -p {threads} \
          -o data/HEK_K562_bulk/h3k27ac.mat.gz

        computeMatrix reference-point \
          --referencePoint 'center' \
          --missingDataAsZero \
          -S {input.plex_ac} {input.plex_me3} {input.plex_rna} \
          {input.mono_ac} {input.mono_me3} {input.mono_rna} \
          -R {input.rna_peaks} \
          --upstream 5000 --downstream 5000 \
          -p {threads} \
          -o data/HEK_K562_bulk/rna.mat.gz
          
        plotHeatmap \
          -m data/HEK_K562_bulk/h3k27me3.mat.gz  \
          -o {output.me3} \
          --whatToShow 'heatmap and colorbar' \
          --colorList '#FFFFFF,#F98401' '#FFFFFF,#D3145A' '#FFFFFF,#036C9A' '#FFFFFF,#676767' '#FFFFFF,#676767' '#FFFFFF,#676767'

        plotHeatmap \
          -m data/HEK_K562_bulk/h3k27ac.mat.gz  \
          -o {output.ac} \
          --zMax 0.5 \
          --whatToShow 'heatmap and colorbar' \
          --colorList '#FFFFFF,#F98401' '#FFFFFF,#D3145A' '#FFFFFF,#036C9A' '#FFFFFF,#676767' '#FFFFFF,#676767' '#FFFFFF,#676767'

        plotHeatmap \
          -m data/HEK_K562_bulk/rna.mat.gz  \
          -o {output.rna} \
          --zMax 0.5 \
          --whatToShow 'heatmap and colorbar' \
          --colorList '#FFFFFF,#F98401' '#FFFFFF,#D3145A' '#FFFFFF,#036C9A' '#FFFFFF,#676767' '#FFFFFF,#676767' '#FFFFFF,#676767'
        """

rule analyze_cellculture:
    input: "objects/hek_k562.rds", "data/k562_peaks/ENCFF031FSF.bed.gz"
    output: "plots/hek_k562/correlation.png"
    threads: 1
    shell:
        """
        Rscript code/HEK_K562_sc/analysis.R
        """

rule scatterplots_cellculture:
    input: "objects/hek_k562.rds"
    output: "plots/hek_k562/scatterplots_bulk.png"
    threads: 1
    shell:
        """
        Rscript code/HEK/K562_bulk/quantify_regions.R
        """
rule scatterplots_multict:
    input:
        "data/mesc/ENCFF008XKX.bed.gz",
        "data/mesc/ENCFF360VIS.bed.gz",
        "data/multict/fragments/H3K27ac-H3K27ac.tsv.gz"
    output:
        "plots/hek_k562/multi_cuttag_frip.png",
        "plots/hek_k562/multi_cuttag_scatter.png"
    threads: 1
    shell:
        """
        Rscript code/HEK_K562_sc/mct_scatter.R
        """

# pbmc/bmmc
rule process_pbmc_atac:
    input: "data/pbmc_atac/10k_pbmc_ATACv2_nextgem_Chromium_X_fragments.tsv.gz"
    output: "objects/pbmc_atac.rds"
    message: "Processing PBMC scATAC-seq"
    threads: 1
    shell:
        """
        Rscript code/pbmc_atac/process.R
        """

rule process_pbmc_ntt:
    input:
        adt="data/pbmc_protein/adt/outs/alevin/quants_mat.gz",
        ac="data/pbmc_protein/ntt/H3K27ac.tsv.gz",
        me="data/pbmc_protein/ntt/H3K27me3.tsv.gz",
    output: "objects/pbmc_protein.rds"
    threads: 10
    shell:
        """
        Rscript code/pbmc_protein/process.R {threads} {input.adt} {input.ac} {input.me} {output}
        """

rule process_pbmc_sc:
    input:
        ac="data/pbmc_sc/ntt/H3K27ac.tsv.gz",
        me="data/pbmc_sc/ntt/H3K27me3.tsv.gz",
    output: "objects/pbmc_sc.rds"
    threads: 10
    shell:
        """
        Rscript code/pbmc_sc/process.R {threads} {input.ac} {input.me} {output}
        """

rule split_pbmc_protein_fragments:
    input:
        obj="objects/pbmc_protein.rds",
        chrom="data/hg38.chrom.sizes"
    output:
        "data/pbmc_protein/bigwig/B_cell_ac.bw",
        "data/pbmc_protein/bigwig/pbmc_ac.bed.gz",
        "data/pbmc_protein/bigwig/pbmc_me3.bed.gz",
        "data/pbmc_protein/bigwig/pbmc_ac.bw",
        "data/pbmc_protein/bigwig/pbmc_me3.bw"
    message: "Creating PBMC bigwig files"
    threads: 1
    shell:
        """
        # split fragment files
        Rscript code/pbmc_protein/split.R
        
        # create bigwig
        cd data/pbmc_protein/bigwig
        for fragfile in $(ls -d *.bed.gz);do
          fname=(${{fragfile//.bed.gz/ }})
          bedtools genomecov -i $fragfile -g ../../../{input.chrom} -bg > "${{fname}}.bg"
          bedGraphToBigWig "${{fname}}.bg" ../../../{input.chrom} "${{fname}}.bw"
        done
        """

rule split_pbmc_sc_fragments:
    input:
        obj="objects/pbmc_sc.rds",
        chrom="data/hg38.chrom.sizes"
    output:
        "data/pbmc_sc/ntt/pbmc_ac.bed.gz",
        "data/pbmc_sc/ntt/pbmc_me3.bed.gz",
        "data/pbmc_sc/ntt/pbmc_ac.bw",
        "data/pbmc_sc/ntt/pbmc_me3.bw"
    message: "Creating PBMC filtered fragment files"
    threads: 1
    shell:
        """
        # split fragment files
        Rscript code/pbmc_sc/split.R
        
        # create bigwig
        cd data/pbmc_sc/ntt
        for fragfile in $(ls -d *.bed.gz);do
          fname=(${{fragfile//.bed.gz/ }})
          bedtools genomecov -i $fragfile -g ../../../{input.chrom} -bg > "${{fname}}.bg"
          bedGraphToBigWig "${{fname}}.bg" ../../../{input.chrom} "${{fname}}.bw"
        done
        """

rule split_bmmc_fragments:
    input: "objects/bmmc_dual.rds"
    output: "data/bmmc_dual/bulk_fragments/me3.bed.gz"
    message: "Splitting BMMC fragment file"
    threads: 1
    shell:
        """
        Rscript code/bmmc_dual/split.R
        """
        
rule pbmc_bulk_pseudobulk_cor:
    input:
        me3_sc="data/pbmc_protein/bigwig/pbmc_me3.bw",
        ac_sc="data/pbmc_protein/bigwig/pbmc_ac.bw",
        ac_chip="data/encode/ENCFF518PSI.bigWig",
        me3_chip="data/encode/ENCFF598EGZ.bigWig",
        me3_bulk_mono="data/pbmc_bulk/mapped/mono/me3.bw",
        me3_bulk_plex="data/pbmc_bulk/mapped/plex/me3.bw",
        ac_bulk_mono="data/pbmc_bulk/mapped/mono/ac.bw",
        ac_bulk_plex="data/pbmc_bulk/mapped/plex/ac.bw",
        all_peaks="data/encode/all_bulk.bed",
        ctp_ac_bw="data/ct_pro/H3K27ac.bw",
        ctp_me3_bw="data/ct_pro/H3K27me3.bw"
    output: "data/pbmc_bulk/encode_cor.tsv"
    threads: 6
    shell:
        """
        multiBigwigSummary BED-file -b \
          {input.me3_sc} \
          {input.ac_sc} \
          {input.me3_bulk_mono} \
          {input.me3_bulk_plex} \
          {input.ac_bulk_mono} \
          {input.ac_bulk_plex} \
          {input.me3_chip} \
          {input.ac_chip} \
          {input.ctp_ac_bw} \
          {input.ctp_me3_bw} \
          --BED {input.all_peaks} \
          -p {threads} \
          --outRawCounts {output} \
          -o data/pbmc_bulk/cor_matrix_all.npz
        """

rule pbmc_replicate_cor:
    input:
        me3_1="data/pbmc_protein/bigwig/pbmc_me3.bw",
        ac_1="data/pbmc_protein/bigwig/pbmc_ac.bw",
        me3_2="data/pbmc_sc/ntt/pbmc_me3.bw",
        ac_2="data/pbmc_sc/ntt/pbmc_ac.bw",
        all_peaks="data/encode/all_bulk.bed"
    output: "data/pbmc_protein/replicate_cor.tsv"
    threads: 6
    shell:
        """
        multiBigwigSummary BED-file -b \
          {input.me3_1} \
          {input.ac_1} \
          {input.me3_2} \
          {input.ac_2} \
          --BED {input.all_peaks} \
          -p {threads} \
          --outRawCounts {output} \
          -o data/pbmc_protein/replicate_cor_matrix_all.npz
        """

rule plot_bulk_cor_pbmc:
    input: "data/pbmc_bulk/encode_cor.tsv"
    output: "plots/pbmc/encode_cor_all.png"
    threads: 1
    shell:
        """
        Rscript code/pbmc_bulk/encode_cor_all.R
        """

rule pbmc_encode_cor:
    input:
        sc="data/pbmc_protein/bigwig/B_cell_ac.bw",
        encode="data/encode/ENCFF842JLZ.bigWig",
        me3="data/encode/h3k27me3.bed",
        ac="data/encode/h3k27ac.bed"
    output:
        "data/pbmc_protein/bigwig/raw_ac.tsv",
        "data/pbmc_protein/bigwig/raw_me3.tsv",
        "data/pbmc_protein/bigwig/raw.tsv"
    message: "Computing ENCODE/NTT correlations"
    threads: 6
    shell:
        """
        # ac cor
        multiBigwigSummary BED-file -b \
          data/pbmc_protein/bigwig/B_cell_ac.bw \
          data/pbmc_protein/bigwig/NK_ac.bw \
          data/pbmc_protein/bigwig/CD14_Mono_ac.bw \
          data/pbmc_protein/bigwig/CD8_T_cell_ac.bw \
          data/pbmc_protein/bigwig/CD4_T_cell_ac.bw \
          data/pbmc_protein/bigwig/Late_erythroid_ac.bw \
          data/encode/ENCFF293ETP.bigWig \
          data/encode/ENCFF611GRL.bigWig \
          data/encode/ENCFF181OXO.bigWig \
          data/encode/ENCFF526VJO.bigWig \
          data/encode/ENCFF601NLG.bigWig \
          data/encode/ENCFF951ZBV.bigWig \
          data/encode/ENCFF064JOI.bigWig \
          --BED {input.ac} \
          -p {threads} \
          --outRawCounts data/pbmc_protein/bigwig/raw_ac.tsv \
          -o data/pbmc_protein/bigwig/cor_matrix_ac.npz
        
        # me3 cor
        multiBigwigSummary BED-file -b \
          data/pbmc_protein/bigwig/B_cell_me3.bw \
          data/pbmc_protein/bigwig/NK_me3.bw \
          data/pbmc_protein/bigwig/CD14_Mono_me3.bw \
          data/pbmc_protein/bigwig/CD8_T_cell_me3.bw \
          data/pbmc_protein/bigwig/CD4_T_cell_me3.bw \
          data/pbmc_protein/bigwig/Late_erythroid_me3.bw \
          data/encode/ENCFF569IPX.bigWig \
          data/encode/ENCFF842JLZ.bigWig \
          data/encode/ENCFF811VAX.bigWig \
          data/encode/ENCFF499VWN.bigWig \
          data/encode/ENCFF777EGG.bigWig \
          data/encode/ENCFF046VLL.bigWig \
          data/encode/ENCFF085YMB.bigWig \
          --BED {input.me3} \
          -p {threads} \
          --outRawCounts data/pbmc_protein/bigwig/raw_me3.tsv \
          -o data/pbmc_protein/bigwig/cor_matrix_me3.npz
        
        # all cor
        multiBigwigSummary BED-file -b \
          data/pbmc_protein/bigwig/B_cell_me3.bw \
          data/pbmc_protein/bigwig/NK_me3.bw \
          data/pbmc_protein/bigwig/CD14_Mono_me3.bw \
          data/pbmc_protein/bigwig/CD8_T_cell_me3.bw \
          data/pbmc_protein/bigwig/CD4_T_cell_me3.bw \
          data/pbmc_protein/bigwig/Late_erythroid_me3.bw \
          data/pbmc_protein/bigwig/B_cell_ac.bw \
          data/pbmc_protein/bigwig/NK_ac.bw \
          data/pbmc_protein/bigwig/CD14_Mono_ac.bw \
          data/pbmc_protein/bigwig/CD8_T_cell_ac.bw \
          data/pbmc_protein/bigwig/CD4_T_cell_ac.bw \
          data/pbmc_protein/bigwig/Late_erythroid_ac.bw \
          data/encode/ENCFF569IPX.bigWig \
          data/encode/ENCFF842JLZ.bigWig \
          data/encode/ENCFF811VAX.bigWig \
          data/encode/ENCFF499VWN.bigWig \
          data/encode/ENCFF777EGG.bigWig \
          data/encode/ENCFF046VLL.bigWig \
          data/encode/ENCFF085YMB.bigWig \
          data/encode/ENCFF293ETP.bigWig \
          data/encode/ENCFF611GRL.bigWig \
          data/encode/ENCFF181OXO.bigWig \
          data/encode/ENCFF526VJO.bigWig \
          data/encode/ENCFF601NLG.bigWig \
          data/encode/ENCFF951ZBV.bigWig \
          data/encode/ENCFF064JOI.bigWig \
          --BED data/encode/all.bed \
          -p {threads} \
          --outRawCounts data/pbmc_protein/bigwig/raw.tsv \
          -o data/pbmc_protein/bigwig/cor_matrix_all.npz
        """

rule plot_pbmc_encode_cor:
    input:
        "data/pbmc_protein/bigwig/raw_ac.tsv",
        "data/pbmc_protein/bigwig/raw_me3.tsv",
        "data/pbmc_protein/bigwig/raw.tsv",
        "data/pbmc_protein/replicate_cor.tsv"
    output:
        "plots/pbmc/encode_cor_spearman.png",
        "plots/pbmc/replicate_correlation.png"
    threads: 1
    shell:
        """
        Rscript code/pbmc_protein/encode_cor.R
        """

rule bulk_pbmc_covplot:
    input:
        "data/pbmc_bulk/mapped/mono/me3.bw",
        "data/pbmc_bulk/mapped/mono/ac.bw",
        "data/pbmc_bulk/mapped/plex/me3.bw",
        "data/pbmc_bulk/mapped/plex/ac.bw"
    output: "plots/pbmc_bulk/bulk_covplot_pbmc.png"
    shell:
        """
        Rscript code/pbmc_bulk/analysis.R
        """

rule pbmc_bulk_frip:
    input: 
        "data/pbmc_bulk/mapped/mono/me3.bed.gz",
        "data/pbmc_bulk/mapped/mono/ac.bed.gz",
        "data/pbmc_bulk/mapped/plex/me3.bed.gz",
        "data/pbmc_bulk/mapped/plex/ac.bed.gz",
        "data/encode/ENCFF832RWT.bed.gz",
        "data/encode/ENCFF291LVP.bed.gz"
    output: "plots/pbmc_bulk/frip_ac.png", "plots/pbmc_bulk/frip_me3.png"
    shell:
        """
        Rscript code/pbmc_bulk/frip.R
        """
        
rule pbmc_sc_frip:
    input: 
        "objects/pbmc_protein.rds",
        "objects/pbmc_sc.rds",
        "data/ct_pro/H3K27ac_updated.rds",
        "data/ct_pro/H3K27me3_updated.rds",
        "objects/henikoff/me3_filt.rds",
        "objects/henikoff/ac_filt.rds",
        "data/encode/ENCFF832RWT.bed.gz",
        "data/encode/ENCFF291LVP.bed.gz"
    output: "plots/pbmc/frip_pbmc.png"
    shell:
        """
        Rscript code/pbmc_protein/frip.R
        """

rule bmmc_frip:
    input:
        "objects/bmmc_dual.rds",
        "data/encode/ENCFF832RWT.bed.gz",
        "data/encode/ENCFF291LVP.bed.gz"
    output: "plots/bmmc/frip_bmmc.png"
    shell:
        """
        Rscript code/bmmc_dual/frip_bmmc.R
        """
        
rule bulk_pbmc_heatmap:
    input:
        mono_ac="data/pbmc_bulk/mapped/mono/ac.bw",
        mono_me3="data/pbmc_bulk/mapped/mono/me3.bw",
        plex_ac="data/pbmc_bulk/mapped/plex/ac.bw",
        plex_me3="data/pbmc_bulk/mapped/plex/me3.bw",
        me3_peaks="data/encode/ENCFF291LVP.bed.gz",
        ac_peaks="data/encode/ENCFF832RWT.bed.gz"
    output:
        me3="plots/pbmc_bulk/heatmap_pbmc_me3.png",
        ac="plots/pbmc_bulk/heatmap_pbmc_ac.png"
    threads: 24
    shell:
        """
        gzip -dc {input.me3_peaks} > data/encode/ENCFF291LVP.bed
        gzip -dc {input.ac_peaks} > data/encode/ENCFF832RWT.bed

        computeMatrix reference-point \
          --referencePoint 'center' \
          --missingDataAsZero \
          -S {input.plex_ac} {input.plex_me3} {input.mono_ac} {input.mono_me3} \
          -R data/encode/ENCFF291LVP.bed \
          --upstream 5000 --downstream 5000 \
          -p {threads} \
          -o data/pbmc_bulk/h3k27me3.mat.gz

        computeMatrix reference-point \
          --referencePoint 'center' \
          --missingDataAsZero \
          -S {input.plex_ac} {input.plex_me3} {input.mono_ac} {input.mono_me3} \
          -R data/encode/ENCFF832RWT.bed \
          --upstream 5000 --downstream 5000 \
          -p {threads} \
          -o data/pbmc_bulk/h3k27ac.mat.gz
        
        plotHeatmap \
          -m data/pbmc_bulk/h3k27me3.mat.gz \
          -o {output.me3} \
          --whatToShow 'heatmap only' \
          --colorList '#FFFFFF,#F98401' '#FFFFFF,#D3145A' '#FFFFFF,#676767' '#FFFFFF,#676767'

        plotHeatmap \
          -m data/pbmc_bulk/h3k27ac.mat.gz \
          -o {output.ac} \
          --whatToShow 'heatmap only' \
          --colorList '#FFFFFF,#F98401' '#FFFFFF,#D3145A' '#FFFFFF,#676767' '#FFFFFF,#676767'
        """

rule pbmc_bmmc_scatter:
    input:
        pbmc_protein_me3="data/pbmc_protein/bigwig/pbmc_me3.bed.gz",
        pbmc_protein_ac="data/pbmc_protein/bigwig/pbmc_ac.bed.gz",
        pbmc_bulk_me3="data/pbmc_bulk/mapped/mono/me3.bed.gz",
        pbmc_bulk_ac="data/pbmc_bulk/mapped/mono/ac.bed.gz",
        bmmc_me3="data/bmmc_dual/bulk_fragments/me3.bed.gz",
        bmmc_ac="data/bmmc_dual/bulk_fragments/ac.bed.gz",
        ac_peaks="data/encode/ENCFF832RWT.bed.gz",
        me3_peaks="data/encode/ENCFF291LVP.bed.gz"
    output: "plots/pbmc/scatterplot_pbmc.png", "plots/bmmc/scatterplot_bmmc.png"
    threads: 1
    shell:
        """
        Rscript code/pbmc_protein/scatterplots.R {input.ac_peaks} {input.me3_peaks}
        """

rule pbmc_sc_scatter:
    input:
        pbmc_sc_me3="data/pbmc_sc/ntt/pbmc_me3.bed.gz",
        pbmc_sc_ac="data/pbmc_sc/ntt/pbmc_ac.bed.gz",
        pbmc_bulk_me3="data/pbmc_bulk/mapped/mono/me3.bed.gz",
        pbmc_bulk_ac="data/pbmc_bulk/mapped/mono/ac.bed.gz",
        ac_peaks="data/encode/ENCFF832RWT.bed.gz",
        me3_peaks="data/encode/ENCFF291LVP.bed.gz"
    output: "plots/pbmc/scatterplot_pbmc_sc_ac_me3.png", "plots/pbmc/scatterplot_pbmc_sc_ac_ac.png"
    threads: 1
    shell:
        """
        Rscript code/pbmc_sc/scatterplots_sc.R {input.ac_peaks} {input.me3_peaks}
        """

rule process_bmmc_dual:
    input:
        ac="data/bmmc_dual/DNA/H3K27ac.bed.gz",
        me="data/bmmc_dual/DNA/H3K27me3.bed.gz"
    output: "objects/bmmc_dual.rds"
    threads: 10
    shell:
        """
        Rscript code/bmmc_dual/process.R {threads} {output}
        """

rule bmmc_pseudotime:
    input: "objects/bmmc_dual.rds", "objects/fullref.Rds"
    output:
        "plots/bmmc/hmap_ac.png",
        "plots/bmmc/trajectory.png",
        "plots/bmmc/rna_repressed.pdf",
        "plots/bmmc/rna_activated.pdf"
    threads: 12
    shell:
        """
        Rscript code/bmmc_dual/trajectory.R
        """

rule transfer_bmmc:
    input: "objects/bmmc_atac.rds", "objects/bmmc_dual.rds"
    output: "data/bmmc_dual_predictions.tsv"
    shell:
        """
        Rscript code/bmmc_dual/transfer.R
        """

rule analyze_bmmc_dual:
    input: "objects/bmmc_dual.rds"
    output: "plots/bmmc/fragments_bmmc.png"
    shell:
        """
        Rscript code/bmmc_dual/analysis.R
        """

rule analyse_pbmc:
    input: 
        ntt="objects/pbmc_protein.rds",
        ct_ac="data/ct_pro/H3K27ac_updated.rds",
        ct_me3="data/ct_pro/H3K27me3_updated.rds",
        ntt_2="data/pbmc_sc.rds"
    output: "plots/pbmc/fragments.png"
    message: "Processing PBMC"
    threads: 1
    shell:
        """
        Rscript code/pbmc_protein/analysis.R {input.ntt} {input.ct_ac} {input.ct_me3} {input.ntt_2}
        """

rule pbmc_bulk:
    input:
        me3_peaks="data/encode/ENCFF291LVP.bed.gz",
        ac_peaks="data/encode/ENCFF832RWT.bed.gz",
        me3_mono="data/pbmc_bulk/mapped/mono/me3.bed.gz",
        me3_plex="data/pbmc_bulk/mapped/plex/me3.bed.gz",
        ac_mono="data/pbmc_bulk/mapped/mono/ac.bed.gz",
        ac_plex="data/pbmc_bulk/mapped/plex/ac.bed.gz"
    output: "plots/pbmc/bulk_scatter.png"
    threads: 1
    shell:
        """
        Rscript code/pbmc_bulk/quantify_regions.R
        """
