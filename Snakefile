
rule all:
    input:
        "plots/pbmc/bulk_scatter.png",
        "plots/pbmc/fragments.png",
        "plots/bmmc/fragments_bmmc.png",
        "plots/bmmc/hmap_ac.png",
        "plots/bmmc/trajectory.png",
        "plots/bmmc/rna_repressed.pdf",
        "plots/hek_k562/scatterplots_bulk.png"
        "plots/hek_k562/correlation.png"

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
        "genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.chromap",
        "genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        "genome/hg38/genome.fa"
    threads: 1
    message: "Download hg38 genome"
    shell:
        """
        cd genome/hg38
        aws s3 sync s3://stuart-genomes/hg38_analysis/bwa-mem2/ .
        aws s3 cp s3://stuart-genomes/hg38_analysis/chromap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.chromap .
        aws s3 cp s3://stuart-genomes/hg38_analysis/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz .
        
        # decompress and index fasta for vartrix
        gzip -dc GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > genome.fa
        samtools faidx genome.fa
        picard CreateSequenceDictionary -R genome.fa -O genome.dict
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

### Map ###

rule map_cellculture_sc:
    input:
        r1="data/HEK_K562_sc/dmx/scCmix_{mark}_R1.barcoded.fastq.gz",
        r2="data/HEK_K562_sc/dmx/scCmix_{mark}_R3.barcoded.fastq.gz",
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output:
        frags="data/HEK_K562_sc/sinto/{mark}.bed.gz",
        bam="data/HEK_K562_sc/mapped/{mark}.bam"
    threads: 18
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
        index="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.chromap",
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output: "data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bed.gz"
    threads: 6
    message: "Mapping {input.r1}"
    shell:
        """
        chromap --preset atac \
          -x {input.index} \
          -r {input.genome} \
          -1 {input.r1} \
          -2 {input.r2} \
          -o data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.tmp \
          -t {threads}
        sort -k1,1 -k2,2n data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.tmp \
          > data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.bed
        bgzip -@ {threads} data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.bed
        tabix -p bed {output}
        rm data/HEK_K562_bulk/mapped/{wildcards.cell}-{wildcards.plex}-{wildcards.mark}.tmp
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
        index="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.chromap",
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output:
        me3_mono="data/pbmc_bulk/mapped/mono/me3.bed.gz",
        me3_plex="data/pbmc_bulk/mapped/plex/me3.bed.gz",
        ac_mono="data/pbmc_bulk/mapped/mono/ac.bed.gz",
        ac_plex="data/pbmc_bulk/mapped/plex/ac.bed.gz"
    threads: 12
    message: "Mapping PBMC bulk NTT-seq"
    shell:
        """
        # me3
        chromap --preset atac \
          -x {input.index} \
          -r {input.genome} \
          -1 {input.me3_r1_mono} \
          -2 {input.me3_r2_mono} \
          -o data/pbmc_bulk/mapped/mono/me3.tmp \
          -t {threads}

        sort -k1,1 -k2,2n data/pbmc_bulk/mapped/mono/me3.tmp \
          > data/pbmc_bulk/mapped/mono/me3.bed
        bgzip -@ {threads} data/pbmc_bulk/mapped/mono/me3.bed
        tabix -p bed {output.me3_mono}
        rm data/pbmc_bulk/mapped/mono/me3.tmp
        
        chromap --preset atac \
          -x {input.index} \
          -r {input.genome} \
          -1 {input.me3_r1_plex} \
          -2 {input.me3_r2_plex} \
          -o data/pbmc_bulk/mapped/plex/me3.tmp \
          -t {threads}

        sort -k1,1 -k2,2n data/pbmc_bulk/mapped/plex/me3.tmp \
          > data/pbmc_bulk/mapped/plex/me3.bed
        bgzip -@ {threads} data/pbmc_bulk/mapped/plex/me3.bed
        tabix -p bed {output.me3_plex}
        rm data/pbmc_bulk/mapped/plex/me3.tmp
        
        chromap --preset atac \
          -x {input.index} \
          -r {input.genome} \
          -1 {input.ac_r1_mono} \
          -2 {input.ac_r2_mono} \
          -o data/pbmc_bulk/mapped/mono/ac.tmp \
          -t {threads}

        sort -k1,1 -k2,2n data/pbmc_bulk/mapped/mono/ac.tmp \
          > data/pbmc_bulk/mapped/mono/ac.bed
        bgzip -@ {threads} data/pbmc_bulk/mapped/mono/ac.bed
        tabix -p bed {output.ac_mono}
        rm data/pbmc_bulk/mapped/mono/ac.tmp
        
        chromap --preset atac \
          -x {input.index} \
          -r {input.genome} \
          -1 {input.ac_r1_plex} \
          -2 {input.ac_r2_plex} \
          -o data/pbmc_bulk/mapped/plex/ac.tmp \
          -t {threads}

        sort -k1,1 -k2,2n data/pbmc_bulk/mapped/plex/ac.tmp \
          > data/pbmc_bulk/mapped/plex/ac.bed
        bgzip -@ {threads} data/pbmc_bulk/mapped/plex/ac.bed
        tabix -p bed {output.ac_plex}
        rm data/pbmc_bulk/mapped/plex/ac.tmp
        """

rule create_bulk_pbmc_bigwig:
    input:
        fragments="data/pbmc_bulk/mapped/{plex}/{mark}.bed.gz",
        chrom="data/hg38.chrom.sizes"
    output:
        bg="data/pbmc_bulk/mapped/{plex}/{mark}.bg",
        bw="data/pbmc_bulk/mapped/{plex}/{mark}.bw"
    threads: 1
    shell:
        """
        bedtools genomecov -i {input.fragments} -g {input.chrom} -bg > {output.bg}
        bedGraphToBigWig {output.bg} {input.chrom} {output.bw}
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
    threads: 1
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
    output: "data/k562_peaks/ENCFF031FSF.bed.gz"
    threads: 1
    shell:
        """
        wget -i {input} -P data/k562_peaks
        """

rule cellculture_bulk_bw:
    input:
        bed="data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bed.gz",
        chrom="data/hg38.chrom.sizes"
    output:
        bg="data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bedgraph",
        bw="data/HEK_K562_bulk/mapped/{cell}-{plex}-{mark}.bw",
    threads: 6
    message: "Creating bigwig for {input}"
    shell:
        """
        bedtools genomecov -i {input.bed} -g {input.chrom} -bg > {output.bg}
        bedGraphToBigWig {output.bg} {input.chrom} {output.bw}
        """

rule demux_pbmc_ntt:
    input:
        read1="data/pbmc_protein/ntt/Undetermined_S0_R1_001.fastq.gz",
        readi5="data/pbmc_protein/ntt/Undetermined_S0_R2_001.fastq.gz",
        read2="data/pbmc_protein/ntt/Undetermined_S0_R3_001.fastq.gz",
        bc="data/pbmc_protein/barcodes.fa"
    output: directory("data/pbmc_protein/ntt/out")
    message: "Demultiplex NTT reads"
    threads: 1
    shell:
        """
        python code/pbmc_protein/demux.py \
           --read1 {input.read1} \
           --read_i5 {input.readi5} \
           --read2 {input.read2} \
           --tn5 {input.bc} \
           --output {output}
        """

rule map_pbmc_ntt:
    input:
        genome="genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        reads=directory("data/pbmc_protein/ntt/out")
    output:
        ac="data/pbmc_protein/ntt/H3K27ac.tsv.gz",
        me="data/pbmc_protein/ntt/H3K27me3.tsv.gz"
    message: "Mapping PBMC scNTT-seq"
    threads: 24
    shell:
        """
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            data/pbmc_protein/ntt/out/H3K27ac.R1.fastq \
            data/pbmc_protein/ntt/out/H3K27ac.R2.fastq \
            | samtools sort -@ {threads} -O bam - \
            > data/pbmc_protein/ntt/out/H3K27ac.bam
        
        bwa-mem2 mem {input.genome} \
            -t {threads} \
            data/pbmc_protein/ntt/out/H3K27me3.R1.fastq \
            data/pbmc_protein/ntt/out/H3K27me3.R2.fastq \
            | samtools sort -@ {threads} -O bam - \
            > data/pbmc_protein/ntt/out/H3K27me3.bam
            
        samtools index data/pbmc_protein/ntt/out/H3K27ac.bam
        samtools index data/pbmc_protein/ntt/out/H3K27me3.bam
        
        sinto fragments -p {threads} \
          -b data/pbmc_protein/ntt/out/H3K27ac.bam \
          --barcode_regex "[^:]*" \
          -f data/pbmc_protein/ntt/H3K27ac.frag
          
        sinto fragments -p {threads} \
          -b data/pbmc_protein/ntt/out/H3K27me3.bam \
          --barcode_regex "[^:]*" \
          -f data/pbmc_protein/ntt/H3K27me3.frag

        sort -k1,1 -k2,2n data/pbmc_protein/ntt/H3K27ac.frag > data/pbmc_protein/ntt/H3K27ac.tsv
        sort -k1,1 -k2,2n data/pbmc_protein/ntt/H3K27me3.frag > data/pbmc_protein/ntt/H3K27me3.tsv

        rm data/pbmc_protein/ntt/H3K27ac.frag data/pbmc_protein/ntt/H3K27me3.frag
        bgzip -@ {threads} data/pbmc_protein/ntt/H3K27ac.tsv
        bgzip -@ {threads} data/pbmc_protein/ntt/H3K27me3.tsv
        tabix -p bed {output.ac}
        tabix -p bed {output.me}
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
    threads: 1
    message: "Combining ADT reads for quantification"
    shell:
        """
        paste <(gunzip -c {input.r1}) <(gunzip -c {input.r3}) \
          | paste - - - - \
          | awk -F'\\t' -v one="data/pbmc_protein/adt/read1.fastq" '{{ OFS="\\n"; print $1,$3$4,$5,$7$8 >> one;}}'
        
        gzip data/pbmc_protein/adt/read1.fastq
        """

rule map_pbmc_adt:
    input:
        index=directory("data/adt_index"),
        reads="data/pbmc_protein/adt/read1.fastq.gz"
    output: "data/pbmc_protein/adt/outs/alevin/quants_mat.gz"
    message: "Quantify PBMC ADTs"
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
          -b 16 \
          
        sinto barcode \
          --barcode_fastq K27me_R2.fastq \
          --read1 scBM_k27me_S1_R1_001.fastq \
          --read2 scBM_k27me_S1_R3_001.fastq \
          -b 16 \
          
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

### Analysis ###

# cellculture

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
        Rscript code/HEK_K562_sc/process.R {threads} {output.obj} "data/HEK_K562_sc/sinto
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

# pbmc/bmmc
rule process_pbmc_ntt:
    input:
        adt="data/pbmc_protein/adt/outs/alevin/quants_mat.gz",
        ac="data/pbmc_protein/ntt/H3K27ac.tsv.gz",
        me="data/pbmc_protein/ntt/H3K27me3.tsv.gz"
    output: "objects/pbmc_protein.rds"
    threads: 10
    shell:
        """
        Rscript code/pbmc_protein/process.R {threads} {input.adt} {input.ac} {input.me} {output}
        """
        
rule split_pbmc_fragments:
    input:
        obj="objects/pbmc_protein.rds",
        chrom="data/hg38.chrom.sizes"
    output: "data/pbmc_protein/bigwig/B_cell_ac.bw"
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
        ac_peaks="data/encode/ENCFF832RWT.bed.gz",
        me3_peaks="data/encode/ENCFF291LVP.bed.gz",
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

rule pbmc_encode_cor:
    input:
        sc="data/pbmc_protein/bigwig/B_cell_ac.bw",
        encode="data/encode/ENCFF842JLZ.bigWig",
        me3="data/encode/h3k27me3.bed",
        ac="data/encode/h3k27ac.bed"
    output: "data/pbmc_protein/bigwig/raw_ac.tsv", "data/pbmc_protein/bigwig/raw_me3.tsv", "data/pbmc_protein/bigwig/raw.tsv"
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
          data/pbmc_protein/bigwig/HSC_ac.bw \
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
          data/pbmc_protein/bigwig/HSC_me3.bw \
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
          data/pbmc_protein/bigwig/HSC_me3.bw \
          data/pbmc_protein/bigwig/B_cell_ac.bw \
          data/pbmc_protein/bigwig/NK_ac.bw \
          data/pbmc_protein/bigwig/CD14_Mono_ac.bw \
          data/pbmc_protein/bigwig/CD8_T_cell_ac.bw \
          data/pbmc_protein/bigwig/CD4_T_cell_ac.bw \
          data/pbmc_protein/bigwig/HSC_ac.bw \
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
        atac="objects/pbmc_atac.rds",
        ct_ac="data/ct_pro/H3K27ac_updated.rds",
        ct_me3="data/ct_pro/H3K27me3_updated.rds"
    output: "plots/pbmc/fragments.png"
    message: "Processing PBMC"
    threads: 1
    shell:
        """
        Rscript code/pbmc_protein/analysis.R {input.ntt} {input.atac} {input.ct_ac} {input.ct_me3}
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
