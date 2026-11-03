# 0. Retrieve the incoming parameters
    while getopts "i:o:t:" opt
    do
        case $opt in
            i)
            indir=$OPTARG # 原始测序数据的绝对路径
            ;;
            o)
            outdir=$OPTARG # 结果输出路径
            ;;
            t)
            threads=$OPTARG # 同时并行跑几个样本
            ;;
            ?)
            echo "Usage: `basename $0` [-i indir] [-o outdir] [-t threads]"
            exit 1
            ;;
        esac
    done

# 1. fastqc and multiqc
    mkdir -p ${outdir}
    shopt -s globstar
    fastqc ${indir}/**/*fq* -t ${threads} -o ${outdir} > /dev/null
    multiqc ${outdir}/* -o ${outdir}
    seqkit stats ${indir}/**/*fq* > ${outdir}/Seqkit_stats.txt

# 2. Extract key information
    workdir=`dirname ${outdir}`
    # 输出质量分数体系
    awk -F '\t' '{print $4}' ${outdir}/multiqc_data/multiqc_fastqc.txt|sed '1d'|sort|uniq > ${workdir}/Encoding     # Sanger / Illumina 1.9为phred33，Illumina 1.5为phred64
    
    # 输出没去除干净接头的run列表
    jq -r '.report_plot_data[].datasets[].[].name' ${outdir}/multiqc_data/multiqc_data.json|awk -F ' - ' '{print $1"\t"$2}'|sort -k 1 > ${workdir}/Adapter.1
    awk -F '\t' '{print $1"\t"$NF}' ${outdir}/multiqc_data/multiqc_fastqc.txt > ${workdir}/Adapter.2
    cat ${workdir}/Adapter.1|while read line;
    do
        file=`echo ${line}|awk '{print $1}'`
        grep ${file} ${workdir}/Adapter.2|sort -k 1 >> ${workdir}/Adapter.3
    done
    paste ${workdir}/Adapter.1 ${workdir}/Adapter.3|awk '{print $1"\t"$2"\t"$4}' > ${workdir}/Adapter
    rm ${workdir}/Adapter.1 ${workdir}/Adapter.2 ${workdir}/Adapter.3
    
    cd ${workdir}
    echo "FastQC:OK." >> Analyze_progress

# 3. Import environment variables / Activate QIIME2 environment
    conda activate qiime2-amplicon-2023.9

# 4. Set up the table of contents and prepare the manifest file
    indir='/mnt/raid9/01_research/04_family_transmission/03_Feces_16S_Bacteria_MbPL2021121899/split'
    workdir='/mnt/raid9/01_research/04_family_transmission/c_Feces_16S'
    ls ${indir}|cut -d "_" -f 1|sort|uniq > ${workdir}/sample.list
    bash ${workdir}/00_script/01_prepare_manifest.sh ${workdir}/sample.list ${indir}

# 5. Quality control, checking the sequencing data quality system
    bash ${workdir}/00_script/00_fastqc.sh -i ${indir} -o ${workdir}/02_rawdata_fastqc/ -t 10
    cat ${workdir}/Encoding # Sanger / Illumina 1.9为phred33，Illumina 1.5为phred64

# 6. Import data
    mkdir -p ${workdir}/03_qiime
    
    # 测序是单端
        qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path ${workdir}/manifest_SE.txt \
        --output-path ${workdir}/03_qiime/01_single-end.qza \
        --input-format SingleEndFastqManifestPhred33V2 # 如果上面Encoding为Illumina 1.5，这里需要改成PairedEndFastqManifestPhred64V2
    
    # 测序是双端
        qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${workdir}/manifest_PE.txt \
        --output-path ${workdir}/03_qiime/01_paired-end.qza \
        --input-format PairedEndFastqManifestPhred33V2 # 如果上面Encoding为Illumina 1.5，这里需要改成PairedEndFastqManifestPhred64V2

# 7. Visualization
    # 单端
        qiime demux summarize  --i-data ${workdir}/03_qiime/01_single-end.qza --o-visualization ${workdir}/03_qiime/01_single-end.qzv
    # 双端
        qiime demux summarize  --i-data ${workdir}/03_qiime/01_paired-end.qza --o-visualization ${workdir}/03_qiime/01_paired-end.qzv

# 8. DADA2

    # 单端
        # 这里的截断参数需要进行记录，之后要上传至数据库
        # --p-trim-left：这里要把上一步可视化生成的02_demux.qzv文件下载到本地电脑，再上传到https://view.qiime2.org/进行可视化，根据序列的质量进行低质量序列的截短，一般设置在0~6之间
        # --p-trunc-len：这里要把上一步可视化生成的02_demux.qzv文件下载到本地电脑，再上传到https://view.qiime2.org/进行可视化，根据序列的质量进行低质量序列的截短，这里的长度需要参考16S测序的区段长度，保留的序列长度应该比测序区段的长度长至少20bp
        qiime dada2 denoise-single \
        --i-demultiplexed-seqs ${workdir}/03_qiime/01_single-end.qza \
        --p-trim-left 0 \
        --p-trunc-len 300 \
        --o-table ${workdir}/03_qiime/02_dada2_feat_table.qza \
        --o-representative-sequences ${workdir}/03_qiime/02_dada2_rep_seqs.qza \
        --o-denoising-stats ${workdir}/03_qiime/02_dada2_denoising_stats.qza \
        --p-n-threads 20  # 同时占多少个cpu，如果设置成20个，则watson最多只能同时跑两个这种任务
    
    # 双端
        # 这里的四个截断参数需要进行记录，之后要上传至数据库
        # --p-trim-left-f、--p-trim-left-r：这里要把上一步可视化生成的02_demux.qzv文件下载到本地电脑，再上传到https://view.qiime2.org/进行可视化，根据序列的质量进行低质量序列的截短，一般设置在0~6之间
        # --p-trunc-len-f、--p-trunc-len-r：这里要把上一步可视化生成的02_demux.qzv文件下载到本地电脑，再上传到https://view.qiime2.org/进行可视化，根据序列的质量进行低质量序列的截短，这里的长度需要参考16S测序的区段长度，双端保留的序列长度加起来应该比测序区段的长度长至少20bp
        qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${workdir}/03_qiime/01_paired-end.qza \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-trunc-len-f 240 \
        --p-trunc-len-r 240 \
        --o-table ${workdir}/03_qiime/02_dada2_feat_table.qza \
        --o-representative-sequences ${workdir}/03_qiime/02_dada2_rep_seqs.qza \
        --o-denoising-stats ${workdir}/03_qiime/02_dada2_denoising_stats.qza \
        --p-n-threads 40

# 9. Deblur

    # 单端数据
        # Q-quality fiter
        qiime quality-filter q-score \
            --i-demux ${workdir}/03_qiime/01_single-end.qza \
            --o-filtered-sequences ${workdir}/03_qiime/01_quality-filtered.qza \
            --o-filter-stats ${workdir}/03_qiime/01_quality-filter-stats.qza
        qiime demux summarize  --i-data ${workdir}/03_qiime/01_quality-filtered.qza --o-visualization ${workdir}/03_qiime/01_quality-filtered.qzv
    
        # Deblur
        # the Deblur developers recommend setting this value to a length where the median quality score begins to drop too low, 要看01_quality-filtered.qzv可视化的结果
        qiime deblur denoise-16S \
            --i-demultiplexed-seqs ${workdir}/03_qiime/01_quality-filtered.qza \
            --p-trim-length 120 \
            --o-representative-sequences ${workdir}/03_qiime/02_deblur_rep_seqs.qza \
            --o-table ${workdir}/03_qiime/02_deblur_feat_table.qza \
            --p-sample-stats \
            --o-stats ${workdir}/03_qiime/02_deblur_denoising_stats.qza


    # 双端数据
        # Join paired reads
        qiime vsearch merge-pairs \
            --i-demultiplexed-seqs ${workdir}/03_qiime/01_paired-end.qza \
            --o-merged-sequences ${workdir}/03_qiime/01_paired-joined.qza
        qiime demux summarize \
            --i-data ${workdir}/03_qiime/01_paired-joined.qza \
            --o-visualization ${workdir}/03_qiime/01_paired-joined.qzv
    
        # Q-quality fiter
        qiime quality-filter q-score \
            --i-demux ${workdir}/03_qiime/01_paired-joined.qza \
            --o-filtered-sequences ${workdir}/03_qiime/01_quality-filtered.qza \
            --o-filter-stats ${workdir}/03_qiime/01_quality-filter-stats.qza
        qiime demux summarize \
            --i-data ${workdir}/03_qiime/01_quality-filtered.qza \
            --o-visualization ${workdir}/03_qiime/01_quality-filtered.qzv
    
        # Deblur
        # the Deblur developers recommend setting this value to a length where the median quality score begins to drop too low, 要看03_quality-filtered.qzv可视化的结果
        qiime deblur denoise-16S \
            --i-demultiplexed-seqs ${workdir}/03_qiime/01_quality-filtered.qza \
            --p-trim-length 460 \
            --p-sample-stats \
            --o-representative-sequences ${workdir}/03_qiime/02_deblur_rep_seqs.qza \
            --o-table ${workdir}/03_qiime/02_deblur_feat_table.qza \
            --o-stats ${workdir}/03_qiime/02_deblur_denoising_stats.qza


# 10. Taxonomy Annotation
    # 有v4测序区域用的是silva的16S数据构建的分类器515F/806R，无v4测序区域建议使用全长参考数据库
    denoising_method='dada2'  # 或者deblur
    # v4
        qiime feature-classifier classify-sklearn \
        --i-classifier /mnt/raid8/jiayingzhu/database/qiime/silva-138-99-515-806-nb-classifier.qza \
        --i-reads ${workdir}/03_qiime/02_${denoising_method}_rep_seqs.qza \
        --o-classification ${workdir}/03_qiime/03_${denoising_method}_taxonomy.qza
    
    # full_length
        qiime feature-classifier classify-sklearn \
        --i-classifier /mnt/raid8/jiayingzhu/database/qiime/silva-138-99-nb-classifier.qza \
        --i-reads ${workdir}/03_qiime/02_${denoising_method}_rep_seqs.qza \
        --o-classification ${workdir}/03_qiime/03_${denoising_method}_taxonomy.qza

# 11. Calculate the relative abundance of each taxa at the corresponding taxonomic level for the sample based on the taxonomy hierarchy of ASVs (feature table obtained from dada2 denoising)
    denoising_method='dada2'  # 或者deblur
    
    # Genus水平
    qiime taxa collapse \
    --i-table ${workdir}/03_qiime/02_${denoising_method}_rep_seqs.qza \
    --i-taxonomy ${workdir}/03_qiime/03_${denoising_method}_taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table ${workdir}/03_qiime/04_${denoising_method}_genus_table.qza
    
    qiime tools export \
    --input-path ${workdir}/03_qiime/04_${denoising_method}_genus_table.qza \
    --output-path ${workdir}/03_qiime/04_${denoising_method}_genus_table
    
    cd ${workdir}/03_qiime/04_${denoising_method}_genus_table
    biom head -i feature-table.biom
    biom summarize-table -i feature-table.biom
    biom convert \
    -i feature-table.biom \
    -o 04_${denoising_method}_genus_table.tsv \
    --to-tsv
    sed -i '1d' 04_${denoising_method}_genus_table.tsv
    
    # species 水平
    qiime taxa collapse \
    --i-table ${workdir}/03_qiime/02_${denoising_method}_rep_seqs.qza \
    --i-taxonomy ${workdir}/03_qiime/03_${denoising_method}_taxonomy.qza \
    --p-level 7 \
    --o-collapsed-table ${workdir}/03_qiime/04_${denoising_method}_species_table.qza
    
    qiime tools export \
    --input-path ${workdir}/03_qiime/04_${denoising_method}_species_table.qza \
    --output-path ${workdir}/03_qiime/04_${denoising_method}_species_table
    
    cd ${workdir}/03_qiime/04_${denoising_method}_species_table
    biom head -i feature-table.biom
    biom summarize-table -i feature-table.biom
    biom convert \
    -i feature-table.biom \
    -o 04_${denoising_method}_species_table.tsv \
    --to-tsv
    sed -i '1d' 04_${denoising_method}_species_table.tsv
    
