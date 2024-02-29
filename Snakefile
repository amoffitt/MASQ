configfile: "config.yaml" # loads info as config dictionary

# one sample name
if isinstance(config["samples"],str):
    config["samples"] = [config["samples"]]

# one fastq directory
if "fastq_dirs" in config:
    if isinstance(config["fastq_dirs"],str):
        newdict={}
        for s in config["samples"]:
            newdict[s] = config["fastq_dirs"]
        config["fastq_dirs"] = newdict

# Allow for locus table or filtered primer design table
# Check header of SNV_table
with open(config['SNV_table']) as f:
    header = f.readline().split()
    if "SNP_ID" in header:
        # primer design table, needs transforming
        config["convert_SNV_table"]=True
    else:
        # already transformed
        config["convert_SNV_table"]=False

# get region number from snp input file 
if "regions" not in config:
    with open(config["SNV_table"],'r') as f:
        num_lines=0
        for line in f:
            num_lines += 1
        newlist = [str(x).zfill(3) for x in range(num_lines-1)]
        config["regions"] = newlist
else:
    # expand region number to region list
    if isinstance(config["regions"],int):
        regnum = config["regions"]
        newlist = [str(x).zfill(3) for x in range(regnum)]
        config["regions"] = newlist


# Allow for locus table or filtered primer design table
# Check header of SNV_table
with open(config['SNV_table']) as f:
    header = f.readline().split()
    if "SNP_ID" in header:
        # primer design table, needs transforming
        config["convert_SNV_table"]=True
    else:
        # already transformed
        config["convert_SNV_table"]=False


# When barcode splitting is needed
if "barcodes" in config:
    # Load barcode table
    barcode_dict=dict()
    with open(config["barcode_file"],'r') as f:
        for line in f:
            x=line.strip().split(":")
            barcode_dict[int(x[0])]=x[1]
    # get sample barcode list string
    sample_bc_string="\""
    for S,bc in config["barcodes"].items():
        sample_bc_string=sample_bc_string+S+":"+barcode_dict[bc]+" "
    sample_bc_string=sample_bc_string+"\""
    # Only works for one barcode group per sample
    config["fastq_dirs"]=dict()
    for S in config["samples"]:
        config["fastq_dirs"][S]="bc_split_fastqs/"+S
else:
    config["presplit_fastqs"]=["",""]
    sample_bc_string=""


# For backwards compatibility, set defaults
if "dna_input_ng" not in config:
    config["dna_input_ng"]=10000000
if "groupname" not in config:
    config["groupname"]=config["samples"][0]

######################################################################

rule all:
    input:
        expand("{sample}/plots/{sample}.WGS_varfreqs.{region}.png", sample=config["samples"], region=config["regions"]),
        expand("{sample}/plots/{sample}.hamming_distance_SP1_SP2.png", sample=config["samples"]),
        expand("{sample}/pickles/{sample}.vt_counter.{region}.pickle", zip, sample=config["samples"], region=config["regions"]),
        expand("{sample}/pickles/{sample}.vt_seq_counter.{region}.pickle", zip, sample=config["samples"], region=config["regions"]),
        expand("{sample}/pickles/{sample}.flip_counter.{region}.pickle", zip, sample=config["samples"], region=config["regions"]),
        expand("{sample}/fastqc/r1_fastqc.html", sample=config["samples"]),
        expand("{sample}/extended_var_table.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.base_count_allbases.perbase.combined.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.alignment_counter.region_{region}.txt", zip, sample=config["samples"], region=config["regions"]),
        expand("{sample}/reports/{sample}.base_count_allbases.region_{region}.txt", zip, sample=config["samples"], region=config["regions"]),
        expand("{sample}/reports/{sample}.withintagerrors.region_{region}.txt", zip, sample=config["samples"], region=config["regions"]),
        expand("{sample}/plots/{sample}.atleastX_reads_per_tag.region_{region}.png", zip, sample=config["samples"], region=config["regions"]),
        expand("{sample}/reports/{sample}.base_count_allbases.perbase.combined.txt", sample=config["samples"]),
        expand("{sample}/plots/{sample}.atleastX_reads_per_tag.allregions.png", sample=config["samples"]),
        expand("{sample}/reports/{sample}.alignment_counter.combined.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.base_count_allbases.combined.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.rollup_results.combined.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.base_count_variantbasesonly.txt", sample=config["samples"]),
        expand("{sample}/plots/{sample}.rollup_results.allregions.png", sample=config["samples"]),
        expand("{sample}/reports/{sample}.within_tag_errors.fractions.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.number_reads_per_tag.combined.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.number_reads_per_tag.allregions.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.final_report.txt", sample=config["samples"]),
        expand("{sample}/reports/{sample}.base_count_allbases.perbase.combined.qcfiltered.txt", sample=config["samples"]),
        "combined/"+config["groupname"]+".masq_QC_plots.png"

################################################################################
rule convert_SNV_table:
    input:
        oldtable=config["SNV_table"]
    output:
        newtable="locus_seq_table.txt"
    params:
        ref_genome = config["ref_genome_fa"]
    shell:
        "masq_primer_table_to_sd_table {input.oldtable} {output.newtable} {params.ref_genome}"


###############################################################################
rule barcode_split:
    input:
        presplit_fastq1=config["presplit_fastqs"][0],
        presplit_fastq2=config["presplit_fastqs"][1]
    output:
        output_fastqs=expand("bc_split_fastqs/{sample}/{read}.fastq.gz",read=["r1","r2"],sample=config["samples"])
    params:
        output_dir="bc_split_fastqs/",
        threads=16,
        sample_bc_list=sample_bc_string
    shell:
        """
        scripts/trim_bcs.sh {params.threads} {params.output_dir} {input.presplit_fastq1} {input.presplit_fastq2} {params.sample_bc_list}
        """

################################################################################
rule fastqc:
    input:
        fastq1=lambda wildcards: config["fastq_dirs"][wildcards.sample]+"/r1.fastq.gz",
        fastq2=lambda wildcards: config["fastq_dirs"][wildcards.sample]+"/r2.fastq.gz"
    output:
        "{sample}/fastqc/r1_fastqc.html",
        "{sample}/fastqc/r2_fastqc.html",
        temp("{sample}/fastqc/r1_fastqc.zip"),
        temp("{sample}/fastqc/r2_fastqc.zip")
    log:
        "{sample}/logs/fastqc.log"
    threads: 2
    shell:
        "(fastqc --nogroup -t {threads} {input} -o {wildcards.sample}/fastqc) >& {log}"

################################################################################
rule check_loci_extend_and_plot:
    input:
        bam=config["WGS_bam"],
        SNV_table="locus_seq_table.txt" if config["convert_SNV_table"] else config["SNV_table"] # allows 2 different format inputs
    output:
        new_SNV_table="{sample}/extended_var_table.txt",
        plots = expand("{{sample}}/plots/{{sample}}.WGS_varfreqs.{region}.png", region=config["regions"])
    params:
        wgs_name = config["WGS_name"],
        wgs_ref = config["WGS_ref"],
    log:
        "{sample}/logs/check_loci_plot_and_extend.log"
    script:
        "scripts/check_loci_plot_and_extend.py"

################################################################################

rule sort_data_by_tag_and_locus:
    input:
        fastq1=lambda wildcards: config["fastq_dirs"][wildcards.sample]+"/r1.fastq.gz",
        fastq2=lambda wildcards: config["fastq_dirs"][wildcards.sample]+"/r2.fastq.gz",
        SNV_table=lambda wildcards: wildcards.sample+"/extended_var_table.txt"
    output:
        vt_counters=expand("{{sample}}/pickles/{{sample}}.vt_counter.{region}.pickle", region=config["regions"]),
        vt_seq_counters=expand("{{sample}}/pickles/{{sample}}.vt_seq_counter.{region}.pickle", region=config["regions"]),
        flip_counters=expand("{{sample}}/pickles/{{sample}}.flip_counter.{region}.pickle", region=config["regions"]),
        ss_hamming_plot="{sample}/plots/{sample}.hamming_distance_SP1_SP2.png",
        counter_report="{sample}/reports/{sample}.primer_counters.txt",
        up2_unmatched_report="{sample}/reports/{sample}.unmatched_UP2_seqs.txt",
        ss1ss2_unmatched_report="{sample}/reports/{sample}.unmatched_SS1_SS2_seqs.txt",
        goodtag_report="{sample}/reports/{sample}.good_tags.txt",
        badtag_report="{sample}/reports/{sample}.bad_tags.txt"
    params:
        min_len= config["min_len"],
        trim_len= config["trim_len"],
        SS_sum_hamming= config["SS_sum_hamming"],
        UP2_hamming= config["UP2_hamming"],
        tag= config["tag"],
        UP2= config["UP2"],
        protocol= config["protocol"],
        quick_run= config["quick_run"],
        quick_run_reads = config["quick_run_reads"],
        mask_lowqual_bases = config["mask_lowqual_bases"],
        qual_cutoff = config["qual_cutoff"],
        max_N_ratio = config["max_N_ratio"],
        use_edit_distance = config["use_edit_distance"]
    log:
        "{sample}/logs/sort_data_by_tag_and_locus.log"
    script:
        "scripts/sort_data_by_tag_and_locus.py"

################################################################################

rule all_base_report:
    input:
        vt_counter="{sample}/pickles/{sample}.rolledup.vt_counter.{region}.pickle" if config["dorollup"] else "{sample}/pickles/{sample}.vt_counter.{region}.pickle",
        vt_seq_counter="{sample}/pickles/{sample}.rolledup.vt_seq_counter.{region}.pickle" if config["dorollup"] else "{sample}/pickles/{sample}.vt_seq_counter.{region}.pickle",
        flip_counters="{sample}/pickles/{sample}.flip_counter.{region}.pickle"
    output:
        base_count_report="{sample}/reports/{sample}.base_count_allbases.region_{region}.txt",
        base_count_report_per_base="{sample}/reports/{sample}.base_count_allbases.perbase.region_{region}.txt",
        alignment_report="{sample}/reports/{sample}.alignment_counter.region_{region}.txt",
        withintagerrors="{sample}/pickles/{sample}.withintagerrors.region_{region}.pickle",
        withintagerrors_table="{sample}/reports/{sample}.withintagerrors.region_{region}.txt",
        unaligned_reads="{sample}/reports/{sample}.unaligned_reads.region_{region}.txt"
    params:
        SNV_table=lambda wildcards: wildcards.sample+"/extended_var_table.txt",
        region = lambda wildcards: wildcards.region,
        target_hamming = config["target_hamming"],
        base_error_rate = config["base_error_rate"],
        coverage_list = config["coverage_list"],
        mask_lowqual_bases = config["mask_lowqual_bases"],
        qual_cutoff = config["qual_cutoff"],
        max_N_ratio = config["max_N_ratio"],
        ref_genome = config["ref_genome"],
        tag= config["tag"],
        UP2= config["UP2"],
        trim_len= config["trim_len"],
        filter_ns = config['filter_ns']
    log:
        "{sample}/logs/all_base_report.{region}.log"
    script:
        "scripts/all_base_report.py" 

################################################################################
rule combine_within_tag_errors:
    input:
        within_tag_error_pickles=expand("{{sample}}/pickles/{{sample}}.withintagerrors.region_{region}.pickle", region=config["regions"])
    output:
        within_tag_table1="{sample}/reports/{sample}.within_tag_errors.counts.txt",
        within_tag_table2="{sample}/reports/{sample}.within_tag_errors.fractions.txt"
    log:
        "{sample}/logs/combine_within_tag_errors.log"
    script:
        "scripts/combine_withintagerr.py"

################################################################################

rule tag_count_graphs_per_region:
    input:
        vt_counter="{sample}/pickles/{sample}.rolledup.vt_counter.{region}.pickle" if config["dorollup"] else "{sample}/pickles/{sample}.vt_counter.{region}.pickle"
    output:
        plot1="{sample}/plots/{sample}.number_reads_per_tag.region_{region}.png",
        plot2="{sample}/plots/{sample}.atleastX_reads_per_tag.region_{region}.png",
        tagcounts="{sample}/reports/{sample}.number_reads_per_tag.region_{region}.txt"
    params:
        region = lambda wildcards: wildcards.region,
        sample = lambda wildcards: wildcards.sample
    log:
        "{sample}/logs/tag_counts_per_region.{region}.log"
    script:
        "scripts/tag_count_graphs.py"

################################################################################

rule tag_count_graphs_combined:
    input:
        vt_counters=expand("{{sample}}/pickles/{{sample}}.rolledup.vt_counter.{region}.pickle", region=config["regions"]) if config["dorollup"] else expand("{{sample}}/pickles/{{sample}}.vt_counter.{region}.pickle", region=config["regions"])
    output:
        plot1="{sample}/plots/{sample}.number_reads_per_tag.allregions.png",
        plot2="{sample}/plots/{sample}.atleastX_reads_per_tag.allregions.png",
        tagcounts="{sample}/reports/{sample}.number_reads_per_tag.allregions.txt"
    params:
        sample = lambda wildcards: wildcards.sample
    log:
        "{sample}/logs/tag_counts_allregions.total.log"
    script:
        "scripts/tag_count_graphs_allregions.py"

################################################################################

rule rollup_tags:
    input:
        vt_counter="{sample}/pickles/{sample}.vt_counter.{region}.pickle",
        vt_seq_counter="{sample}/pickles/{sample}.vt_seq_counter.{region}.pickle",
        flip_counter="{sample}/pickles/{sample}.flip_counter.{region}.pickle"
    output:
        vt_counter="{sample}/pickles/{sample}.rolledup.vt_counter.{region}.pickle",
        vt_seq_counter="{sample}/pickles/{sample}.rolledup.vt_seq_counter.{region}.pickle",
        flip_counter="{sample}/pickles/{sample}.rolledup.flip_counter.{region}.pickle",
        tagcounts="{sample}/reports/{sample}.rollup_results.region_{region}.txt"
    params:
        region = lambda wildcards: wildcards.region,
        sample = lambda wildcards: wildcards.sample,
        dna_input_ng = config["dna_input_ng"]
    log:
        "{sample}/logs/rollup_tags.region_{region}.log"
    script:
        "scripts/collapse_tags.py"

################################################################################

rule combine_reports_alignment:
    input:
        region_reports=expand("{{sample}}/reports/{{sample}}.alignment_counter.region_{region}.txt", region=config["regions"])
    output:
        combined_report="{sample}/reports/{sample}.alignment_counter.combined.txt"
    log:
        "{sample}/logs/combine_reports.alignment.log"
    script:
        "scripts/combine_reports.py"


################################################################################

rule combine_reports_basecount:
    input:
        region_reports=expand("{{sample}}/reports/{{sample}}.base_count_allbases.region_{region}.txt", region=config["regions"])
    output:
        combined_report="{sample}/reports/{sample}.base_count_allbases.combined.txt"
    log:
        "{sample}/logs/combine_reports.basecount.log"
    script:
        "scripts/combine_reports.py"

################################################################################

rule combine_reports_per_base:
    input:
        region_reports=expand("{{sample}}/reports/{{sample}}.base_count_allbases.perbase.region_{region}.txt", region=config["regions"])
    output:
        combined_report="{sample}/reports/{sample}.base_count_allbases.perbase.combined.txt"
    log:
        "{sample}/logs/combine_reports.perbase.log"
    script:
        "scripts/combine_reports.py"

################################################################################

rule combine_reports_rollup:
    input:
        region_reports=expand("{{sample}}/reports/{{sample}}.rollup_results.region_{region}.txt", region=config["regions"])
    output:
        combined_report="{sample}/reports/{sample}.rollup_results.combined.txt"
    log:
        "{sample}/logs/combine_reports.rollup.log"
    script:
        "scripts/combine_reports.py"

################################################################################

rule combine_reports_readspertag:
    input:
        region_reports=expand("{{sample}}/reports/{{sample}}.number_reads_per_tag.region_{region}.txt", region=config["regions"])+
            ["{sample}/reports/{sample}.number_reads_per_tag.allregions.txt"]
    output:
        combined_report="{sample}/reports/{sample}.number_reads_per_tag.combined.txt"
    log:
        "{sample}/logs/combine_reports.readspertag.log"
    script:
        "scripts/combine_reports.py"

################################################################################

rule extract_variant_bases_from_report:
    input:
        combined_report="{sample}/reports/{sample}.base_count_allbases.combined.txt"
    output:
        variant_report="{sample}/reports/{sample}.base_count_variantbasesonly.txt"
    log:
        "{sample}/logs/extract_variant_info.log"
    script:
        "scripts/extract_variant_info.py"

################################################################################

rule plot_rollup_results:
    input:
        combined_report="{sample}/reports/{sample}.rollup_results.combined.txt"
    output:
        plot1="{sample}/plots/{sample}.rollup_results.allregions.png"
    params:
        sample = lambda wildcards: wildcards.sample
    log:
        "{sample}/logs/plot_rollup_results.log"
    script:
        "scripts/plot_rollup_results.py"

################################################################################

rule final_report:
    input:
        input_snv_table = "{sample}/extended_var_table.txt",
        report_primers = "{sample}/reports/{sample}.primer_counters.txt",
        report_rollup = "{sample}/reports/{sample}.rollup_results.combined.txt",
        report_alignment = "{sample}/reports/{sample}.alignment_counter.combined.txt",
        report_variants = "{sample}/reports/{sample}.base_count_variantbasesonly.txt",
    output:
        combined_report="{sample}/reports/{sample}.final_report.txt"
    log:
        "{sample}/logs/final_report.log"
    script:
        "scripts/final_report.py"

################################################################################
rule bad_loci:
    output: "combined/bad_loci_list.txt"
    run:
        with open(output[0],'w') as f:
            if "bad_loci" in config:
                for loc in config["bad_loci"]:
                    f.write(str(loc)+"\n")

rule qc_plots:
    input:
        reports = expand("{sample}/reports/{sample}.final_report.txt",sample=config["samples"]),
        bad_loci = "combined/bad_loci_list.txt"
    output:
        plot = "combined/"+config["groupname"]+".masq_QC_plots.png",
        qcfail = "combined/"+config["groupname"]+".qc_fail_loci.txt"
    run:
        import os
        if config["groupname"] == "test_examples_standardPCR":
            for f in [output.plot, output.qcfail]:
                os.makedirs(os.path.dirname(f), exist_ok=True)
                with open(f, 'w') as empty_file:
                    pass
        else:
            shell("R_LIBS=""; Rscript scripts/masq_QC_plots.R {output.plot} {output.qcfail} {input.bad_loci} {input.reports}")
################################################################################

rule filter_base_report:
    input:
        qcfail = "combined/"+config["groupname"]+".qc_fail_loci.txt",
        base_report = "{sample}/reports/{sample}.base_count_allbases.perbase.combined.txt"
    output:
        filtered_base_report = "{sample}/reports/{sample}.base_count_allbases.perbase.combined.qcfiltered.txt"
    log:
        "{sample}/logs/filter_base_report.log"
    script:
        "scripts/qcfilter_report.py"


