import os
import pysam
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

if config == {}:
    configfile: "%s/config.yaml" % SNAKEMAKE_DIR

if not os.path.exists("log"):
    os.makedirs("log")

localrules: extract_query_sequence, merge_regions, tile_regions, filter_low_mapq_and_splits, get_length_and_mismatches, summary_stats, bgzip_regions, get_het_differences, get_vcf_mapping_overlap, get_corrected_mismatches, get_fixed_stats

rule summary_stats:
    input: "results/stats.tab", "results/fixed_regions.tab"
    output: "results/summary.tab"
    params: query_assembly_error=config["query_assembly_error_rate"], query_name=config["query_assembly_name"]
    run:
        from math import log10
        hets = pd.read_table(input[1])
        corrections = hets.groupby(["chr", "start", "end", "ed"])["new_ed"].agg("min").reset_index()
        het_corrections = corrections.ed.sum() - corrections.new_ed.sum()

        dat = pd.read_table(input[0])
        total_bases = dat.length.sum()
        differences = dat.edist.sum()
        error_rate = (differences - het_corrections) / total_bases - params.query_assembly_error
        phred_score = -10 * log10(error_rate)
        with open(output[0], "w") as outfile:
            print("Total bases:", total_bases, file=outfile)
            print("Total diff bases:", differences, file=outfile)
            print("Heterozygous differences:", het_corrections, file=outfile)
            print(params.query_name, "error rate:", params.query_assembly_error, file=outfile)
            print("Initial:", -10 * log10(differences/total_bases), file=outfile)
            print("Het corrected:", -10 * log10((differences-het_corrections)/total_bases), file=outfile)
            print("query_assembly+Het corrected:", -10*log10(error_rate), file=outfile)
            print("equal QV corrected:", -10*log10(error_rate/2), file=outfile)

rule get_length_and_mismatches:
    input: "mapping/hits.filt.bam"
    output: "results/stats.tab"
    run:
        samfile = pysam.AlignmentFile(input[0], "rb")
        with open(output[0], "w") as outfile:
            print("chr", "start", "end", "length", "edist", "cigar", "md", sep="\t", file=outfile)
            for record in samfile:
                print(record.reference_name,
                      record.reference_start,
                      record.reference_end,
                      record.infer_query_length(),
                      record.get_tag("NM"),
                      record.cigarstring,
                      record.get_tag("MD"),
                      sep="\t",
                      file=outfile)

rule get_corrected_mismatches:
    input: "results/het_differences.stats.tab"
    output: "results/corrected_hets.tab"
    run:
        dat = pd.read_table(input[0])
        out = dat.groupby(["chr", "start", "end", "edit_distance", "pos", "ref", "alt_aln", "cigar", "md"])["alt"].agg("count").reset_index()
        out.to_csv(output[0], sep="\t", index=False)

rule get_het_differences:
    input: "results/overlap_vcf_mapping.stats.tab"
    output: "results/het_differences.stats.tab"
    run:
        def get_tuples_from_cigar(cigar):
            tuples = []
            annot = [x for x in cigar if x not in "0123456789"]
            for i, x in enumerate(annot):
                index = cigar.index(x)
                tuples.append((int(cigar[:index]), x))
                if i < len(annot):
                    cigar = cigar[index+1:]
            return(tuples)

        def format_md(md):
            tuples = []
            current_num = 0
            is_del = False
            current_del = ""
            for i, char in enumerate(md):
                if char == "^":
                    is_del = True
                    if current_num != 0:
                        tuples.append((current_num, "-"))
                        current_num = 0
                elif char in "0123456789":
                    if is_del:
                        tuples.append((current_num, current_del))
                        is_del = False
                        current_num = 0
                    current_num = current_num * 10 + int(char)
                else:
                    if is_del:
                        current_num -= 1
                        current_del += char
                    elif current_num > 0:
                        tuples.append((current_num, "-"))
                        current_num = 0
                    tuples.append((1, char))
            if current_num > 0:
                tuples.append((current_num, "-"))
            if is_del:
                tuples.append(current_num, current_del)
            return tuples

        dat = pd.read_table(input[0])
        dat["cigar_tuples"] = dat.cigar.apply(lambda x: get_tuples_from_cigar(x))
        dat["md_tuples"] = dat.md.apply(lambda x: format_md(x))
        subset = dat[(dat.edit_distance > 0) & (dat.alt_aln == dat.alt)]
        subset.to_csv(output[0], sep="\t", index=False)

rule get_fixed_stats:
    input: vcf="vcf/reference_variants.vcf.gz", bam="mapping/hits.filt.bam"
    output: "results/fixed_regions.tab"
    params: indiv=config["vcf_ref_indiv"]
    run:
        def vcf_filter_pass(record):
            return record.samples[params.indiv]["GT"] == (0, 1) and record.info["TYPE"][0] in ["snp", "ins", "del"]
        samfile = pysam.AlignmentFile(input.bam, "rb")
        vcf_reader = pysam.VariantFile(input.vcf)
        outfile = open(output[0], "w")
        print("chr", "start", "end", "ed", "new_ed", "pos", "ref", "alt", sep="\t", file=outfile)
        for read in samfile:
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            ref_seq = read.get_reference_sequence()
            edit_distance = read.get_tag("NM")
            if edit_distance < 1:
                continue
            temp_seq = ref_seq
            temp_ed = edit_distance
            for record in vcf_reader.fetch(chrom, start, end):
                if not vcf_filter_pass(record):
                    continue
                rec_start = record.pos
                rec_end = record.pos + len(record.ref)
                new_seq = temp_seq[:rec_start - start - 1] + record.alts[0] + temp_seq[rec_end - start - 1:]
                new_alt = new_seq[rec_start - start - 1:rec_start - start - 1 + len(record.alts[0])]
                new_dist = sum([1 for (query, new_ref) in zip(read.seq, new_seq) if query != new_ref])
                if new_dist < temp_ed:
                    temp_seq = new_seq
                    temp_ed = new_dist
                    print(chrom, start, end, edit_distance, new_dist, record.pos, record.ref, record.alts[0], sep="\t", file=outfile)
        outfile.close()
        samfile.close()

rule get_vcf_mapping_overlap:
    input: vcf="vcf/reference_variants.vcf.gz", bam="mapping/hits.filt.bam"
    output: "results/overlap_vcf_mapping.stats.tab"
    params: indiv=config["vcf_ref_indiv"]
    run:
        def vcf_filter_pass(record):
            "Confirm indiv is heterozygous and variant is a snp or indel"
            return record.samples[params.indiv]["GT"] == (0, 1) and record.info["TYPE"][0] in ["snp", "ins", "del"]

        samfile = pysam.AlignmentFile(input.bam, "rb")
        vcf_reader = pysam.VariantFile(input.vcf)
        with open(output[0], "w") as outfile:
            print("chr", "start", "end", "length", "edit_distance", "cigar", "md", "hg38_region", "pos", "ref", "alt", "ref_aln", "alt_aln", "context", "query_seq", sep="\t", file=outfile)
            for read in samfile:
                chrom = read.reference_name
                region_start = read.reference_start
                region_end = read.reference_end
                region_string = "{}:{}-{}".format(chrom, region_start, region_end)
                for record in vcf_reader.fetch(chrom, region_start, region_end):
                    if not vcf_filter_pass(record):
                        continue
                    diff = len(record.alts[0]) - len(record.ref)
                    start = min(record.pos + diff, record.pos)
                    end = max(record.pos + diff, record.pos) + 1
                    ref_length = record.rlen
                    alt_length = len(record.alts[0])
                    query_start = start - read.reference_start - 1
                    seq_ref_aln = "".join([base for base in read.seq[query_start:query_start + ref_length] if base is not None])
                    seq_alt_aln = "".join([base for base in read.seq[query_start:query_start + alt_length] if base is not None])

                    context = "".join([base for base in read.seq[query_start-3:query_start+ref_length+3] if base is not None])
                    print(chrom,
                          region_start,
                          region_end,
                          read.infer_query_length(),
                          read.get_tag("NM"),
                          read.cigarstring,
                          read.get_tag("MD"),
                          read.query_name,
                          record.pos,
                          record.ref,
                          record.alts[0],
                          seq_ref_aln,
                          seq_alt_aln,
                          context,
                          read.seq,
                          sep="\t",
                          file=outfile)
        samfile.close()

rule get_reference_vcf:
    input: vcf=config["reference_vcf"], regions="regions/clint.hits.filt.bed.gz"
    output: "vcf/reference_variants.vcf.gz"
    params: indiv=config["vcf_ref_indiv"], sge_opts="-l mfree=8G -l h_rt=2:0:0:0"
    shell:
        """tabix -h {input.vcf} -R {input.regions} | vcffilter -f "QUAL > 20" |\
           vcfkeepsamples - {params.indiv} | vcffixup - | vcffilter -g "GT = 0/1" -f "AC = 1" |\
           vcfkeepinfo - AC AB TYPE > $TMPDIR/temp.vcf
           vcfstreamsort -a $TMPDIR/temp.vcf | bgzip -c > {output}
           tabix {output}"""

rule get_clint_bed:
    input: "mapping/hits.filt.bam"
    output: "regions/clint.hits.filt.bed", "regions/clint.hits.filt.bed.gz"
    params: sge_opts="-l mfree=1G -l h_rt=1:0:0"
    shell:
        "bedtools bamtobed -i {input} | sort -k 1,1 -k 2,2n -k 3,3n > {output[0]}; "
        "bedtools bamtobed -i {input} | sort -k 1,1 -k 2,2n -k 3,3n | bgzip -c > {output[1]}; "
        "tabix {output[1]}"

rule filter_low_mapq_and_splits:
    input: "mapping/hits.bam"
    output: "mapping/hits.filt.bam", "mapping/hits.filt.log"
    params: min_qual=40
    run:
        def is_split(read):
            clip_indices = [4, 5] # 4 is soft clip, 5 is hard clip
            return read.cigartuples[0][0] in clip_indices or read.cigartuples[-1][0] in clip_indices

        samfile = pysam.AlignmentFile(input[0], "rb")
        outfile = pysam.AlignmentFile(output[0], "wb", template=samfile)
        logfile = open(output[1], "w")
        for record in samfile:
            if record.mapping_quality >= params.min_qual and not is_split(record):
                outfile.write(record)
            else:
                print("Record excluded: {} Qual: {} Cigar: {}".format(record.query_name, record.mapping_quality, record.cigarstring), file=logfile)

        samfile.close()
        outfile.close()
        logfile.close()

rule map_query_to_reference:
    input: seq="sequence/query_regions.fasta", ref=config["reference_sequence"]
    output: "mapping/hits.bam"
    params: sge_opts="-l mfree=4G -l h_rt=1:0:0 -pe serial 4"
    shell:
        "bwa mem -t 4 -x intractg {input.ref} {input.seq} | samblaster | samtools sort -@ 4 -O bam -T $TMPDIR/hits -o {output}"

rule extract_query_sequence:
    input: regions="regions/regions.tiled.bed", seq=config["query_sequence"]
    output: "sequence/query_regions.fasta"
    params: sge_opts="-l mfree=4G -l h_rt=1:0:0"
    shell:
        "bedtools getfasta -fi {input.seq} -bed {input.regions} -fo {output}"

rule bgzip_regions:
    input: "regions/regions.tiled.bed"
    output: "regions/regions.tiled.bed.gz"
    shell:
        "bgzip -c {input} > {output}; "
        "tabix {output}"

rule tile_regions:
    input: regions="regions/regions.merged.bed"
    output: "regions/regions.tiled.bed"
    params: sge_opts="-l mfree=4G -l h_rt=1:0:0"
    shell:
        "bedtools makewindows -w 500 -b {input.regions} > {output}"

rule merge_regions:
    input: regions=config["query_regions"]
    output: "regions/regions.merged.bed"
    shell:
        "bedtools merge -i {input} > {output}"
