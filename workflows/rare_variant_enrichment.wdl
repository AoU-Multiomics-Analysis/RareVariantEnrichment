version 1.0

task PrepareVcfIndex {
    input {
        File rare_variant_vcf
        File? rare_variant_vcf_tbi
        Array[String] chromosomes
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    command <<<
        set -euo pipefail

        ln -s "~{rare_variant_vcf}" variants.vcf.gz
        if [[ "~{defined(rare_variant_vcf_tbi)}" == "true" ]]; then
            ln -s "~{select_first([rare_variant_vcf_tbi, rare_variant_vcf])}" variants.vcf.gz.tbi
        else
            tabix -p vcf variants.vcf.gz
        fi

        tabix -l variants.vcf.gz > indexed_contigs.txt
        for chromosome in ~{sep=" " chromosomes}; do
            if ! grep -F -x -q -- "$chromosome" indexed_contigs.txt; then
                echo "Requested chromosome is absent from VCF index: $chromosome" >&2
                exit 1
            fi
        done

        tabix -H variants.vcf.gz \
            | awk -F '\t' '$1 == "#CHROM" { for (column = 10; column <= NF; column++) print $column }' \
            > vcf_samples.txt
    >>>

    output {
        File vcf = "variants.vcf.gz"
        File vcf_tbi = "variants.vcf.gz.tbi"
        File vcf_samples = "vcf_samples.txt"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}

task PreparePhenotypes {
    input {
        File phenotype_bed
        File vcf_samples
        Array[String] chromosomes
        Array[Float] z_thresholds
        String outlier_tail
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    command <<<
        set -euo pipefail

        rare-variant-enrichment prepare-phenotypes \
            --phenotype-bed "~{phenotype_bed}" \
            --vcf-samples "~{vcf_samples}" \
            --chromosomes "~{sep="," chromosomes}" \
            --z-thresholds "~{sep="," z_thresholds}" \
            --tail "~{outlier_tail}" \
            --feature-output features.tsv \
            --sample-output shared_samples.txt \
            --qc-output phenotype_qc.json
    >>>

    output {
        File features_tsv = "features.tsv"
        File shared_samples = "shared_samples.txt"
        File phenotype_qc_json = "phenotype_qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}

task DetermineMaximumDistance {
    input {
        Array[Int] distance_thresholds_bp
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    command <<<
        set -euo pipefail

        maximum_distance_bp=-1
        for distance_bp in ~{sep=" " distance_thresholds_bp}; do
            if (( distance_bp > maximum_distance_bp )); then
                maximum_distance_bp="$distance_bp"
            fi
        done
        if (( maximum_distance_bp < 0 )); then
            echo "At least one non-negative distance threshold is required" >&2
            exit 1
        fi
        printf '%s\n' "$maximum_distance_bp" > maximum_distance_bp.txt
    >>>

    output {
        Int maximum_distance_bp = read_int("maximum_distance_bp.txt")
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}

task ClassifyChromosome {
    input {
        File rare_variant_vcf
        File rare_variant_vcf_tbi
        File features_tsv
        File shared_samples
        String chromosome
        Array[Int] exact_allele_counts
        Array[Int] cumulative_allele_count_maxima
        Int maximum_distance_bp
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    command <<<
        set -euo pipefail

        ln -s "~{rare_variant_vcf}" variants.vcf.gz
        ln -s "~{rare_variant_vcf_tbi}" variants.vcf.gz.tbi

        rare-variant-enrichment classify-chromosome \
            --vcf variants.vcf.gz \
            --features "~{features_tsv}" \
            --shared-samples "~{shared_samples}" \
            --chromosome "~{chromosome}" \
            --exact-ac "~{sep="," exact_allele_counts}" \
            --cumulative-ac-max "~{sep="," cumulative_allele_count_maxima}" \
            --max-distance "~{maximum_distance_bp}" \
            --carrier-output carrier_pairs.tsv \
            --regions-output query_regions.bed \
            --qc-output chromosome_qc.json
    >>>

    output {
        File carrier_pairs_tsv = "carrier_pairs.tsv"
        File query_regions_bed = "query_regions.bed"
        File chromosome_qc_json = "chromosome_qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}

task GatherCarrierPairs {
    input {
        Array[File] carrier_pairs
        Array[File] chromosome_qc
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    command <<<
        set -euo pipefail

        rare-variant-enrichment gather \
            ~{sep=" " prefix("--carrier-input ", carrier_pairs)} \
            ~{sep=" " prefix("--qc-input ", chromosome_qc)} \
            --carrier-output carrier_minimum_distances.tsv \
            --qc-output chromosome_qc.tsv
    >>>

    output {
        File carrier_minimum_distances_tsv = "carrier_minimum_distances.tsv"
        File chromosome_qc_tsv = "chromosome_qc.tsv"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}

task CalculateEnrichment {
    input {
        File phenotype_bed
        File shared_samples
        File carrier_minimum_distances_tsv
        Array[Int] exact_allele_counts
        Array[Int] cumulative_allele_count_maxima
        Array[Float] z_thresholds
        Array[Int] distance_thresholds_bp
        String outlier_tail
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    command <<<
        set -euo pipefail

        rare-variant-enrichment calculate \
            --phenotype-bed "~{phenotype_bed}" \
            --shared-samples "~{shared_samples}" \
            --carriers "~{carrier_minimum_distances_tsv}" \
            --exact-ac "~{sep="," exact_allele_counts}" \
            --cumulative-ac-max "~{sep="," cumulative_allele_count_maxima}" \
            --z-thresholds "~{sep="," z_thresholds}" \
            --distance-thresholds "~{sep="," distance_thresholds_bp}" \
            --tail "~{outlier_tail}" \
            --output-tsv enrichment.tsv \
            --output-json enrichment.json
    >>>

    output {
        File enrichment_tsv = "enrichment.tsv"
        File enrichment_json = "enrichment.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}

workflow RareVariantEnrichment {
    input {
        File phenotype_bed
        File rare_variant_vcf
        File? rare_variant_vcf_tbi
        Array[String] chromosomes
        Array[Float] z_thresholds = [2.0, 3.0, 4.0, 5.0]
        Array[Int] exact_allele_counts = [1, 2, 3, 4, 5]
        Array[Int] cumulative_allele_count_maxima = [1, 2, 3, 5, 10]
        Array[Int] distance_thresholds_bp = [1000, 10000, 100000, 1000000]
        String outlier_tail = "absolute"
        String docker_image = "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main"
        Int prepare_cpu = 2
        Int prepare_memory_gb = 8
        Int prepare_disk_gb = 50
        Int scatter_cpu = 2
        Int scatter_memory_gb = 8
        Int scatter_disk_gb = 20
        Int gather_cpu = 2
        Int gather_memory_gb = 16
        Int gather_disk_gb = 50
    }

    call PrepareVcfIndex {
        input:
            rare_variant_vcf = rare_variant_vcf,
            rare_variant_vcf_tbi = rare_variant_vcf_tbi,
            chromosomes = chromosomes,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = prepare_disk_gb
    }

    call PreparePhenotypes {
        input:
            phenotype_bed = phenotype_bed,
            vcf_samples = PrepareVcfIndex.vcf_samples,
            chromosomes = chromosomes,
            z_thresholds = z_thresholds,
            outlier_tail = outlier_tail,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = prepare_disk_gb
    }

    call DetermineMaximumDistance {
        input:
            distance_thresholds_bp = distance_thresholds_bp,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = prepare_disk_gb
    }

    scatter (chromosome in chromosomes) {
        call ClassifyChromosome {
            input:
                rare_variant_vcf = PrepareVcfIndex.vcf,
                rare_variant_vcf_tbi = PrepareVcfIndex.vcf_tbi,
                features_tsv = PreparePhenotypes.features_tsv,
                shared_samples = PreparePhenotypes.shared_samples,
                chromosome = chromosome,
                exact_allele_counts = exact_allele_counts,
                cumulative_allele_count_maxima = cumulative_allele_count_maxima,
                maximum_distance_bp = DetermineMaximumDistance.maximum_distance_bp,
                docker_image = docker_image,
                cpu = scatter_cpu,
                memory_gb = scatter_memory_gb,
                disk_gb = scatter_disk_gb
        }
    }

    call GatherCarrierPairs {
        input:
            carrier_pairs = ClassifyChromosome.carrier_pairs_tsv,
            chromosome_qc = ClassifyChromosome.chromosome_qc_json,
            docker_image = docker_image,
            cpu = gather_cpu,
            memory_gb = gather_memory_gb,
            disk_gb = gather_disk_gb
    }

    call CalculateEnrichment {
        input:
            phenotype_bed = phenotype_bed,
            shared_samples = PreparePhenotypes.shared_samples,
            carrier_minimum_distances_tsv = GatherCarrierPairs.carrier_minimum_distances_tsv,
            exact_allele_counts = exact_allele_counts,
            cumulative_allele_count_maxima = cumulative_allele_count_maxima,
            z_thresholds = z_thresholds,
            distance_thresholds_bp = distance_thresholds_bp,
            outlier_tail = outlier_tail,
            docker_image = docker_image,
            cpu = gather_cpu,
            memory_gb = gather_memory_gb,
            disk_gb = gather_disk_gb
    }

    output {
        File enrichment_tsv = CalculateEnrichment.enrichment_tsv
        File enrichment_json = CalculateEnrichment.enrichment_json
        File chromosome_qc_tsv = GatherCarrierPairs.chromosome_qc_tsv
        File carrier_minimum_distances_tsv = GatherCarrierPairs.carrier_minimum_distances_tsv
        File generated_or_validated_vcf_tbi = PrepareVcfIndex.vcf_tbi
        Array[File] chromosome_query_regions = ClassifyChromosome.query_regions_bed
    }
}
