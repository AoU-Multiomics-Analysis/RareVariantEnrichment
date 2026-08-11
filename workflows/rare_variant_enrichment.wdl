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
        Int max_retries
    }

    File chromosomes_file = write_lines(chromosomes)

    command <<<
        set -euo pipefail

        ln -s "~{rare_variant_vcf}" variants.vcf.gz
        if [[ "~{defined(rare_variant_vcf_tbi)}" == "true" ]]; then
            ln -s "~{select_first([rare_variant_vcf_tbi, rare_variant_vcf])}" variants.vcf.gz.tbi
            printf 'supplied\n' > index_provenance.txt
        else
            tabix -p vcf variants.vcf.gz
            printf 'generated\n' > index_provenance.txt
        fi

        tabix -l variants.vcf.gz > indexed_contigs.txt
        while IFS= read -r chromosome; do
            if ! grep -F -x -q -- "$chromosome" indexed_contigs.txt; then
                echo "Requested chromosome is absent from VCF index: $chromosome" >&2
                exit 1
            fi
        done < "~{chromosomes_file}"

        tabix -H variants.vcf.gz \
            | awk -F '\t' '$1 == "#CHROM" { for (column = 10; column <= NF; column++) print $column }' \
            > vcf_samples.txt
    >>>

    output {
        File vcf = rare_variant_vcf
        File vcf_tbi = if defined(rare_variant_vcf_tbi) then select_first([rare_variant_vcf_tbi]) else "variants.vcf.gz.tbi"
        File vcf_samples = "vcf_samples.txt"
        String index_provenance = read_string("index_provenance.txt")
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task PrepareVatIndex {
    input {
        File variant_annotation_table
        File? variant_annotation_table_tbi
        Array[String] chromosomes
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File chromosomes_file = write_lines(chromosomes)

    command <<<
        set -euo pipefail

        join_by_comma() {
            local IFS=,
            printf '%s' "$*"
        }

        chromosomes=()
        while IFS= read -r chromosome; do
            chromosomes+=("$chromosome")
        done < "~{chromosomes_file}"
        chromosomes_csv="$(join_by_comma "${chromosomes[@]}")"

        ln -s "~{variant_annotation_table}" annotations.tsv.bgz
        if [[ "~{defined(variant_annotation_table_tbi)}" == "true" ]]; then
            ln -s "~{select_first([variant_annotation_table_tbi, variant_annotation_table])}" annotations.tsv.bgz.tbi
            printf 'supplied\n' > index_provenance.txt
        else
            printf 'generated\n' > index_provenance.txt
        fi

        rare-variant-enrichment prepare-vat \
            --vat annotations.tsv.bgz \
            --chromosomes "$chromosomes_csv" \
            --schema-output vat_schema.json \
            --loftee-enabled-output loftee_enabled.txt
    >>>

    output {
        File vat = variant_annotation_table
        File vat_tbi = if defined(variant_annotation_table_tbi) then select_first([variant_annotation_table_tbi]) else "annotations.tsv.bgz.tbi"
        File vat_schema_json = "vat_schema.json"
        Boolean loftee_enabled = read_boolean("loftee_enabled.txt")
        String index_provenance = read_string("index_provenance.txt")
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
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
        Int max_retries
    }

    File chromosomes_file = write_lines(chromosomes)
    File z_thresholds_file = write_lines(z_thresholds)
    File outlier_tail_file = write_lines([outlier_tail])

    command <<<
        set -euo pipefail

        join_by_comma() {
            local IFS=,
            printf '%s' "$*"
        }

        chromosomes=()
        while IFS= read -r value; do
            chromosomes+=("$value")
        done < "~{chromosomes_file}"
        z_thresholds=()
        while IFS= read -r value; do
            z_thresholds+=("$value")
        done < "~{z_thresholds_file}"
        IFS= read -r outlier_tail < "~{outlier_tail_file}"

        chromosomes_csv="$(join_by_comma "${chromosomes[@]}")"
        z_thresholds_csv="$(join_by_comma "${z_thresholds[@]}")"

        rare-variant-enrichment prepare-phenotypes \
            --phenotype-bed "~{phenotype_bed}" \
            --vcf-samples "~{vcf_samples}" \
            --chromosomes "$chromosomes_csv" \
            --z-thresholds "$z_thresholds_csv" \
            --tail "$outlier_tail" \
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
        maxRetries: max_retries
    }
}

task DetermineMaximumDistance {
    input {
        Array[Int] distance_thresholds_bp
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File distance_thresholds_file = write_lines(distance_thresholds_bp)

    command <<<
        set -euo pipefail

        distance_thresholds_bp=()
        while IFS= read -r distance_bp; do
            distance_thresholds_bp+=("$distance_bp")
        done < "~{distance_thresholds_file}"

        maximum_distance_bp=-1
        for distance_bp in "${distance_thresholds_bp[@]}"; do
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
        maxRetries: max_retries
    }
}

task ClassifyChromosome {
    input {
        File rare_variant_vcf
        File rare_variant_vcf_tbi
        File variant_annotation_table
        File variant_annotation_table_tbi
        File vat_schema_json
        File features_tsv
        File shared_samples
        String chromosome
        Array[Int] exact_allele_counts
        Array[Int] cumulative_allele_count_maxima
        Array[String] consequence_classes
        Float maximum_gvs_maf
        Int maximum_distance_bp
        Int annotation_chunk_size_bp
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File chromosome_file = write_lines([chromosome])
    File exact_allele_counts_file = write_lines(exact_allele_counts)
    File cumulative_allele_count_maxima_file = write_lines(cumulative_allele_count_maxima)
    File consequence_classes_file = write_lines(consequence_classes)

    command <<<
        set -euo pipefail

        join_by_comma() {
            local IFS=,
            printf '%s' "$*"
        }

        IFS= read -r chromosome < "~{chromosome_file}"
        exact_allele_counts=()
        while IFS= read -r value; do
            exact_allele_counts+=("$value")
        done < "~{exact_allele_counts_file}"
        cumulative_allele_count_maxima=()
        while IFS= read -r value; do
            cumulative_allele_count_maxima+=("$value")
        done < "~{cumulative_allele_count_maxima_file}"
        consequence_classes=()
        while IFS= read -r value; do
            consequence_classes+=("$value")
        done < "~{consequence_classes_file}"

        exact_allele_counts_csv="$(join_by_comma "${exact_allele_counts[@]}")"
        cumulative_allele_count_maxima_csv="$(join_by_comma "${cumulative_allele_count_maxima[@]}")"
        consequence_classes_csv="$(join_by_comma "${consequence_classes[@]}")"

        ln -s "~{rare_variant_vcf}" variants.vcf.gz
        ln -s "~{rare_variant_vcf_tbi}" variants.vcf.gz.tbi
        ln -s "~{variant_annotation_table}" annotations.tsv.bgz
        ln -s "~{variant_annotation_table_tbi}" annotations.tsv.bgz.tbi

        rare-variant-enrichment classify-chromosome \
            --vcf variants.vcf.gz \
            --vat annotations.tsv.bgz \
            --vat-schema "~{vat_schema_json}" \
            --features "~{features_tsv}" \
            --shared-samples "~{shared_samples}" \
            --chromosome "$chromosome" \
            --exact-ac "$exact_allele_counts_csv" \
            --cumulative-ac-max "$cumulative_allele_count_maxima_csv" \
            --consequence-classes "$consequence_classes_csv" \
            --maximum-gvs-maf "~{maximum_gvs_maf}" \
            --max-distance "~{maximum_distance_bp}" \
            --annotation-chunk-size-bp "~{annotation_chunk_size_bp}" \
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
        maxRetries: max_retries
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
        Int max_retries
    }

    File carrier_pairs_file = write_lines(carrier_pairs)
    File chromosome_qc_file = write_lines(chromosome_qc)

    command <<<
        set -euo pipefail

        gather_arguments=()
        while IFS= read -r carrier_path; do
            gather_arguments+=(--carrier-input "$carrier_path")
        done < "~{carrier_pairs_file}"
        while IFS= read -r qc_path; do
            gather_arguments+=(--qc-input "$qc_path")
        done < "~{chromosome_qc_file}"

        rare-variant-enrichment gather \
            "${gather_arguments[@]}" \
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
        maxRetries: max_retries
    }
}

task PublishCarrierAudit {
    input {
        File source_carrier_minimum_distances_tsv
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    command <<<
        set -euo pipefail

        cp "~{source_carrier_minimum_distances_tsv}" carrier_minimum_distances.tsv
    >>>

    output {
        File carrier_minimum_distances_tsv = "carrier_minimum_distances.tsv"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task CalculateEnrichment {
    input {
        File phenotype_bed
        File shared_samples
        File carrier_minimum_distances_tsv
        File features_tsv
        File phenotype_qc_json
        File chromosome_qc_tsv
        Array[String] selected_chromosomes
        Array[Int] exact_allele_counts
        Array[Int] cumulative_allele_count_maxima
        Array[Float] z_thresholds
        Array[Int] distance_thresholds_bp
        String outlier_tail
        Array[String] consequence_classes
        Boolean loftee_enabled
        String vcf_index_provenance
        String vat_index_provenance
        Float maximum_gvs_maf
        Int annotation_chunk_size_bp
        String workflow_version
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File exact_allele_counts_file = write_lines(exact_allele_counts)
    File cumulative_allele_count_maxima_file = write_lines(cumulative_allele_count_maxima)
    File z_thresholds_file = write_lines(z_thresholds)
    File distance_thresholds_bp_file = write_lines(distance_thresholds_bp)
    File outlier_tail_file = write_lines([outlier_tail])
    File selected_chromosomes_file = write_lines(selected_chromosomes)
    File container_image_file = write_lines([docker_image])
    File consequence_classes_file = write_lines(consequence_classes)
    File loftee_enabled_file = write_lines([loftee_enabled])
    File vcf_index_provenance_file = write_lines([vcf_index_provenance])
    File vat_index_provenance_file = write_lines([vat_index_provenance])
    File maximum_gvs_maf_file = write_lines([maximum_gvs_maf])
    File annotation_chunk_size_bp_file = write_lines([annotation_chunk_size_bp])
    File workflow_version_file = write_lines([workflow_version])

    command <<<
        set -euo pipefail

        join_by_comma() {
            local IFS=,
            printf '%s' "$*"
        }

        exact_allele_counts=()
        while IFS= read -r value; do
            exact_allele_counts+=("$value")
        done < "~{exact_allele_counts_file}"
        cumulative_allele_count_maxima=()
        while IFS= read -r value; do
            cumulative_allele_count_maxima+=("$value")
        done < "~{cumulative_allele_count_maxima_file}"
        z_thresholds=()
        while IFS= read -r value; do
            z_thresholds+=("$value")
        done < "~{z_thresholds_file}"
        distance_thresholds_bp=()
        while IFS= read -r value; do
            distance_thresholds_bp+=("$value")
        done < "~{distance_thresholds_bp_file}"
        IFS= read -r outlier_tail < "~{outlier_tail_file}"
        selected_chromosomes=()
        while IFS= read -r value; do
            selected_chromosomes+=("$value")
        done < "~{selected_chromosomes_file}"
        consequence_classes=()
        while IFS= read -r value; do
            consequence_classes+=("$value")
        done < "~{consequence_classes_file}"
        IFS= read -r container_image < "~{container_image_file}"
        IFS= read -r loftee_enabled < "~{loftee_enabled_file}"
        IFS= read -r vcf_index_provenance < "~{vcf_index_provenance_file}"
        IFS= read -r vat_index_provenance < "~{vat_index_provenance_file}"
        IFS= read -r maximum_gvs_maf < "~{maximum_gvs_maf_file}"
        IFS= read -r annotation_chunk_size_bp < "~{annotation_chunk_size_bp_file}"
        IFS= read -r workflow_version < "~{workflow_version_file}"

        exact_allele_counts_csv="$(join_by_comma "${exact_allele_counts[@]}")"
        cumulative_allele_count_maxima_csv="$(join_by_comma "${cumulative_allele_count_maxima[@]}")"
        z_thresholds_csv="$(join_by_comma "${z_thresholds[@]}")"
        distance_thresholds_bp_csv="$(join_by_comma "${distance_thresholds_bp[@]}")"
        selected_chromosomes_csv="$(join_by_comma "${selected_chromosomes[@]}")"
        consequence_classes_csv="$(join_by_comma "${consequence_classes[@]}")"

        rare-variant-enrichment calculate \
            --phenotype-bed "~{phenotype_bed}" \
            --shared-samples "~{shared_samples}" \
            --carriers "~{carrier_minimum_distances_tsv}" \
            --features "~{features_tsv}" \
            --exact-ac "$exact_allele_counts_csv" \
            --cumulative-ac-max "$cumulative_allele_count_maxima_csv" \
            --z-thresholds "$z_thresholds_csv" \
            --distance-thresholds "$distance_thresholds_bp_csv" \
            --tail "$outlier_tail" \
            --output-tsv enrichment.tsv \
            --output-json enrichment.json \
            --phenotype-qc "~{phenotype_qc_json}" \
            --chromosome-qc "~{chromosome_qc_tsv}" \
            --selected-chromosomes "$selected_chromosomes_csv" \
            --container-image "$container_image" \
            --workflow-version "$workflow_version" \
            --max-retries "~{max_retries}" \
            --index-provenance "$vcf_index_provenance" \
            --consequence-classes "$consequence_classes_csv" \
            --loftee-enabled "$loftee_enabled" \
            --vat-index-provenance "$vat_index_provenance" \
            --maximum-gvs-maf "$maximum_gvs_maf" \
            --annotation-chunk-size-bp "$annotation_chunk_size_bp"
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
        maxRetries: max_retries
    }
}

workflow RareVariantEnrichment {
    input {
        File phenotype_bed
        File rare_variant_vcf
        File? rare_variant_vcf_tbi
        File variant_annotation_table
        File? variant_annotation_table_tbi
        Array[String] chromosomes = [
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
            "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
            "chr19", "chr20", "chr21", "chr22"
        ]
        Array[Float] z_thresholds = [2.0, 3.0, 4.0, 5.0]
        Array[Int] exact_allele_counts = [1, 2, 3, 4, 5]
        Array[Int] cumulative_allele_count_maxima = [1, 2, 3, 5, 10]
        Array[String] consequence_classes = [
            "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
            "frameshift_variant", "stop_lost", "start_lost",
            "inframe_insertion", "inframe_deletion", "missense_variant",
            "protein_altering_variant", "splice_region_variant",
            "synonymous_variant", "coding_sequence_variant"
        ]
        Float maximum_gvs_maf = 0.01
        Int annotation_chunk_size_bp = 10000000
        Array[Int] distance_thresholds_bp = [10000, 100000]
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
        Int max_retries = 1
        Boolean publish_carrier_audit = false
    }

    String workflow_version = "0.3.0"

    call PrepareVcfIndex {
        input:
            rare_variant_vcf = rare_variant_vcf,
            rare_variant_vcf_tbi = rare_variant_vcf_tbi,
            chromosomes = chromosomes,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = prepare_disk_gb,
            max_retries = max_retries
    }

    call PrepareVatIndex {
        input:
            variant_annotation_table = variant_annotation_table,
            variant_annotation_table_tbi = variant_annotation_table_tbi,
            chromosomes = chromosomes,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = prepare_disk_gb,
            max_retries = max_retries
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
            disk_gb = prepare_disk_gb,
            max_retries = max_retries
    }

    call DetermineMaximumDistance {
        input:
            distance_thresholds_bp = distance_thresholds_bp,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = prepare_disk_gb,
            max_retries = max_retries
    }

    scatter (chromosome in chromosomes) {
        call ClassifyChromosome {
            input:
                rare_variant_vcf = PrepareVcfIndex.vcf,
                rare_variant_vcf_tbi = PrepareVcfIndex.vcf_tbi,
                variant_annotation_table = PrepareVatIndex.vat,
                variant_annotation_table_tbi = PrepareVatIndex.vat_tbi,
                vat_schema_json = PrepareVatIndex.vat_schema_json,
                features_tsv = PreparePhenotypes.features_tsv,
                shared_samples = PreparePhenotypes.shared_samples,
                chromosome = chromosome,
                exact_allele_counts = exact_allele_counts,
                cumulative_allele_count_maxima = cumulative_allele_count_maxima,
                consequence_classes = consequence_classes,
                maximum_gvs_maf = maximum_gvs_maf,
                maximum_distance_bp = DetermineMaximumDistance.maximum_distance_bp,
                annotation_chunk_size_bp = annotation_chunk_size_bp,
                docker_image = docker_image,
                cpu = scatter_cpu,
                memory_gb = scatter_memory_gb,
                disk_gb = scatter_disk_gb,
                max_retries = max_retries
        }
    }

    call GatherCarrierPairs {
        input:
            carrier_pairs = ClassifyChromosome.carrier_pairs_tsv,
            chromosome_qc = ClassifyChromosome.chromosome_qc_json,
            docker_image = docker_image,
            cpu = gather_cpu,
            memory_gb = gather_memory_gb,
            disk_gb = gather_disk_gb,
            max_retries = max_retries
    }

    if (publish_carrier_audit) {
        call PublishCarrierAudit {
            input:
                source_carrier_minimum_distances_tsv = GatherCarrierPairs.carrier_minimum_distances_tsv,
                docker_image = docker_image,
                cpu = gather_cpu,
                memory_gb = gather_memory_gb,
                disk_gb = gather_disk_gb,
                max_retries = max_retries
        }
    }

    call CalculateEnrichment {
        input:
            phenotype_bed = phenotype_bed,
            shared_samples = PreparePhenotypes.shared_samples,
            carrier_minimum_distances_tsv = GatherCarrierPairs.carrier_minimum_distances_tsv,
            features_tsv = PreparePhenotypes.features_tsv,
            phenotype_qc_json = PreparePhenotypes.phenotype_qc_json,
            chromosome_qc_tsv = GatherCarrierPairs.chromosome_qc_tsv,
            selected_chromosomes = chromosomes,
            exact_allele_counts = exact_allele_counts,
            cumulative_allele_count_maxima = cumulative_allele_count_maxima,
            z_thresholds = z_thresholds,
            distance_thresholds_bp = distance_thresholds_bp,
            outlier_tail = outlier_tail,
            consequence_classes = consequence_classes,
            loftee_enabled = PrepareVatIndex.loftee_enabled,
            vcf_index_provenance = PrepareVcfIndex.index_provenance,
            vat_index_provenance = PrepareVatIndex.index_provenance,
            maximum_gvs_maf = maximum_gvs_maf,
            annotation_chunk_size_bp = annotation_chunk_size_bp,
            workflow_version = workflow_version,
            docker_image = docker_image,
            cpu = gather_cpu,
            memory_gb = gather_memory_gb,
            disk_gb = gather_disk_gb,
            max_retries = max_retries
    }

    output {
        File enrichment_tsv = CalculateEnrichment.enrichment_tsv
        File enrichment_json = CalculateEnrichment.enrichment_json
        File chromosome_qc_tsv = GatherCarrierPairs.chromosome_qc_tsv
        File? carrier_minimum_distances_tsv = PublishCarrierAudit.carrier_minimum_distances_tsv
        File generated_or_validated_vat_tbi = PrepareVatIndex.vat_tbi
        File generated_or_validated_vcf_tbi = PrepareVcfIndex.vcf_tbi
        Array[File] chromosome_query_regions = ClassifyChromosome.query_regions_bed
        File phenotype_qc_json = PreparePhenotypes.phenotype_qc_json
        String vat_index_provenance = PrepareVatIndex.index_provenance
        File vat_schema_json = PrepareVatIndex.vat_schema_json
        String vcf_index_provenance = PrepareVcfIndex.index_provenance
    }
}
