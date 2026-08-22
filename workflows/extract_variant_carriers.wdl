version 1.0

task PrepareCarrierInputs {
    input {
        File filtered_vcf
        File filtered_vcf_tbi
        File transcript_annotations
        File? transcript_annotations_tbi
        Array[String] chromosomes
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File chromosome_file = write_lines(chromosomes)
    String transcript_index_source = if defined(transcript_annotations_tbi) then "~{select_first([transcript_annotations_tbi])}" else ""

    command <<<
        set -euo pipefail

        echo "Starting carrier input preparation"
        echo "Input byte count: VCF=$(wc -c < '~{filtered_vcf}'), transcripts=$(wc -c < '~{transcript_annotations}')"

        ln -s "~{filtered_vcf}" filtered_variants.vcf.gz
        ln -s "~{filtered_vcf_tbi}" filtered_variants.vcf.gz.tbi
        ln -s "~{transcript_annotations}" transcript_annotations.tsv.bgz

        if [[ -n "~{transcript_index_source}" ]]; then
            ln -s "~{transcript_index_source}" transcript_annotations.tsv.bgz.tbi
            transcript_index_provenance="supplied"
        else
            tabix -f -S 1 -s 1 -b 2 -e 2 transcript_annotations.tsv.bgz
            transcript_index_provenance="generated"
        fi
        printf '%s\n' "$transcript_index_provenance" > transcript_index_provenance.txt

        requested_count=$(grep -cve '^[[:space:]]*$' "~{chromosome_file}")
        while IFS= read -r chromosome; do
            [[ -z "$chromosome" ]] && continue
            tabix -l filtered_variants.vcf.gz | grep -Fxq "$chromosome"
            tabix -l transcript_annotations.tsv.bgz | grep -Fxq "$chromosome"
        done < "~{chromosome_file}"
        chromosomes_csv=$(paste -sd, "~{chromosome_file}")
        echo "Selected chromosome count: $requested_count"
        echo "Transcript index provenance: $transcript_index_provenance"

        rare-variant-enrichment prepare-carrier-inputs \
            --vcf filtered_variants.vcf.gz \
            --annotations transcript_annotations.tsv.bgz \
            --chromosomes "$chromosomes_csv" \
            --vcf-index-provenance supplied \
            --transcript-index-provenance "$transcript_index_provenance" \
            --schema-output transcript.schema.json \
            --qc-output transcript.prepare.qc.json

        echo "Completed carrier input preparation; selected chromosome count: $requested_count"
    >>>

    output {
        File prepared_vcf = filtered_vcf
        File prepared_vcf_tbi = filtered_vcf_tbi
        File prepared_transcript_annotations = transcript_annotations
        File prepared_transcript_annotations_tbi = if defined(transcript_annotations_tbi) then select_first([transcript_annotations_tbi]) else "transcript_annotations.tsv.bgz.tbi"
        File transcript_schema_json = "transcript.schema.json"
        File preparation_qc_json = "transcript.prepare.qc.json"
        String transcript_index_provenance = read_string("transcript_index_provenance.txt")
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task ExtractChromosomeCarriers {
    input {
        File filtered_vcf
        File filtered_vcf_tbi
        File transcript_annotations
        File transcript_annotations_tbi
        File transcript_schema_json
        String chromosome
        Int annotation_chunk_size_bp
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
        Int preemptible
    }

    command <<<
        set -euo pipefail

        echo "Starting chromosome carrier extraction: ~{chromosome}"
        echo "Chunk base-pair count: ~{annotation_chunk_size_bp}"
        ln -s "~{filtered_vcf}" filtered_variants.vcf.gz
        ln -s "~{filtered_vcf_tbi}" filtered_variants.vcf.gz.tbi
        ln -s "~{transcript_annotations}" transcript_annotations.tsv.bgz
        ln -s "~{transcript_annotations_tbi}" transcript_annotations.tsv.bgz.tbi

        rare-variant-enrichment extract-gene-carriers \
            --vcf filtered_variants.vcf.gz \
            --annotations transcript_annotations.tsv.bgz \
            --schema "~{transcript_schema_json}" \
            --chromosome "~{chromosome}" \
            --chunk-size-bp "~{annotation_chunk_size_bp}" \
            --audit-output chromosome.audit.tsv.gz \
            --qc-output chromosome.qc.json

        audit_count=$(gzip -cd chromosome.audit.tsv.gz | tail -n +2 | wc -l | tr -d ' ')
        echo "Completed chromosome carrier extraction: ~{chromosome}; audit row count: $audit_count"
    >>>

    output {
        File audit_tsv_gz = "chromosome.audit.tsv.gz"
        File qc_json = "chromosome.qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
        preemptible: preemptible
    }
}

task GatherVariantCarriers {
    input {
        Array[File] audit_inputs
        Array[File] qc_inputs
        File preparation_qc_json
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    command <<<
        set -euo pipefail

        echo "Starting carrier gather; shard count: ~{length(audit_inputs)}"
        rare-variant-enrichment gather-gene-carriers \
            --audit-input "~{sep="\" --audit-input \"" audit_inputs}" \
            --qc-input "~{sep="\" --qc-input \"" qc_inputs}" \
            --preparation-qc "~{preparation_qc_json}" \
            --audit-output variant_carrier_audit.tsv.gz \
            --carrier-output variant_carriers.tsv.gz \
            --qc-output variant_carriers.qc.json

        audit_count=$(gzip -cd variant_carrier_audit.tsv.gz | tail -n +2 | wc -l | tr -d ' ')
        carrier_count=$(gzip -cd variant_carriers.tsv.gz | tail -n +2 | wc -l | tr -d ' ')
        echo "Completed carrier gather; audit row count: $audit_count; carrier row count: $carrier_count"
    >>>

    output {
        File audit_tsv_gz = "variant_carrier_audit.tsv.gz"
        File carrier_tsv_gz = "variant_carriers.tsv.gz"
        File qc_json = "variant_carriers.qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

workflow ExtractVariantCarriers {
    input {
        File filtered_vcf
        File filtered_vcf_tbi
        File transcript_annotations
        File? transcript_annotations_tbi
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
        Int annotation_chunk_size_bp = 10000000
        String docker_image = "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main"
        Int prepare_cpu = 2
        Int prepare_memory_gb = 8
        Int prepare_disk_gb = 50
        Int scatter_cpu = 2
        Int scatter_memory_gb = 16
        Int scatter_disk_gb = 50
        Int gather_cpu = 2
        Int gather_memory_gb = 32
        Int gather_disk_gb = 100
        Int max_retries = 1
        Int scatter_preemptible = 2
    }

    Int calculated_prepare_disk_gb = ceil((size(filtered_vcf, "GiB") + size(transcript_annotations, "GiB")) * 2.0 + 20.0)
    Int dynamic_prepare_disk_gb = if calculated_prepare_disk_gb > prepare_disk_gb then calculated_prepare_disk_gb else prepare_disk_gb

    call PrepareCarrierInputs {
        input:
            filtered_vcf = filtered_vcf,
            filtered_vcf_tbi = filtered_vcf_tbi,
            transcript_annotations = transcript_annotations,
            transcript_annotations_tbi = transcript_annotations_tbi,
            chromosomes = chromosomes,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = dynamic_prepare_disk_gb,
            max_retries = max_retries
    }

    Int calculated_scatter_disk_gb = ceil((size(filtered_vcf, "GiB") + size(transcript_annotations, "GiB")) * 2.0 + 20.0)
    Int dynamic_scatter_disk_gb = if calculated_scatter_disk_gb > scatter_disk_gb then calculated_scatter_disk_gb else scatter_disk_gb

    scatter (chromosome in chromosomes) {
        call ExtractChromosomeCarriers {
            input:
                filtered_vcf = PrepareCarrierInputs.prepared_vcf,
                filtered_vcf_tbi = PrepareCarrierInputs.prepared_vcf_tbi,
                transcript_annotations = PrepareCarrierInputs.prepared_transcript_annotations,
                transcript_annotations_tbi = PrepareCarrierInputs.prepared_transcript_annotations_tbi,
                transcript_schema_json = PrepareCarrierInputs.transcript_schema_json,
                chromosome = chromosome,
                annotation_chunk_size_bp = annotation_chunk_size_bp,
                docker_image = docker_image,
                cpu = scatter_cpu,
                memory_gb = scatter_memory_gb,
                disk_gb = dynamic_scatter_disk_gb,
                max_retries = max_retries,
                preemptible = scatter_preemptible
        }
    }

    Int calculated_gather_disk_gb = ceil((size(ExtractChromosomeCarriers.audit_tsv_gz, "GiB") + size(ExtractChromosomeCarriers.qc_json, "GiB")) * 2.0 + 20.0)
    Int dynamic_gather_disk_gb = if calculated_gather_disk_gb > gather_disk_gb then calculated_gather_disk_gb else gather_disk_gb

    call GatherVariantCarriers {
        input:
            audit_inputs = ExtractChromosomeCarriers.audit_tsv_gz,
            qc_inputs = ExtractChromosomeCarriers.qc_json,
            preparation_qc_json = PrepareCarrierInputs.preparation_qc_json,
            docker_image = docker_image,
            cpu = gather_cpu,
            memory_gb = gather_memory_gb,
            disk_gb = dynamic_gather_disk_gb,
            max_retries = max_retries
    }

    output {
        File variant_carrier_audit_tsv_gz = GatherVariantCarriers.audit_tsv_gz
        File variant_carriers_tsv_gz = GatherVariantCarriers.carrier_tsv_gz
        File variant_carriers_qc_json = GatherVariantCarriers.qc_json
        Array[File] chromosome_qc_jsons = ExtractChromosomeCarriers.qc_json
        File transcript_schema_json = PrepareCarrierInputs.transcript_schema_json
        File generated_or_validated_transcript_annotations_tbi = PrepareCarrierInputs.prepared_transcript_annotations_tbi
        String transcript_index_provenance = PrepareCarrierInputs.transcript_index_provenance
    }
}
