version development

##########################################################################################
# A workflow that runs the Guppy basecaller on ONT FAST5 files.
# - The docker tag number will match the version of Guppy that is being run. You can change
#   this value to run a different version of Guppy. Currently supports... [3.5.2, 3.6.0, 4.0.14]
# - All fast5 files within the given GCS dir, gcs_fast5_dir, will be processed
# - Takes a few hours to process 130GB. Best guess is that the processing time scales
#   linearly but untested.
##########################################################################################

import "https://raw.githubusercontent.com/garlandcrow/helix_wdl/master/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/garlandcrow/helix_wdl/master/tasks/ONTUtils.wdl" as ONT
import "https://raw.githubusercontent.com/garlandcrow/helix_wdl/master/tasks/Structs.wdl"

workflow Guppy {
    input {
        Array [File] gcs_fast5_files

        String config
        String? barcode_kit

        String instrument = "unknown"
        String flow_cell_id = "unknown"
        String? protocol_run_id
        String? sample_name

        String gcs_out_root_dir
    }

    call Basecall {
        input:
            fast5_files  = gcs_fast5_files,
            config       = config,
            barcode_kit  = barcode_kit
    }

    call Utils.Uniq as UniqueBarcodes { input: strings = Basecall.barcodes }

    call FinalizeBasecalls {
        input:
            pass_fastqs        = Basecall.pass_fastqs,
            sequencing_summary = Basecall.sequencing_summary,
            barcodes           = UniqueBarcodes.unique_strings,
            outdir             = gcs_out_root_dir
    }

    output {
        String gcs_dir = FinalizeBasecalls.gcs_dir
        Array[String] barcodes = UniqueBarcodes.unique_strings
        Int num_pass_fastqs = Basecall.num_pass_fastqs
        Int num_fail_fastqs = Basecall.num_fail_fastqs
    }
}

task Basecall {
    input {
        Array[File] fast5_files
        String config = "dna_r9.4.1_450bps_hac_prom.cfg"
        String? barcode_kit

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(fast5_files, "GB"))

    String barcode_arg = if defined(barcode_kit) then "--barcode_kits \"~{barcode_kit}\" --trim_barcodes" else ""

    command <<<
        set -x

        ls -l /cromwell_root/ > input_files_list.txt

        guppy_basecaller \
            -r \
            -i /cromwell_root/ \
            -s guppy_output/ \
            -x "cuda:all" \
            -c ~{config} \
            ~{barcode_arg} \
            --compress_fastq

        # Make a list of the barcodes that were seen in the data
        find guppy_output/ -name '*fastq*' -not -path '*fail*' -type f | \
            awk -F"/" '{ a=NF-1; a=$a; gsub(/pass/, "unclassified", a); print a }' | \
            sort -n | \
            uniq | tee barcodes.txt

        # Reorganize and rename the passing filter data to include the barcode in the filename
        mkdir pass
        find guppy_output/ -name '*fastq*' -not -path '*fail*' -type f | \
            awk -F"/" '{ a=NF-1; a=$a; b=$NF; gsub(/pass/, "unclassified", a); c=$NF; for (i = NF-1; i > 0; i--) { c=$i"/"c }; system("mv " c " pass/" a "." b); }'

        # Reorganize and rename the failing filter data to include the barcode in the filename
        mkdir fail
        find guppy_output/ -name '*fastq*' -not -path '*pass*' -type f | \
            awk -F"/" '{ a=NF-1; a=$a; b=$NF; gsub(/pass/, "unclassified", a); c=$NF; for (i = NF-1; i > 0; i--) { c=$i"/"c }; system("mv " c " fail/" a "." b); }'

        # Count passing and failing files
        find pass -name '*fastq.gz' | wc -l | tee num_pass.txt
        find fail -name '*fastq.gz' | wc -l | tee num_fail.txt

        # run barcoder 
        guppy_barcoder \
            --recursive \
            --input_path guppy_output \
            --save_path guppy_barcode_output \
            ~{barcode_arg} \
            --device "cuda:all"
    >>>

    output {
        File input_files_list = "input_files_list.txt"
        Array[File] pass_fastqs = glob("pass/*.fastq.gz")
        File sequencing_summary = "guppy_output/sequencing_summary.txt"
        Array[String] barcodes = read_lines("barcodes.txt")

        Directory guppy_barcode_output = "guppy_barcode_output"

        Int num_pass_fastqs = read_int("num_pass.txt")
        Int num_fail_fastqs = read_int("num_fail.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-guppy:4.5.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-central1-c", "us-central1-f", "us-east1-b", "us-east1-c", "us-west1-a", "us-west1-b"]
        cpuPlatform:            "Intel Haswell"
    }
}

task MakeSequencingSummary {
    input {
        Array[File] sequencing_summaries

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(sequencing_summaries, "GB"))

    command <<<
        set -euxo pipefail

        head -1 ~{sequencing_summaries[0]} > sequencing_summary.txt

        while read p; do
            awk 'NR > 1 { print }' "$p" >> sequencing_summary.txt
        done <~{write_lines(sequencing_summaries)}
    >>>

    output {
        File sequencing_summary = "sequencing_summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FinalizeBasecalls {
    input {
        Array[String] pass_fastqs
        File sequencing_summary
        File final_summary
        Array[String] barcodes

        String outdir

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir + "/", "/+$", "")

    command <<<
        set -euxo pipefail

        PASS_FASTQ="~{write_lines(pass_fastqs)}"

        while read b; do
            OUT_DIR="~{gcs_output_dir}/$b"
            PASS_DIR="$OUT_DIR/fastq_pass/"

            grep -w $b $PASS_FASTQ | gsutil -m cp -I $PASS_DIR

            if [ ~{length(barcodes)} -eq 1 ]; then
                cp ~{sequencing_summary} sequencing_summary.$b.txt
                cp ~{final_summary} final_summary.$b.txt
            else
                grep -w -e filename -e $b ~{sequencing_summary} > sequencing_summary.$b.txt
                sed "s/sample_id=/sample_id=$b./" ~{final_summary} > final_summary.$b.txt
            fi

            gsutil cp sequencing_summary.$b.txt $OUT_DIR/
            gsutil cp final_summary.$b.txt $OUT_DIR/
        done <~{write_lines(barcodes)}
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}