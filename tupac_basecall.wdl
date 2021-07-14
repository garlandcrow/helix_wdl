version 1.0

import "https://raw.githubusercontent.com/garlandcrow/helix_wdl/master/tasks/Guppy.wdl" as Guppy

workflow TupacBasecall {
    input {
        String gcs_fast5_dir
        String config = "dna_r9.4.1_450bps_hac.cfg"
        String? barcode_kit
        String gcs_out_root_dir
        String prefix
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/TupacBasecall/~{prefix}"

    call Guppy.Guppy {
        input:
            gcs_fast5_dir    = gcs_fast5_dir,
            config           = config,
            barcode_kit      = barcode_kit,
            gcs_out_root_dir = outdir
    }

    output {
        String gcs_basecall_dir = Guppy.gcs_dir
    }
}

