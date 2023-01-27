#!/usr/bin/env cwl-runner
{
    "$graph": [
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "hints": [
                {
                    "class": "ResourceRequirement",
                    "coresMin": 4,
                    "ramMin": 18000
                }
            ],
            "baseCommand": [
                "mist_add_to_bam.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--bamIO-threads",
                        "valueFrom": "$(runtime.coresMin)"
                    },
                    "id": "#AddtoBam.cwl/Bam_IO_Threads"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--cell-order"
                    },
                    "id": "#AddtoBam.cwl/Cell_Order"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#AddtoBam.cwl/Generate_Bam"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--annot-mol-file"
                    },
                    "id": "#AddtoBam.cwl/Molecular_Annotation"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--r2-bam"
                    },
                    "id": "#AddtoBam.cwl/R2_Bam"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#AddtoBam.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--tag-calls"
                    },
                    "id": "#AddtoBam.cwl/Tag_Calls"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--target-gene-mapping"
                    },
                    "id": "#AddtoBam.cwl/Target_Gene_Mapping"
                }
            ],
            "id": "#AddtoBam.cwl",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "Annotated_mapping_R2.BAM"
                    },
                    "id": "#AddtoBam.cwl/Annotated_Bam"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AddtoBam.cwl/output"
                }
            ]
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "doc": "AlignmentAnalysis stage of the Rain pipeline annotates aligned reads and collects a myriad of metrics on the aligned reads. Additional annotation is performed to the reads\n",
            "hints": [
                {
                    "class": "ResourceRequirement",
                    "coresMin": 8,
                    "ramMin": 24000
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "AlignmentAnalysisAndCountCB.sh"
            ],
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--exclude-intronic-reads"
                    },
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Exclude_Intronic_Reads"
                },
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/GTF"
                },
                {
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Maximum_Threads"
                },
                {
                    "inputBinding": {
                        "prefix": "--r2-bam",
                        "itemSeparator": ","
                    },
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#AlignmentAnalysis.cwl/R2_BAM"
                },
                {
                    "inputBinding": {
                        "prefix": "--quality-metrics"
                    },
                    "type": "File",
                    "id": "#AlignmentAnalysis.cwl/ReadQualityMetrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AlignmentAnalysis.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--transcript-length"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Transcript_Length"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*.annotated.*.bam"
                    },
                    "id": "#AlignmentAnalysis.cwl/Annotated_Bam_Files"
                },
                {
                    "outputBinding": {
                        "glob": "*logs.tar.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Logs"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_SeqMetrics.csv"
                    },
                    "id": "#AlignmentAnalysis.cwl/Seq_Metrics"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*Sorted_Valid_Reads.csv.*"
                    },
                    "id": "#AlignmentAnalysis.cwl/Sorted_Valid_Reads_CSV"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "count_estimates.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 30000; } return parseInt(JSON.parse(self[0].contents).num_bioproducts); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_bioproducts"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "count_estimates.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 10000; } return parseInt(JSON.parse(self[0].contents).num_cell_estimate); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_cell_estimate"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).BCR); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_valid_ig_reads"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).TCR); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_valid_tcr_reads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_BCR_Valid_Reads.fastq.gz"
                    },
                    "id": "#AlignmentAnalysis.cwl/validIgReads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_TCR_Valid_Reads.fastq.gz"
                    },
                    "id": "#AlignmentAnalysis.cwl/validTcrReads"
                }
            ],
            "id": "#AlignmentAnalysis.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_annotate_molecules.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--umi-option"
                    },
                    "id": "#AnnotateMolecules.cwl/AbSeq_UMI"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#AnnotateMolecules.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "prefix": "--use-dbec"
                    },
                    "id": "#AnnotateMolecules.cwl/Use_DBEC"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--valid-annot"
                    },
                    "id": "#AnnotateMolecules.cwl/Valids"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_GeneStatus.csv.*"
                    },
                    "id": "#AnnotateMolecules.cwl/Gene_Status_List"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "stats.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).max_count)\n"
                    },
                    "id": "#AnnotateMolecules.cwl/Max_Count"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_Annotation_Molecule.csv.*"
                    },
                    "id": "#AnnotateMolecules.cwl/Mol_Annot_List"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "stats.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).total_molecules)\n"
                    },
                    "id": "#AnnotateMolecules.cwl/Total_Molecules"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AnnotateMolecules.cwl/output"
                }
            ],
            "id": "#AnnotateMolecules.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "id": "#Assay_Settings.cwl/_Assay"
                }
            ],
            "outputs": [
                {
                    "type": "string",
                    "id": "#Assay_Settings.cwl/Assay"
                }
            ],
            "expression": "${\n  return ({Assay: inputs._Assay})\n}",
            "id": "#Assay_Settings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#BamSettings.cwl/_Generate_Bam"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#BamSettings.cwl/Generate_Bam"
                }
            ],
            "expression": "${\n  // the create bam flag defaults to false\n  var generateBam = false;\n  // the user can set this flag to true, to enable creation of the bam file.\n  if (inputs._Generate_Bam != null) {\n    generateBam = inputs._Generate_Bam;\n  }\n  return ({\n    Generate_Bam: generateBam,\n  });\n}",
            "id": "#BamSettings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#BundleLogs.cwl/log_files"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#BundleLogs.cwl/logs_dir"
                }
            ],
            "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    /*\n      Checking here for null in case a Node was skipped because of conditional execution.\n      For e.g. Generate_Bam is used to skip the AddToBam, MergeBam and IndexBam nodes\n    */\n    if (log_file != null) {\n      log_file.basename = uuid() + \"-\" + log_file.basename;\n      listing.push(log_file);\n    }\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
            "id": "#BundleLogs.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_check_references.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--abseq-reference"
                    },
                    "id": "#CheckReference.cwl/AbSeq_Reference"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "prefix": "--assay"
                    },
                    "id": "#CheckReference.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "id": "#CheckReference.cwl/Putative_Cell_Call"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--reference"
                    },
                    "id": "#CheckReference.cwl/Reference"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--sample-tags-version"
                    },
                    "id": "#CheckReference.cwl/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--supplemental-reference"
                    },
                    "id": "#CheckReference.cwl/Supplemental_Reference"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-version"
                    },
                    "id": "#CheckReference.cwl/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "combined_extra_seq.fasta"
                    },
                    "id": "#CheckReference.cwl/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "full-gene-list.json"
                    },
                    "id": "#CheckReference.cwl/Full_Genes"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*gtf",
                        "outputEval": "${\n    // get the WTA modified GTF with extra seqs\n    if (self.length == 1) {\n        return self;\n    // there is no modified GTF\n    } else if (self.length == 0) {\n        // if Reference is null (i.e. AbSeq_Reference only), return no GTF\n        if (inputs.Reference === null) {\n            return null;\n        } else {\n            // get the original WTA GTF without extra seqs\n            for (var i = 0; i < inputs.Reference.length; i++) {\n                if (inputs.Reference[i].basename.toLowerCase().indexOf('gtf') !== -1) {\n                    return inputs.Reference[i];\n                }\n            }\n            // return no GTF for Targeted\n            return null\n        }\n    }\n}\n"
                    },
                    "id": "#CheckReference.cwl/GTF"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*-annot.*",
                        "outputEval": "${\n    if (self.length == 1) { // Targeted\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or targets\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('tar.gz') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
                    },
                    "id": "#CheckReference.cwl/Index"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "target-gene.json"
                    },
                    "id": "#CheckReference.cwl/Target_Gene_Mapping"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "transcript_length.json"
                    },
                    "id": "#CheckReference.cwl/Transcript_Length"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#CheckReference.cwl/output"
                }
            ],
            "id": "#CheckReference.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_dense_to_sparse.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--cell-order"
                    },
                    "id": "#DensetoSparse.cwl/Cell_Order"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--dense-data-table"
                    },
                    "id": "#DensetoSparse.cwl/Dense_Data_Table"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--gene-list"
                    },
                    "id": "#DensetoSparse.cwl/Gene_List"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#DensetoSparse.cwl/Run_Metadata"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.csv.gz"
                    },
                    "id": "#DensetoSparse.cwl/Data_Tables"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#DensetoSparse.cwl/output"
                }
            ],
            "id": "#DensetoSparse.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": "cat",
            "stdout": "cell_order.json",
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#DensetoSparseFile.cwl/GDT_cell_order"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "cell_order.json"
                    },
                    "id": "#DensetoSparseFile.cwl/Cell_Order"
                }
            ],
            "id": "#DensetoSparseFile.cwl"
        },
        {
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 48000
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:2.0b2"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_get_datatables.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--full-gene-list"
                    },
                    "id": "#GetDataTable.cwl/Full_Genes"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--gene-status"
                    },
                    "id": "#GetDataTable.cwl/Gene_Status_List"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--max-count"
                    },
                    "id": "#GetDataTable.cwl/Max_Count"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--mol-annot"
                    },
                    "id": "#GetDataTable.cwl/Molecule_Annotation_List"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "id": "#GetDataTable.cwl/Putative_Cell_Call"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#GetDataTable.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--tag-names"
                    },
                    "id": "#GetDataTable.cwl/Sample_Tag_Names"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--seq-metrics"
                    },
                    "id": "#GetDataTable.cwl/Seq_Metrics"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "id": "#GetDataTable.cwl/Total_Molecules"
                },
                {
                    "type": "int",
                    "id": "#GetDataTable.cwl/num_bioproducts"
                },
                {
                    "type": "int",
                    "id": "#GetDataTable.cwl/num_cell_estimate"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "metrics-files.tar.gz"
                    },
                    "id": "#GetDataTable.cwl/Annot_Files"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Annotations/*_Bioproduct_Stats.csv"
                    },
                    "id": "#GetDataTable.cwl/Bioproduct_Stats"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*_Cell_Label_Graph_Data.json"
                    },
                    "id": "#GetDataTable.cwl/Cell_Label_Graph_Data"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "cell_order.json"
                    },
                    "id": "#GetDataTable.cwl/Cell_Order"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*cell_type_experimental.csv"
                    },
                    "id": "#GetDataTable.cwl/Cell_Type_Predictions"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Annotation_Molecule_corrected.csv.gz"
                    },
                    "id": "#GetDataTable.cwl/Corrected_Molecular_Annotation"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*PerCell_Dense.csv.gz"
                    },
                    "id": "#GetDataTable.cwl/Dense_Data_Tables"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*PerCell_Unfiltered_Dense.csv.gz"
                    },
                    "id": "#GetDataTable.cwl/Dense_Data_Tables_Unfiltered"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_coordinates.csv"
                    },
                    "id": "#GetDataTable.cwl/Dim_Reduction_Coord"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_Expression_Data.st.gz"
                    },
                    "id": "#GetDataTable.cwl/Expression_Data"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Expression_Data_Unfiltered.st.gz"
                    },
                    "id": "#GetDataTable.cwl/Expression_Data_Unfiltered"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "gene_list.json"
                    },
                    "id": "#GetDataTable.cwl/Gene_List"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "Annotations/*_Annotation_Molecule.csv.gz"
                    },
                    "id": "#GetDataTable.cwl/Molecular_Annotation"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*_Protein_Aggregates_Experimental.csv"
                    },
                    "id": "#GetDataTable.cwl/Protein_Aggregates_Experimental"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*_Putative_Cells_Origin.csv"
                    },
                    "id": "#GetDataTable.cwl/Putative_Cells_Origin"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "SampleTag/*csv"
                    },
                    "id": "#GetDataTable.cwl/SampleTag_out"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "SampleTag/*zip"
                    },
                    "id": "#GetDataTable.cwl/SampleTag_zip"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Annotations/*_Annotation_Molecule_SampleTag.csv"
                    },
                    "id": "#GetDataTable.cwl/Tag_Annotation"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "SampleTag/*_Calls.csv"
                    },
                    "id": "#GetDataTable.cwl/Tag_Calls"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Annotations/*_UMI_Adjusted_CellLabel_Stats.csv"
                    },
                    "id": "#GetDataTable.cwl/UMI_Adjusted_CellLabel_Stats"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Visual_Metrics_MolsPerCell.csv.gz"
                    },
                    "id": "#GetDataTable.cwl/Visual_Metrics_MolsPerCell"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#GetDataTable.cwl/output"
                }
            ],
            "id": "#GetDataTable.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "samtools",
                "index"
            ],
            "stdout": "samtools_index.log",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#IndexBAM.cwl/BamFile"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#IndexBAM.cwl/Generate_Bam"
                }
            ],
            "arguments": [
                {
                    "position": 2,
                    "valueFrom": "${\n    return inputs.BamFile.basename + \".bai\"\n}"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.bai"
                    },
                    "id": "#IndexBAM.cwl/Index"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#IndexBAM.cwl/log"
                }
            ],
            "id": "#IndexBAM.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/_AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/_Barcode_Num"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#InternalSettings.cwl/_Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/_MinChunkSize"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#InternalSettings.cwl/_NumRecordsPerSplit"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/_Read_Filter_Off"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#InternalSettings.cwl/_Seq_Run"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/_Subsample_Sample_Tags"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/_Target_analysis"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/_Use_DBEC"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/_VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/_VDJ_VGene_Evalue"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/Barcode_Num"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#InternalSettings.cwl/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/Label_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/MinChunkSize"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#InternalSettings.cwl/NumRecordsPerSplit"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Read_Filter_Off"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#InternalSettings.cwl/Seq_Run"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/Subsample_Sample_Tags"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Target_analysis"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Use_DBEC"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/VDJ_VGene_Evalue"
                }
            ],
            "expression": "${\n  var internalInputs = [\n    '_Label_Version',\n    '_Read_Filter_Off',\n    '_Barcode_Num',\n    '_Seq_Run',\n    '_AbSeq_UMI',\n    '_Use_DBEC',\n    '_Extra_Seqs',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n    '_Target_analysis',\n    '_Subsample_Sample_Tags',\n    '_VDJ_VGene_Evalue',\n    '_VDJ_JGene_Evalue',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
            "id": "#InternalSettings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#IntronicReadsSettings.cwl/_Exclude_Intronic_Reads"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#IntronicReadsSettings.cwl/Exclude_Intronic_Reads"
                }
            ],
            "expression": "${\n  // the exclude intronic reads flag defaults to false\n  var excludeIntronicReads = false;\n  // the user can set the flag to exclude intronic reads\n  if (inputs._Exclude_Intronic_Reads) {\n    excludeIntronicReads = inputs._Exclude_Intronic_Reads;\n  }\n  return ({\n    Exclude_Intronic_Reads: excludeIntronicReads,\n  });\n}",
            "id": "#IntronicReadsSettings.cwl"
        },
        {
            "class": "Workflow",
            "label": "BD Rhapsody\u2122 WTA Analysis Pipeline - 2.0beta2",
            "doc": "The BD Rhapsody\u2122 WTA Analysis Pipeline is used to create sequencing libraries from single cell transcriptomes without having to specify a targeted panel.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files, a reference genome file and a transcriptome annotation file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "AbSeq Reference",
                    "id": "#main/AbSeq_Reference"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/Barcode_Num"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Disable Refined Putative Cell Calling",
                    "doc": "Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.",
                    "id": "#main/Basic_Algo_Only"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Exact Cell Count",
                    "doc": "Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count",
                    "id": "#main/Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Exclude Intronic Reads",
                    "doc": "By default, reads aligned to exons and introns are considered and represented in molecule counts. Including intronic reads may increase sensitivity, resulting in an increase in molecule counts and the number of genes per cell for both cellular and nuclei samples. Intronic reads may indicate unspliced mRNAs and are also useful, for example, in the study of nuclei and RNA velocity. When set to true, intronic reads will be excluded.",
                    "id": "#main/Exclude_Intronic_Reads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Generate Bam Output",
                    "doc": "Default: false.  A Bam read alignment file contains reads from all the input libraries, but creating it can consume a lot of compute and disk resources. By setting this field to true, the Bam file will be created.\n",
                    "id": "#main/Generate_Bam"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Maximum Number of Threads",
                    "doc": "The maximum number of threads to use in the pipeline. By default, all available cores are used.",
                    "id": "#main/Maximum_Threads"
                },
                {
                    "label": "Putative Cell Calling",
                    "doc": "Specify the data to be used for putative cell calling. mRNA is the default selected option.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Putative_Cell_Call/Putative_Cell_Call",
                            "symbols": [
                                "#main/Putative_Cell_Call/Putative_Cell_Call/mRNA",
                                "#main/Putative_Cell_Call/Putative_Cell_Call/AbSeq_Experimental",
                                "#main/Putative_Cell_Call/Putative_Cell_Call/AbSeq"
                            ]
                        }
                    ],
                    "id": "#main/Putative_Cell_Call"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/Read_Filter_Off"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "label": "Reads",
                    "id": "#main/Reads"
                },
                {
                    "type": "File",
                    "label": "Reference Genome",
                    "id": "#main/Reference_Genome"
                },
                {
                    "label": "Run Name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "This is a name for output files, for example Experiment1_Metrics_Summary.csv. Default if left empty is to name run based on a library. Any non-alpha numeric characters will be changed to a hyphen.",
                    "id": "#main/Run_Name"
                },
                {
                    "label": "Sample Tags Version",
                    "doc": "The sample multiplexing kit version.  This option should only be set for a multiplexed experiment.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Sample_Tags_Version/Sample_Tags_Version",
                            "symbols": [
                                "#main/Sample_Tags_Version/Sample_Tags_Version/human",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/hs",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mouse",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mm",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/flex",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/custom"
                            ]
                        }
                    ],
                    "id": "#main/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#main/Seq_Run"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Supplemental Reference",
                    "id": "#main/Supplemental_Reference"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "label": "Sample Tag Names",
                    "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Should be alpha numeric, with + - and _ allowed. Any special characters: &, (), [], {}, <>, ?, | will be corrected to underscores. \n",
                    "id": "#main/Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/Target_analysis"
                },
                {
                    "type": "File",
                    "label": "Transcriptome Annotation",
                    "id": "#main/Transcriptome_Annotation"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/Use_DBEC"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "e-value threshold for J gene",
                    "doc": "The e-value threshold for J gene call by IgBlast/PyIR, default is set as 0.001\n",
                    "id": "#main/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "e-value threshold for V gene",
                    "doc": "The e-value threshold for V gene call by IgBlast/PyIR, default is set as 0.001\n",
                    "id": "#main/VDJ_VGene_Evalue"
                },
                {
                    "label": "VDJ Species Version",
                    "doc": "The VDJ species and chain types.  This option should only be set for VDJ experiment.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/VDJ_Version/VDJ_Version",
                            "symbols": [
                                "#main/VDJ_Version/VDJ_Version/human",
                                "#main/VDJ_Version/VDJ_Version/hs",
                                "#main/VDJ_Version/VDJ_Version/mouse",
                                "#main/VDJ_Version/VDJ_Version/mm",
                                "#main/VDJ_Version/VDJ_Version/humanBCR",
                                "#main/VDJ_Version/VDJ_Version/humanTCR",
                                "#main/VDJ_Version/VDJ_Version/mouseBCR",
                                "#main/VDJ_Version/VDJ_Version/mouseTCR"
                            ]
                        }
                    ],
                    "id": "#main/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "label": "BAM File",
                    "type": "File",
                    "outputSource": "#main/MergeBAM/Bam",
                    "id": "#main/Bam"
                },
                {
                    "label": "BAM Index",
                    "type": "File",
                    "outputSource": "#main/IndexBAM/Index",
                    "id": "#main/Bam_Index"
                },
                {
                    "label": "Bioproduct Statistics",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Bioproduct_Stats",
                    "id": "#main/Bioproduct_Stats"
                },
                {
                    "label": "Data Tables",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/Uncompress_Datatables/Uncompressed_Data_Tables",
                    "id": "#main/Data_Tables"
                },
                {
                    "label": "Unfiltered Data Tables",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/Dense_to_Sparse_Datatable_Unfiltered/Data_Tables",
                    "id": "#main/Data_Tables_Unfiltered"
                },
                {
                    "label": "Dimensionality Reduction Coordinates",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Dim_Reduction_Coord",
                    "id": "#main/Dim_Reduction_Coord"
                },
                {
                    "label": "Expression Matrix",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/Uncompress_Datatables/Uncompressed_Expression_Matrix",
                    "id": "#main/Expression_Data"
                },
                {
                    "label": "Unfiltered Expression Matrix",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Expression_Data_Unfiltered",
                    "id": "#main/Expression_Data_Unfiltered"
                },
                {
                    "label": "Immune Cell Classification (Experimental)",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Cell_Type_Predictions",
                    "id": "#main/Immune_Cell_Classification(Experimental)"
                },
                {
                    "label": "Pipeline Logs",
                    "type": "Directory",
                    "outputSource": "#main/BundleLogs/logs_dir",
                    "id": "#main/Logs"
                },
                {
                    "label": "Metrics Summary",
                    "type": "File",
                    "outputSource": "#main/Metrics/Metrics_Summary",
                    "id": "#main/Metrics_Summary"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/MergeMultiplex/Multiplex_out",
                    "id": "#main/Multiplex"
                },
                {
                    "label": "Protein Aggregates (Experimental)",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Protein_Aggregates_Experimental",
                    "id": "#main/Protein_Aggregates_Experimental"
                },
                {
                    "label": "Putative Cells Origin",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Putative_Cells_Origin",
                    "id": "#main/Putative_Cells_Origin"
                },
                {
                    "label": "Pipeline Report HTML",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/Metrics/Visual_Metrics_html",
                    "id": "#main/Visual_Metrics_html"
                },
                {
                    "label": "vdjCellsDatatable",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjCellsDatatable",
                    "id": "#main/vdjCellsDatatable"
                },
                {
                    "label": "vdjCellsDatatableUncorrected",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjCellsDatatableUncorrected",
                    "id": "#main/vdjCellsDatatableUncorrected"
                },
                {
                    "label": "vdjDominantContigsAIRR",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjDominantContigsAIRR",
                    "id": "#main/vdjDominantContigsAIRR"
                },
                {
                    "label": "vdjMetricsCsv",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjMetricsCsv",
                    "id": "#main/vdjMetricsCsv"
                },
                {
                    "label": "vdjUnfilteredContigsAIRR",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjUnfilteredContigsAIRR",
                    "id": "#main/vdjUnfilteredContigsAIRR"
                }
            ],
            "steps": [
                {
                    "run": "#AddtoBam.cwl",
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_File/Cell_Order",
                            "id": "#main/AddtoBam/Cell_Order"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/AddtoBam/Generate_Bam"
                        },
                        {
                            "source": "#main/GetDataTable/Corrected_Molecular_Annotation",
                            "id": "#main/AddtoBam/Molecular_Annotation"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Annotated_Bam_Files",
                            "id": "#main/AddtoBam/R2_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AddtoBam/Run_Metadata"
                        },
                        {
                            "source": "#main/GetDataTable/Tag_Calls",
                            "id": "#main/AddtoBam/Tag_Calls"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AddtoBam/Target_Gene_Mapping"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true )",
                    "out": [
                        "#main/AddtoBam/Annotated_Bam",
                        "#main/AddtoBam/output"
                    ],
                    "scatter": [
                        "#main/AddtoBam/R2_Bam"
                    ],
                    "id": "#main/AddtoBam"
                },
                {
                    "run": "#AlignmentAnalysis.cwl",
                    "in": [
                        {
                            "source": "#main/Intronic_Reads_Settings/Exclude_Intronic_Reads",
                            "id": "#main/AlignmentAnalysis/Exclude_Intronic_Reads"
                        },
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AlignmentAnalysis/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/AlignmentAnalysis/GTF"
                        },
                        {
                            "source": "#main/Maximum_Threads",
                            "id": "#main/AlignmentAnalysis/Maximum_Threads"
                        },
                        {
                            "source": "#main/QualCLAlign/R2BAMFiles",
                            "id": "#main/AlignmentAnalysis/R2_BAM"
                        },
                        {
                            "source": "#main/QualCLAlign/QualCLAlignMetrics",
                            "id": "#main/AlignmentAnalysis/ReadQualityMetrics"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AlignmentAnalysis/Run_Metadata"
                        },
                        {
                            "source": "#main/CheckReference/Transcript_Length",
                            "id": "#main/AlignmentAnalysis/Transcript_Length"
                        }
                    ],
                    "out": [
                        "#main/AlignmentAnalysis/Seq_Metrics",
                        "#main/AlignmentAnalysis/Annotated_Bam_Files",
                        "#main/AlignmentAnalysis/Sorted_Valid_Reads_CSV",
                        "#main/AlignmentAnalysis/num_valid_ig_reads",
                        "#main/AlignmentAnalysis/num_valid_tcr_reads",
                        "#main/AlignmentAnalysis/validIgReads",
                        "#main/AlignmentAnalysis/validTcrReads",
                        "#main/AlignmentAnalysis/num_cell_estimate",
                        "#main/AlignmentAnalysis/num_bioproducts"
                    ],
                    "id": "#main/AlignmentAnalysis"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 32000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#AnnotateMolecules.cwl",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/AbSeq_UMI",
                            "id": "#main/AnnotateMolecules/AbSeq_UMI"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateMolecules/Run_Metadata"
                        },
                        {
                            "source": "#main/Internal_Settings/Use_DBEC",
                            "id": "#main/AnnotateMolecules/Use_DBEC"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Sorted_Valid_Reads_CSV",
                            "id": "#main/AnnotateMolecules/Valids"
                        }
                    ],
                    "out": [
                        "#main/AnnotateMolecules/Mol_Annot_List",
                        "#main/AnnotateMolecules/Gene_Status_List",
                        "#main/AnnotateMolecules/Max_Count",
                        "#main/AnnotateMolecules/Total_Molecules",
                        "#main/AnnotateMolecules/output"
                    ],
                    "scatter": [
                        "#main/AnnotateMolecules/Valids"
                    ],
                    "id": "#main/AnnotateMolecules"
                },
                {
                    "run": "#Assay_Settings.cwl",
                    "in": [
                        {
                            "valueFrom": "WTA",
                            "id": "#main/Assay_Settings/_Assay"
                        }
                    ],
                    "out": [
                        "#main/Assay_Settings/Assay"
                    ],
                    "id": "#main/Assay_Settings"
                },
                {
                    "label": "Bam Settings",
                    "run": "#BamSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Generate_Bam",
                            "id": "#main/Bam_Settings/_Generate_Bam"
                        }
                    ],
                    "out": [
                        "#main/Bam_Settings/Generate_Bam"
                    ],
                    "id": "#main/Bam_Settings"
                },
                {
                    "run": "#BundleLogs.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/CheckReference/output",
                                "#main/GetDataTable/output",
                                "#main/Metrics/output",
                                "#main/AddtoBam/output",
                                "#main/AnnotateMolecules/output",
                                "#main/MergeBAM/log",
                                "#main/Dense_to_Sparse_Datatable/output",
                                "#main/Dense_to_Sparse_Datatable_Unfiltered/output",
                                "#main/IndexBAM/log"
                            ],
                            "pickValue": "all_non_null",
                            "linkMerge": "merge_flattened",
                            "id": "#main/BundleLogs/log_files"
                        }
                    ],
                    "out": [
                        "#main/BundleLogs/logs_dir"
                    ],
                    "id": "#main/BundleLogs"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 10000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#CheckReference.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "id": "#main/CheckReference/AbSeq_Reference"
                        },
                        {
                            "source": "#main/Assay_Settings/Assay",
                            "id": "#main/CheckReference/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/CheckReference/Putative_Cell_Call"
                        },
                        {
                            "source": [
                                "#main/Transcriptome_Annotation",
                                "#main/Reference_Genome"
                            ],
                            "id": "#main/CheckReference/Reference"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/CheckReference/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Supplemental_Reference",
                            "id": "#main/CheckReference/Supplemental_Reference"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/CheckReference/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/CheckReference/Index",
                        "#main/CheckReference/Extra_Seqs",
                        "#main/CheckReference/Full_Genes",
                        "#main/CheckReference/output",
                        "#main/CheckReference/Transcript_Length",
                        "#main/CheckReference/GTF",
                        "#main/CheckReference/Target_Gene_Mapping"
                    ],
                    "id": "#main/CheckReference"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 16000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#DensetoSparse.cwl",
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_File/Cell_Order",
                            "id": "#main/Dense_to_Sparse_Datatable/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Dense_Data_Tables",
                            "id": "#main/Dense_to_Sparse_Datatable/Dense_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Gene_List",
                            "id": "#main/Dense_to_Sparse_Datatable/Gene_List"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Dense_to_Sparse_Datatable/Run_Metadata"
                        }
                    ],
                    "out": [
                        "#main/Dense_to_Sparse_Datatable/Data_Tables",
                        "#main/Dense_to_Sparse_Datatable/output"
                    ],
                    "scatter": [
                        "#main/Dense_to_Sparse_Datatable/Dense_Data_Table"
                    ],
                    "id": "#main/Dense_to_Sparse_Datatable"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 16000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#DensetoSparse.cwl",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Dense_Data_Tables_Unfiltered",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Dense_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Gene_List",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Gene_List"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Run_Metadata"
                        }
                    ],
                    "out": [
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/Data_Tables",
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/output"
                    ],
                    "scatter": [
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/Dense_Data_Table"
                    ],
                    "id": "#main/Dense_to_Sparse_Datatable_Unfiltered"
                },
                {
                    "run": "#DensetoSparseFile.cwl",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/Dense_to_Sparse_File/GDT_cell_order"
                        }
                    ],
                    "out": [
                        "#main/Dense_to_Sparse_File/Cell_Order"
                    ],
                    "id": "#main/Dense_to_Sparse_File"
                },
                {
                    "run": "#GetDataTable.cwl",
                    "in": [
                        {
                            "source": "#main/CheckReference/Full_Genes",
                            "id": "#main/GetDataTable/Full_Genes"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Gene_Status_List",
                            "id": "#main/GetDataTable/Gene_Status_List"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Max_Count",
                            "id": "#main/GetDataTable/Max_Count"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Mol_Annot_List",
                            "id": "#main/GetDataTable/Molecule_Annotation_List"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/GetDataTable/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/GetDataTable/Run_Metadata"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tag_Names",
                            "id": "#main/GetDataTable/Sample_Tag_Names"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Seq_Metrics",
                            "id": "#main/GetDataTable/Seq_Metrics"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Total_Molecules",
                            "id": "#main/GetDataTable/Total_Molecules"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_bioproducts",
                            "id": "#main/GetDataTable/num_bioproducts"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_cell_estimate",
                            "id": "#main/GetDataTable/num_cell_estimate"
                        }
                    ],
                    "out": [
                        "#main/GetDataTable/Tag_Calls",
                        "#main/GetDataTable/Molecular_Annotation",
                        "#main/GetDataTable/Corrected_Molecular_Annotation",
                        "#main/GetDataTable/Tag_Annotation",
                        "#main/GetDataTable/Annot_Files",
                        "#main/GetDataTable/Cell_Label_Graph_Data",
                        "#main/GetDataTable/Dense_Data_Tables",
                        "#main/GetDataTable/Dense_Data_Tables_Unfiltered",
                        "#main/GetDataTable/Expression_Data",
                        "#main/GetDataTable/Expression_Data_Unfiltered",
                        "#main/GetDataTable/Bioproduct_Stats",
                        "#main/GetDataTable/UMI_Adjusted_CellLabel_Stats",
                        "#main/GetDataTable/Putative_Cells_Origin",
                        "#main/GetDataTable/Protein_Aggregates_Experimental",
                        "#main/GetDataTable/SampleTag_out",
                        "#main/GetDataTable/SampleTag_zip",
                        "#main/GetDataTable/output",
                        "#main/GetDataTable/Cell_Order",
                        "#main/GetDataTable/Gene_List",
                        "#main/GetDataTable/Dim_Reduction_Coord",
                        "#main/GetDataTable/Visual_Metrics_MolsPerCell",
                        "#main/GetDataTable/Cell_Type_Predictions"
                    ],
                    "id": "#main/GetDataTable"
                },
                {
                    "run": "#IndexBAM.cwl",
                    "in": [
                        {
                            "source": "#main/MergeBAM/Bam",
                            "id": "#main/IndexBAM/BamFile"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/IndexBAM/Generate_Bam"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true )",
                    "out": [
                        "#main/IndexBAM/Index",
                        "#main/IndexBAM/log"
                    ],
                    "id": "#main/IndexBAM"
                },
                {
                    "label": "Internal Settings",
                    "run": "#InternalSettings.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_UMI",
                            "id": "#main/Internal_Settings/_AbSeq_UMI"
                        },
                        {
                            "source": "#main/Barcode_Num",
                            "id": "#main/Internal_Settings/_Barcode_Num"
                        },
                        {
                            "source": "#main/Extra_Seqs",
                            "id": "#main/Internal_Settings/_Extra_Seqs"
                        },
                        {
                            "source": "#main/Read_Filter_Off",
                            "id": "#main/Internal_Settings/_Read_Filter_Off"
                        },
                        {
                            "source": "#main/Seq_Run",
                            "id": "#main/Internal_Settings/_Seq_Run"
                        },
                        {
                            "source": "#main/Target_analysis",
                            "id": "#main/Internal_Settings/_Target_analysis"
                        },
                        {
                            "source": "#main/Use_DBEC",
                            "id": "#main/Internal_Settings/_Use_DBEC"
                        },
                        {
                            "source": "#main/VDJ_JGene_Evalue",
                            "id": "#main/Internal_Settings/_VDJ_JGene_Evalue"
                        },
                        {
                            "source": "#main/VDJ_VGene_Evalue",
                            "id": "#main/Internal_Settings/_VDJ_VGene_Evalue"
                        }
                    ],
                    "out": [
                        "#main/Internal_Settings/Read_Filter_Off",
                        "#main/Internal_Settings/Barcode_Num",
                        "#main/Internal_Settings/Seq_Run",
                        "#main/Internal_Settings/AbSeq_UMI",
                        "#main/Internal_Settings/Use_DBEC",
                        "#main/Internal_Settings/Extra_Seqs",
                        "#main/Internal_Settings/Target_analysis",
                        "#main/Internal_Settings/VDJ_VGene_Evalue",
                        "#main/Internal_Settings/VDJ_JGene_Evalue"
                    ],
                    "id": "#main/Internal_Settings"
                },
                {
                    "label": "Intronic Reads Settings",
                    "run": "#IntronicReadsSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Exclude_Intronic_Reads",
                            "id": "#main/Intronic_Reads_Settings/_Exclude_Intronic_Reads"
                        }
                    ],
                    "out": [
                        "#main/Intronic_Reads_Settings/Exclude_Intronic_Reads"
                    ],
                    "id": "#main/Intronic_Reads_Settings"
                },
                {
                    "run": "#MergeBAM.cwl",
                    "in": [
                        {
                            "source": "#main/AddtoBam/Annotated_Bam",
                            "id": "#main/MergeBAM/BamFiles"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/MergeBAM/Generate_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/MergeBAM/Run_Metadata"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true )",
                    "out": [
                        "#main/MergeBAM/Bam",
                        "#main/MergeBAM/log"
                    ],
                    "id": "#main/MergeBAM"
                },
                {
                    "run": {
                        "cwlVersion": "v1.2",
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": {
                                    "items": [
                                        "null",
                                        "File"
                                    ],
                                    "type": "array"
                                },
                                "id": "#main/MergeMultiplex/run/SampleTag_Files"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    {
                                        "type": "array",
                                        "items": "File"
                                    }
                                ],
                                "id": "#main/MergeMultiplex/run/Multiplex_out"
                            }
                        ],
                        "expression": "${\n  var fp_array = [];\n  for (var i = 0; i < inputs.SampleTag_Files.length; i++) {\n    var fp = inputs.SampleTag_Files[i];\n    if (fp != null) {\n      fp_array.push(fp);\n    }\n  }\n  return({\"Multiplex_out\": fp_array});\n}"
                    },
                    "in": [
                        {
                            "source": [
                                "#main/GetDataTable/SampleTag_out",
                                "#main/Metrics/Sample_Tag_Out"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/MergeMultiplex/SampleTag_Files"
                        }
                    ],
                    "out": [
                        "#main/MergeMultiplex/Multiplex_out"
                    ],
                    "id": "#main/MergeMultiplex"
                },
                {
                    "run": "#Metadata.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "id": "#main/Metadata_Settings/AbSeq_Reference"
                        },
                        {
                            "source": "#main/Assay_Settings/Assay",
                            "id": "#main/Metadata_Settings/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Basic_Algo_Only",
                            "id": "#main/Metadata_Settings/Basic_Algo_Only"
                        },
                        {
                            "source": "#main/QualCLAlign/Bead_Version",
                            "id": "#main/Metadata_Settings/Bead_Version"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                            "id": "#main/Metadata_Settings/Exact_Cell_Count"
                        },
                        {
                            "source": "#main/Intronic_Reads_Settings/Exclude_Intronic_Reads",
                            "id": "#main/Metadata_Settings/Exclude_Intronic_Reads"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/Metadata_Settings/Generate_Bam"
                        },
                        {
                            "source": "#main/QualCLAlign/Libraries",
                            "id": "#main/Metadata_Settings/Libraries"
                        },
                        {
                            "valueFrom": "BD Rhapsody WTA Analysis Pipeline",
                            "id": "#main/Metadata_Settings/Pipeline_Name"
                        },
                        {
                            "source": "#main/Version/version",
                            "id": "#main/Metadata_Settings/Pipeline_Version"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/Metadata_Settings/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/QualCLAlign/ReadsList",
                            "id": "#main/Metadata_Settings/Reads"
                        },
                        {
                            "source": [
                                "#main/Transcriptome_Annotation",
                                "#main/Reference_Genome"
                            ],
                            "id": "#main/Metadata_Settings/Reference"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Name",
                            "id": "#main/Metadata_Settings/Run_Name"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tag_Names",
                            "id": "#main/Metadata_Settings/Sample_Tag_Names"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/Metadata_Settings/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Start_Time/Start_Time",
                            "id": "#main/Metadata_Settings/Start_Time"
                        },
                        {
                            "source": "#main/Supplemental_Reference",
                            "id": "#main/Metadata_Settings/Supplemental_Reference"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/Metadata_Settings/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/Metadata_Settings/Run_Metadata"
                    ],
                    "id": "#main/Metadata_Settings"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#Metrics.cwl",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Annot_Files",
                            "id": "#main/Metrics/Annot_Files"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Label_Graph_Data",
                            "id": "#main/Metrics/Cell_Label_Graph_Data"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Type_Predictions",
                            "id": "#main/Metrics/Cell_Type_Predictions"
                        },
                        {
                            "source": "#main/Dense_to_Sparse_Datatable/Data_Tables",
                            "id": "#main/Metrics/Data_Tables"
                        },
                        {
                            "source": "#main/GetDataTable/Dim_Reduction_Coord",
                            "id": "#main/Metrics/Dim_Reduction_Coord"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Metrics/Run_Metadata"
                        },
                        {
                            "source": "#main/GetDataTable/SampleTag_zip",
                            "id": "#main/Metrics/Sample_Tag_Archives"
                        },
                        {
                            "source": "#main/GetDataTable/SampleTag_out",
                            "id": "#main/Metrics/Sample_Tag_Files"
                        },
                        {
                            "source": "#main/Internal_Settings/Seq_Run",
                            "id": "#main/Metrics/Seq_Run"
                        },
                        {
                            "source": "#main/GetDataTable/UMI_Adjusted_CellLabel_Stats",
                            "id": "#main/Metrics/UMI_Adjusted_Stats"
                        },
                        {
                            "source": "#main/GetDataTable/Visual_Metrics_MolsPerCell",
                            "id": "#main/Metrics/Visual_Metrics_MolsPerCell"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/vdjCellsDatatable",
                            "id": "#main/Metrics/vdjCellsDatatable"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/vdjMetricsJson",
                            "id": "#main/Metrics/vdjMetricsJson"
                        }
                    ],
                    "out": [
                        "#main/Metrics/Metrics_Summary",
                        "#main/Metrics/Metrics_Archive",
                        "#main/Metrics/Visual_Metrics_json",
                        "#main/Metrics/Visual_Metrics_html",
                        "#main/Metrics/output",
                        "#main/Metrics/Sample_Tag_Out"
                    ],
                    "id": "#main/Metrics"
                },
                {
                    "label": "Multiplexing Settings",
                    "run": "#MultiplexingSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Tag_Names",
                            "id": "#main/Multiplexing_Settings/_Sample_Tag_Names"
                        },
                        {
                            "source": "#main/Sample_Tags_Version",
                            "id": "#main/Multiplexing_Settings/_Sample_Tags_Version"
                        }
                    ],
                    "out": [
                        "#main/Multiplexing_Settings/Sample_Tag_Names",
                        "#main/Multiplexing_Settings/Sample_Tags_Version"
                    ],
                    "id": "#main/Multiplexing_Settings"
                },
                {
                    "label": "Name Settings",
                    "run": "#NameSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Run_Name",
                            "id": "#main/Name_Settings/_Run_Name"
                        }
                    ],
                    "out": [
                        "#main/Name_Settings/Run_Name"
                    ],
                    "id": "#main/Name_Settings"
                },
                {
                    "label": "Putative Cell Calling Settings",
                    "run": "#PutativeCellSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Basic_Algo_Only",
                            "id": "#main/Putative_Cell_Calling_Settings/_Basic_Algo_Only"
                        },
                        {
                            "source": "#main/Exact_Cell_Count",
                            "id": "#main/Putative_Cell_Calling_Settings/_Exact_Cell_Count"
                        },
                        {
                            "source": "#main/Putative_Cell_Call",
                            "id": "#main/Putative_Cell_Calling_Settings/_Putative_Cell_Call"
                        }
                    ],
                    "out": [
                        "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                        "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                        "#main/Putative_Cell_Calling_Settings/Basic_Algo_Only"
                    ],
                    "id": "#main/Putative_Cell_Calling_Settings"
                },
                {
                    "run": "#QualCLAlign.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay",
                            "id": "#main/QualCLAlign/Assay"
                        },
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/QualCLAlign/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/Index",
                            "id": "#main/QualCLAlign/Index"
                        },
                        {
                            "source": "#main/Maximum_Threads",
                            "id": "#main/QualCLAlign/Maximum_Threads"
                        },
                        {
                            "source": "#main/Reads",
                            "id": "#main/QualCLAlign/Reads"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Name",
                            "id": "#main/QualCLAlign/Run_Name"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/QualCLAlign/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/QualCLAlign/Bead_Version",
                        "#main/QualCLAlign/Libraries",
                        "#main/QualCLAlign/ReadsList",
                        "#main/QualCLAlign/R2BAMFiles",
                        "#main/QualCLAlign/QualCLAlignMetrics",
                        "#main/QualCLAlign/Logs",
                        "#main/QualCLAlign/Run_Base_Name"
                    ],
                    "id": "#main/QualCLAlign"
                },
                {
                    "run": {
                        "cwlVersion": "v1.2",
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "inputs": [],
                        "outputs": [
                            {
                                "type": "string",
                                "id": "#main/Start_Time/run/Start_Time"
                            }
                        ],
                        "expression": "${   \n  var today = new Date();\n  var date = today.toString()\n  return ({Start_Time: date});\n} "
                    },
                    "in": [],
                    "out": [
                        "#main/Start_Time/Start_Time"
                    ],
                    "id": "#main/Start_Time"
                },
                {
                    "run": "#UncompressDatatables.cwl",
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_Datatable/Data_Tables",
                            "id": "#main/Uncompress_Datatables/Compressed_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Expression_Data",
                            "id": "#main/Uncompress_Datatables/Compressed_Expression_Matrix"
                        }
                    ],
                    "out": [
                        "#main/Uncompress_Datatables/Uncompressed_Data_Tables",
                        "#main/Uncompress_Datatables/Uncompressed_Expression_Matrix"
                    ],
                    "id": "#main/Uncompress_Datatables"
                },
                {
                    "run": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl",
                    "in": [
                        {
                            "source": "#main/VDJ_Preprocess_Reads_IG/RSEC_Reads_Fastq",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/VDJ_Version"
                        },
                        {
                            "source": "#main/VDJ_Preprocess_Reads_IG/num_cores",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/num_cores"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Assemble_and_Annotate_Contigs_IG/igCalls"
                    ],
                    "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG"
                },
                {
                    "run": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl",
                    "in": [
                        {
                            "source": "#main/VDJ_Preprocess_Reads_TCR/RSEC_Reads_Fastq",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/VDJ_Version"
                        },
                        {
                            "source": "#main/VDJ_Preprocess_Reads_TCR/num_cores",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/num_cores"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/tcrCalls"
                    ],
                    "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR"
                },
                {
                    "run": "#VDJ_Compile_Results.cwl",
                    "in": [
                        {
                            "source": "#main/AlignmentAnalysis/Seq_Metrics",
                            "id": "#main/VDJ_Compile_Results/Seq_Metrics"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Compile_Results/VDJ_Version"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Type_Predictions",
                            "id": "#main/VDJ_Compile_Results/cellTypeMapping"
                        },
                        {
                            "valueFrom": "$([])",
                            "id": "#main/VDJ_Compile_Results/chainsToIgnore"
                        },
                        {
                            "source": "#main/Internal_Settings/VDJ_JGene_Evalue",
                            "id": "#main/VDJ_Compile_Results/evalueJgene"
                        },
                        {
                            "source": "#main/Internal_Settings/VDJ_VGene_Evalue",
                            "id": "#main/VDJ_Compile_Results/evalueVgene"
                        },
                        {
                            "source": "#main/VDJ_GatherIGCalls/gatheredCalls",
                            "id": "#main/VDJ_Compile_Results/igCalls"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/VDJ_Compile_Results/metadata"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/VDJ_Compile_Results/putativeCells"
                        },
                        {
                            "source": "#main/VDJ_GatherTCRCalls/gatheredCalls",
                            "id": "#main/VDJ_Compile_Results/tcrCalls"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Compile_Results/vdjCellsDatatable",
                        "#main/VDJ_Compile_Results/vdjCellsDatatableUncorrected",
                        "#main/VDJ_Compile_Results/vdjDominantContigsAIRR",
                        "#main/VDJ_Compile_Results/vdjUnfilteredContigsAIRR",
                        "#main/VDJ_Compile_Results/vdjMetricsJson",
                        "#main/VDJ_Compile_Results/vdjMetricsCsv",
                        "#main/VDJ_Compile_Results/vdjReadsPerCellByChainTypeFigure"
                    ],
                    "id": "#main/VDJ_Compile_Results"
                },
                {
                    "run": "#VDJ_GatherCalls.cwl",
                    "in": [
                        {
                            "source": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/igCalls",
                            "id": "#main/VDJ_GatherIGCalls/theCalls"
                        }
                    ],
                    "out": [
                        "#main/VDJ_GatherIGCalls/gatheredCalls"
                    ],
                    "id": "#main/VDJ_GatherIGCalls"
                },
                {
                    "run": "#VDJ_GatherCalls.cwl",
                    "in": [
                        {
                            "source": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/tcrCalls",
                            "id": "#main/VDJ_GatherTCRCalls/theCalls"
                        }
                    ],
                    "out": [
                        "#main/VDJ_GatherTCRCalls/gatheredCalls"
                    ],
                    "id": "#main/VDJ_GatherTCRCalls"
                },
                {
                    "run": "#VDJ_Preprocess_Reads.cwl",
                    "in": [
                        {
                            "source": "#main/AlignmentAnalysis/validIgReads",
                            "id": "#main/VDJ_Preprocess_Reads_IG/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_valid_ig_reads",
                            "id": "#main/VDJ_Preprocess_Reads_IG/num_valid_reads"
                        },
                        {
                            "valueFrom": "BCR",
                            "id": "#main/VDJ_Preprocess_Reads_IG/vdj_type"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Preprocess_Reads_IG/RSEC_Reads_Fastq",
                        "#main/VDJ_Preprocess_Reads_IG/num_splits",
                        "#main/VDJ_Preprocess_Reads_IG/num_cores"
                    ],
                    "id": "#main/VDJ_Preprocess_Reads_IG"
                },
                {
                    "run": "#VDJ_Preprocess_Reads.cwl",
                    "in": [
                        {
                            "source": "#main/AlignmentAnalysis/validTcrReads",
                            "id": "#main/VDJ_Preprocess_Reads_TCR/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_valid_tcr_reads",
                            "id": "#main/VDJ_Preprocess_Reads_TCR/num_valid_reads"
                        },
                        {
                            "valueFrom": "TCR",
                            "id": "#main/VDJ_Preprocess_Reads_TCR/vdj_type"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Preprocess_Reads_TCR/RSEC_Reads_Fastq",
                        "#main/VDJ_Preprocess_Reads_TCR/num_splits",
                        "#main/VDJ_Preprocess_Reads_TCR/num_cores"
                    ],
                    "id": "#main/VDJ_Preprocess_Reads_TCR"
                },
                {
                    "label": "VDJ Settings",
                    "run": "#VDJ_Settings.cwl",
                    "in": [
                        {
                            "source": "#main/VDJ_Version",
                            "id": "#main/VDJ_Settings/_VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Settings/VDJ_Version"
                    ],
                    "id": "#main/VDJ_Settings"
                },
                {
                    "run": "#Version.cwl",
                    "in": [],
                    "out": [
                        "#main/Version/version"
                    ],
                    "id": "#main/Version"
                }
            ],
            "id": "#main"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "class": "ResourceRequirement",
                    "coresMin": 4
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "samtools",
                "merge"
            ],
            "stdout": "samtools_merge.log",
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#MergeBAM.cwl/BamFiles"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#MergeBAM.cwl/Generate_Bam"
                },
                {
                    "type": "File",
                    "loadContents": true,
                    "id": "#MergeBAM.cwl/Run_Metadata"
                }
            ],
            "arguments": [
                {
                    "prefix": "-@",
                    "valueFrom": "$(runtime.cores)"
                },
                {
                    "position": 0,
                    "valueFrom": "${\n    var st_version = JSON.parse(inputs.Run_Metadata.contents).Sample_Tags_Version\n    var run_name = JSON.parse(inputs.Run_Metadata.contents).Run_Base_Name\n    if (st_version) {\n        return \"Combined_\" + run_name + \".BAM\"\n    } else {\n        return run_name + \".BAM\"\n    }\n}"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "${ return \"*\" + JSON.parse(inputs.Run_Metadata.contents).Run_Base_Name + \".BAM\" }"
                    },
                    "id": "#MergeBAM.cwl/Bam"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#MergeBAM.cwl/log"
                }
            ],
            "id": "#MergeBAM.cwl"
        },
        {
            "class": "CommandLineTool",
            "baseCommand": "echo",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#Metadata.cwl/AbSeq_Reference"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Basic_Algo_Only"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#Metadata.cwl/Bead_Version/Library",
                                    "type": "string"
                                },
                                {
                                    "name": "#Metadata.cwl/Bead_Version/bead_version",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "id": "#Metadata.cwl/Bead_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Exclude_Intronic_Reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Generate_Bam"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Label_Version"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Libraries"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Pipeline_Name"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Pipeline_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Putative_Cell_Call"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Read_Filter_Off"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#Metadata.cwl/Reads"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#Metadata.cwl/Reference"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Run_Name"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#Metadata.cwl/Sample_Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Start_Time"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#Metadata.cwl/Subsample_Reads"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#Metadata.cwl/Subsample_Sample_Tags"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Subsample_Seed"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#Metadata.cwl/Supplemental_Reference"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "outputBinding": {
                        "outputEval": "${  \n  var name = inputs.Run_Name;\n  if (name == null){\n    var libraries = inputs.Libraries;\n    name = libraries.split(',')[0];\n  }   \n  return(name)\n}   \n"
                    },
                    "id": "#Metadata.cwl/Run_Base_Name"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "run_metadata.json"
                    },
                    "id": "#Metadata.cwl/Run_Metadata"
                }
            ],
            "stdout": "run_metadata.json",
            "arguments": [
                {
                    "prefix": ""
                },
                {
                    "shellQuote": true,
                    "valueFrom": "${\n  var metadata = inputs;\n  var all_bv = {};\n  var customer_bv = \"Original (V1)\";\n  var detected_bv = \"V1\";\n  for (var i = 0; i < inputs.Bead_Version.length; i++) {\n      var BeadVer = inputs.Bead_Version[i];\n      var Library = BeadVer[\"Library\"];\n      var bead_version = BeadVer[\"bead_version\"];\n      all_bv[Library] = bead_version  \n      var short_bv =  bead_version.substring(0, 5);\n      if (short_bv == \"Enh\") {\n        customer_bv = \"Enhanced\";\n        detected_bv = \"Enh\";\n      }\n      else if (short_bv == \"EnhV2\") {\n        customer_bv = \"Enhanced V2\";\n        detected_bv = \"EnhV2\";\n      }\n  }\n  metadata[\"Bead_Version\"] = all_bv;\n  metadata[\"Bead_Version_Detected\"] = detected_bv;\n\n  var pipeline_name = inputs.Pipeline_Name;\n  var assay = inputs.Assay;\n  var version = inputs.Pipeline_Version;\n  var time = inputs.Start_Time;\n  var libraries = inputs.Libraries.split(\",\");\n  var i = 0;\n  var reference_list = []\n  if(inputs.Reference != null){\n      reference_list = reference_list.concat(inputs.Reference);\n  }\n  if(inputs.AbSeq_Reference != null){\n      reference_list = reference_list.concat(inputs.AbSeq_Reference);\n  }\n\n  var supplemental = \"\"\n  if(inputs.Supplemental_Reference != null){\n      supplemental = \"; Supplemental_Reference - \" + inputs.Supplemental_Reference[0][\"basename\"];\n  }\n  var references = [];\n  for (i = 0; i< reference_list.length; i++) {\n      if(reference_list[i] != null){\n          references.push(reference_list[i][\"basename\"]);\n      }\n  }\n  var parameters = [];\n  if(inputs.Sample_Tags_Version != null){\n      var tags = \"Sample Tag Version: \" + inputs.Sample_Tags_Version;\n  } else{ \n      var tags = \"Sample Tag Version: None\";\n  }\n  parameters.push(tags);\n\n  if(inputs.Sample_Tag_Names != null){\n      var tag_names = inputs.Sample_Tag_Names.join(\" ; \")\n      var tag_list = \"Sample Tag Names: \" + tag_names;\n  } else{\n      var tag_list = \"Sample Tag Names: None\";\n  }\n  parameters.push(tag_list);\n \n  if(inputs.VDJ_Version != null){\n      var vdj = \"VDJ Version: \" + inputs.VDJ_Version;\n  } else{ \n      var vdj = \"VDJ Version: None\";\n  }\n  parameters.push(vdj)\n\n  if(inputs.Subsample_Reads != null){\n      var subsample = \"Subsample Reads: \" + inputs.Subsample_Reads;\n  } else{ \n      var subsample = \"Subsample Reads: None\";\n  }   \n  parameters.push(subsample);\n\n  if(inputs.Putative_Cell_Call == 1){\n      var call = \"Putative Cell Calling Type: AbSeq\";\n  } else{ \n      var call = \"Putative Cell Calling Type: mRNA\";\n  }   \n  parameters.push(call)\n\n  if(inputs.Basic_Algo_Only){\n      var basic = \"Refined Putative Cell Calling: Off\";\n  } else{ \n      var basic = \"Refined Putative Cell Calling: On\";\n  }   \n  parameters.push(basic)\n\n  if(inputs.Exclude_Intronic_Reads){\n      var introns = \"Exclude Intronic Reads: On\";\n  } else{\n      var introns = \"Exclude Intronic Reads: Off\";\n  }\n  parameters.push(introns)\n\n  if(inputs.Generate_Bam){\n      var generateBam = \"Generate Bam: On\";\n  } else{\n      var generateBam = \"Generate Bam: Off\";\n  }\n  parameters.push(generateBam)\n\n  if(inputs.Exact_Cell_Count != null){\n      var cells = \"Exact Cell Count: \" + inputs.Exact_Cell_Count;\n  } else{\n      var cells = \"Exact Cell Count: None\";\n  }\n  parameters.push(cells)\n\n  var name = inputs.Run_Name;\n  if (name == null){\n    var libraries = inputs.Libraries.split(',');\n    name = libraries[0];\n  }        \n\n  var header = [\"####################\"];\n  header.push(\"## \" + pipeline_name + \" Version \" + version);\n  header.push(\"## Analysis Date - \" + time);\n  header.push(\"## Libraries - \" + libraries.join(' | ') + \" - Bead version detected: \" + customer_bv);\n  header.push(\"## References - \" + references.join(' | ') + supplemental);\n  header.push(\"## Parameters - \" + parameters.join(' | '));\n  header.push(\"####################\");\n  metadata[\"Output_Header\"] = header;\n  metadata[\"Run_Base_Name\"] = name;\n  var metadata_json = JSON.stringify(metadata);\n  return metadata_json;\n}\n"
                }
            ],
            "id": "#Metadata.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_metrics.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--annot-files"
                    },
                    "id": "#Metrics.cwl/Annot_Files"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--cell-label-graph-data-file"
                    },
                    "id": "#Metrics.cwl/Cell_Label_Graph_Data"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--cell-type-file"
                    },
                    "id": "#Metrics.cwl/Cell_Type_Predictions"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--data-tables"
                    },
                    "id": "#Metrics.cwl/Data_Tables"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--dim-reduction-coord-file"
                    },
                    "id": "#Metrics.cwl/Dim_Reduction_Coord"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#Metrics.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--sample-tag-archives"
                    },
                    "id": "#Metrics.cwl/Sample_Tag_Archives"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--sample-tag-files"
                    },
                    "id": "#Metrics.cwl/Sample_Tag_Files"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--seq-run"
                    },
                    "id": "#Metrics.cwl/Seq_Run"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--umi-adjusted-stats"
                    },
                    "id": "#Metrics.cwl/UMI_Adjusted_Stats"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--visual-metrics-mols-per-cell-file"
                    },
                    "id": "#Metrics.cwl/Visual_Metrics_MolsPerCell"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-per-cell-file"
                    },
                    "id": "#Metrics.cwl/vdjCellsDatatable"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-metrics-fp"
                    },
                    "id": "#Metrics.cwl/vdjMetricsJson"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "internal-metrics-archive.tar.gz"
                    },
                    "id": "#Metrics.cwl/Metrics_Archive"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_Metrics_Summary.csv"
                    },
                    "id": "#Metrics.cwl/Metrics_Summary"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "*.zip"
                    },
                    "id": "#Metrics.cwl/Sample_Tag_Out"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Pipeline_Report.html"
                    },
                    "id": "#Metrics.cwl/Visual_Metrics_html"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Visual_Metrics.json"
                    },
                    "id": "#Metrics.cwl/Visual_Metrics_json"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#Metrics.cwl/output"
                }
            ],
            "id": "#Metrics.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "default": "Targeted",
                    "id": "#MultiplexingSettings.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/_Sample_Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MultiplexingSettings.cwl/_Sample_Tags_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/Sample_Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MultiplexingSettings.cwl/Sample_Tags_Version"
                }
            ],
            "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version.indexOf('flex') >= 0)\n  {\n    enumifiedSampleTagsVersion = 'flex';\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  var listTagNames = inputs._Sample_Tag_Names\n  var newTagNames = []\n  for (var num in listTagNames) {\n    var tag = listTagNames[num].replace(/[^A-Za-z0-9-+]/g,\"_\");\n    newTagNames.push(tag);    \n  }  \n  return ({\n  Sample_Tag_Names: newTagNames,\n  Sample_Tags_Version: enumifiedSampleTagsVersion\n  });\n}",
            "id": "#MultiplexingSettings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#NameSettings.cwl/_Run_Name"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#NameSettings.cwl/Run_Name"
                }
            ],
            "expression": "${ var name = inputs._Run_Name;\n   if (name != null) {\n     name = name.trim().replace(/[\\W_]+/g,\"-\");\n     name = name.replace(/^-+/, '_');\n   }\n   if (name == '') {\n     name = null;\n   }\n   return({'Run_Name' : name });\n }  ",
            "id": "#NameSettings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Basic_Algo_Only"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Putative_Cell_Call"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#PutativeCellSettings.cwl/Basic_Algo_Only"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/Putative_Cell_Call"
                }
            ],
            "expression": "${\n  // the basic algorithm flag defaults to false\n  var basicAlgOnlyFlag = false;\n  // the user can set the basic algorithm flag\n  if (inputs._Basic_Algo_Only) {\n    basicAlgOnlyFlag = inputs._Basic_Algo_Only;\n  }\n  // convert the Putative_Cell_Call from a string to an integer\n  var putativeCellCallInt = 0;\n  if (inputs._Putative_Cell_Call) {\n    if (inputs._Putative_Cell_Call === \"mRNA\") {\n      putativeCellCallInt = 0;\n    }\n    else if (inputs._Putative_Cell_Call == \"AbSeq_Experimental\" || inputs._Putative_Cell_Call == \"AbSeq\") {\n      putativeCellCallInt = 1;\n      // for protein-only cell calling, we only have the basic algorithm\n      basicAlgOnlyFlag = true;\n    }\n    else if (inputs._Putative_Cell_Call == \"mRNA_and_AbSeq\") {\n      putativeCellCallInt = 2;\n    }\n  }\n  // check the exact cell count\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      throw(\"Illogical value for exact cell count: \" + inputs._Exact_Cell_Count);\n    }\n  }\n  return ({\n    Putative_Cell_Call: putativeCellCallInt,\n    Exact_Cell_Count: inputs._Exact_Cell_Count,\n    Basic_Algo_Only: basicAlgOnlyFlag,\n  });\n}",
            "id": "#PutativeCellSettings.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "doc": "CheckFastqs does several quality control routines including: ensuring that read pair file names are formatted correctly and contain a read pair mate;  QualCLAlign stage of the Rain pipeline overlaps read pairs and then performs a series of filters and mappings to reduce valid reads into a single FastQ file to be fed into the aligner. The R2 reads are annotated with cell index and UMI information derived from the R1 read.\n",
            "hints": [
                {
                    "class": "ResourceRequirement",
                    "coresMin": 8,
                    "ramMin": 48000
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "CheckFastqsAndQualCLAlign.sh"
            ],
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--alignment-compression-threads"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#QualCLAlign.cwl/Alignment_Compression_threads"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "prefix": "--assay"
                    },
                    "id": "#QualCLAlign.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--bgzf-threads"
                    },
                    "id": "#QualCLAlign.cwl/BGZF_Threads"
                },
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualCLAlign.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--index"
                    },
                    "type": "File",
                    "id": "#QualCLAlign.cwl/Index"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#QualCLAlign.cwl/Maximum_Threads"
                },
                {
                    "inputBinding": {
                        "prefix": "--reader-annotation-threads"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#QualCLAlign.cwl/Reader_Annotation_Threads"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "prefix": "--reads",
                        "itemSeparator": ","
                    },
                    "id": "#QualCLAlign.cwl/Reads"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--run-name"
                    },
                    "id": "#QualCLAlign.cwl/Run_Name"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-version"
                    },
                    "id": "#QualCLAlign.cwl/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#QualCLAlign.cwl/Bead_Version/Library",
                                    "type": "string"
                                },
                                {
                                    "name": "#QualCLAlign.cwl/Bead_Version/bead_version",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "outputBinding": {
                        "glob": "bead_version.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).BeadVersion)\n"
                    },
                    "id": "#QualCLAlign.cwl/Bead_Version"
                },
                {
                    "outputBinding": {
                        "glob": "fastq_read_pairs.json"
                    },
                    "type": "File",
                    "id": "#QualCLAlign.cwl/Fastq_read_pairs"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "outputBinding": {
                        "glob": "bead_version.json",
                        "loadContents": true,
                        "outputEval": "${\n  var obj = JSON.parse(self[0].contents);\n  var libraries = [];\n  var beadLibs = obj.BeadVersion\n  for (var i in beadLibs){\n    if (libraries.indexOf(beadLibs[i][\"Library\"]) == -1){ \n      libraries.push(beadLibs[i][\"Library\"]);\n    }\n  }\n  libraries.sort();\n  return(libraries.toString())\n}\n"
                    },
                    "id": "#QualCLAlign.cwl/Libraries"
                },
                {
                    "outputBinding": {
                        "glob": "*logs.tar.gz"
                    },
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#QualCLAlign.cwl/Logs"
                },
                {
                    "outputBinding": {
                        "glob": "ReadQualityMetrics.json"
                    },
                    "type": "File",
                    "id": "#QualCLAlign.cwl/QualCLAlignMetrics"
                },
                {
                    "outputBinding": {
                        "glob": "*_Aligned.out.bam"
                    },
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#QualCLAlign.cwl/R2BAMFiles"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "outputBinding": {
                        "outputEval": "${  \n  var reads = []; \n  var files = inputs.Reads\n  for (var i in files){\n      reads.push(files[i][\"basename\"]);\n  }\n  reads.sort();\n  return(reads)\n}\n"
                    },
                    "id": "#QualCLAlign.cwl/ReadsList"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "outputBinding": {
                        "glob": "run_base_name.json",
                        "loadContents": true,
                        "outputEval": "${\n  return (JSON.parse(self[0].contents));\n} \n"
                    },
                    "id": "#QualCLAlign.cwl/Run_Base_Name"
                }
            ],
            "id": "#QualCLAlign.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#UncompressDatatables.cwl/Compressed_Data_Table"
                },
                {
                    "type": "File",
                    "id": "#UncompressDatatables.cwl/Compressed_Expression_Matrix"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompressed_File",
                    "id": "#UncompressDatatables.cwl/Uncompressed_Data_Tables"
                },
                {
                    "type": "File",
                    "outputSource": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompressed_File",
                    "id": "#UncompressDatatables.cwl/Uncompressed_Expression_Matrix"
                }
            ],
            "steps": [
                {
                    "run": {
                        "hints": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:2.0b2",
                                "class": "DockerRequirement"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "cwlVersion": "v1.2",
                        "class": "CommandLineTool",
                        "baseCommand": [
                            "gunzip"
                        ],
                        "arguments": [
                            {
                                "position": 0,
                                "valueFrom": "-c"
                            }
                        ],
                        "inputs": [
                            {
                                "type": "File",
                                "inputBinding": {
                                    "position": 1
                                },
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/run/Uncompress_Datatable_Inner/Compressed_File"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "type": "File",
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/run/Uncompress_Datatable_Inner/Uncompressed_File"
                            }
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Datatable/run/Uncompress_Datatable_Inner"
                    },
                    "in": [
                        {
                            "source": "#UncompressDatatables.cwl/Compressed_Data_Table",
                            "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Compressed_File"
                        }
                    ],
                    "scatter": [
                        "#UncompressDatatables.cwl/Uncompress_Datatable/Compressed_File"
                    ],
                    "out": [
                        "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompressed_File"
                    ],
                    "id": "#UncompressDatatables.cwl/Uncompress_Datatable"
                },
                {
                    "run": {
                        "hints": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:2.0b2",
                                "class": "DockerRequirement"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "cwlVersion": "v1.2",
                        "class": "CommandLineTool",
                        "baseCommand": [
                            "gunzip"
                        ],
                        "arguments": [
                            {
                                "position": 0,
                                "valueFrom": "-c"
                            }
                        ],
                        "inputs": [
                            {
                                "type": "File",
                                "inputBinding": {
                                    "position": 1
                                },
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/run/Uncompress_Expression_Matrix_Inner/Compressed_File"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "type": "File",
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/run/Uncompress_Expression_Matrix_Inner/Uncompressed_File"
                            }
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/run/Uncompress_Expression_Matrix_Inner"
                    },
                    "in": [
                        {
                            "source": "#UncompressDatatables.cwl/Compressed_Expression_Matrix",
                            "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Compressed_File"
                        }
                    ],
                    "out": [
                        "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompressed_File"
                    ],
                    "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix"
                }
            ],
            "id": "#UncompressDatatables.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "class": "ResourceRequirement",
                    "coresMin": 1,
                    "ramMin": 3200
                }
            ],
            "baseCommand": [
                "AssembleAndAnnotate.sh"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/RSEC_Reads_Fastq"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/Read_Limit"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 3
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_pruned.csv.gz"
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/PyirCall"
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        {
                            "type": "array",
                            "items": [
                                "null",
                                "File"
                            ]
                        }
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/RSEC_Reads_Fastq"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/num_cores"
                }
            ],
            "outputs": [
                {
                    "type": [
                        {
                            "type": "array",
                            "items": [
                                "null",
                                "File"
                            ]
                        }
                    ],
                    "outputSource": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/PyirCall",
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/igCalls"
                }
            ],
            "steps": [
                {
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG",
                    "hints": [
                        {
                            "class": "ResourceRequirement",
                            "coresMin": "$(inputs.num_cores)"
                        }
                    ],
                    "run": "#VDJ_Assemble_and_Annotate_Contigs.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/RSEC_Reads_Fastq",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                        },
                        {
                            "valueFrom": "75000",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/Read_Limit"
                        },
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Version",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/PyirCall"
                    ],
                    "scatter": [
                        "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                    ]
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        {
                            "type": "array",
                            "items": [
                                "null",
                                "File"
                            ]
                        }
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/RSEC_Reads_Fastq"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/num_cores"
                }
            ],
            "outputs": [
                {
                    "type": [
                        {
                            "type": "array",
                            "items": [
                                "null",
                                "File"
                            ]
                        }
                    ],
                    "outputSource": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/PyirCall",
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/tcrCalls"
                }
            ],
            "steps": [
                {
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR",
                    "hints": [
                        {
                            "class": "ResourceRequirement",
                            "coresMin": "$(inputs.num_cores)"
                        }
                    ],
                    "run": "#VDJ_Assemble_and_Annotate_Contigs.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/RSEC_Reads_Fastq",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                        },
                        {
                            "valueFrom": "75000",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/Read_Limit"
                        },
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Version",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/PyirCall"
                    ],
                    "scatter": [
                        "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                    ]
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "baseCommand": [
                "mist_vdj_compile_results.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--seq-metrics",
                        "position": 10
                    },
                    "id": "#VDJ_Compile_Results.cwl/Seq_Metrics"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-version",
                        "position": 2
                    },
                    "id": "#VDJ_Compile_Results.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--cell-type-mapping-fp"
                    },
                    "id": "#VDJ_Compile_Results.cwl/cellTypeMapping"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--ignore",
                        "itemSeparator": ","
                    },
                    "id": "#VDJ_Compile_Results.cwl/chainsToIgnore"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 8,
                        "prefix": "--e-value-for-j"
                    },
                    "id": "#VDJ_Compile_Results.cwl/evalueJgene"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 7,
                        "prefix": "--e-value-for-v"
                    },
                    "id": "#VDJ_Compile_Results.cwl/evalueVgene"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 5
                    },
                    "id": "#VDJ_Compile_Results.cwl/igCalls"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 9,
                        "prefix": "--metadata-fp"
                    },
                    "id": "#VDJ_Compile_Results.cwl/metadata"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--putative-cells-json-fp",
                        "position": 3
                    },
                    "id": "#VDJ_Compile_Results.cwl/putativeCells"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 6
                    },
                    "id": "#VDJ_Compile_Results.cwl/tcrCalls"
                }
            ],
            "outputs": [
                {
                    "doc": "VDJ data per cell, with distribution based error correction",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_perCell.csv"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjCellsDatatable"
                },
                {
                    "doc": "VDJ data per cell, including non-putative cells, no error correction applied",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_perCell_uncorrected.csv.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjCellsDatatableUncorrected"
                },
                {
                    "doc": "AIRR compatible output that only reports the Dominant contigs, counts are DBEC corrected",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_Dominant_Contigs_AIRR.tsv"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjDominantContigsAIRR"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_metrics.csv"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjMetricsCsv"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_metrics.json"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjMetricsJson"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*_DBEC_cutoff.png"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjReadsPerCellByChainTypeFigure"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "AIRR compatible output that reports all the congits, counts are not DBEC corrected",
                    "outputBinding": {
                        "glob": "*_VDJ_Unfiltered_Contigs_AIRR.tsv.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjUnfilteredContigsAIRR"
                }
            ],
            "id": "#VDJ_Compile_Results.cwl"
        },
        {
            "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n",
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_GatherCalls.cwl/VDJ_Version"
                },
                {
                    "type": [
                        {
                            "type": "array",
                            "items": [
                                "null",
                                "File"
                            ]
                        }
                    ],
                    "id": "#VDJ_GatherCalls.cwl/theCalls"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls",
                    "id": "#VDJ_GatherCalls.cwl/gatheredCalls"
                }
            ],
            "steps": [
                {
                    "in": [
                        {
                            "source": "#VDJ_GatherCalls.cwl/theCalls",
                            "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/theCalls"
                        }
                    ],
                    "out": [
                        "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls"
                    ],
                    "run": {
                        "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/run/gather_PyIR",
                        "cwlVersion": "v1.2",
                        "class": "CommandLineTool",
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:2.0b2",
                                "class": "DockerRequirement"
                            },
                            {
                                "class": "InlineJavascriptRequirement"
                            },
                            {
                                "class": "ShellCommandRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": [
                                    {
                                        "type": "array",
                                        "items": [
                                            "null",
                                            "File"
                                        ]
                                    }
                                ],
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/run/gather_PyIR/theCalls"
                            }
                        ],
                        "arguments": [
                            {
                                "shellQuote": false,
                                "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from PyIR detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PrunePyIR\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PrunePyIR\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "outputBinding": {
                                    "glob": "*_constant_region_called_pruned.csv.gz",
                                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from PyIR detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                                },
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/run/gather_PyIR/gatheredCalls"
                            }
                        ]
                    },
                    "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls"
                }
            ],
            "id": "#VDJ_GatherCalls.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/Valid_Reads_Fastq"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/num_valid_reads"
                },
                {
                    "type": "string",
                    "id": "#VDJ_Preprocess_Reads.cwl/vdj_type"
                }
            ],
            "outputs": [
                {
                    "type": [
                        {
                            "type": "array",
                            "items": [
                                "null",
                                "File"
                            ]
                        }
                    ],
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/RSEC_Reads_Fastq",
                    "id": "#VDJ_Preprocess_Reads.cwl/RSEC_Reads_Fastq"
                },
                {
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_cores",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/num_cores"
                },
                {
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/num_splits"
                }
            ],
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "EnvVarRequirement",
                    "envDef": [
                        {
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                            "envValue": "8"
                        }
                    ]
                }
            ],
            "steps": [
                {
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads",
                    "requirements": [
                        {
                            "class": "ResourceRequirement",
                            "ramMin": "${ var est_ram = 0.0006 * parseInt(inputs.num_valid_reads) + 2000; var buffer = 1.25; est_ram *= buffer; if (est_ram < 2000) return 2000; if (est_ram > 370000) return 370000; return parseInt(est_ram); }",
                            "coresMin": 8
                        }
                    ],
                    "run": "#VDJ_RSEC_Reads.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/Valid_Reads"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/num_splits"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/num_valid_reads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/num_valid_reads"
                        }
                    ],
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/RSEC_Reads_Fastq"
                    ]
                },
                {
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads",
                    "hints": [
                        {
                            "class": "ResourceRequirement",
                            "coresMin": 8
                        }
                    ],
                    "run": "#VDJ_Trim_Reads.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/Valid_Reads_Fastq",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads_Fastq"
                        }
                    ],
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads",
                        "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Trim_Report"
                    ]
                },
                {
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits",
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/num_valid_reads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_valid_reads"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/vdj_type",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/vdj_type"
                        }
                    ],
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                        "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_cores"
                    ],
                    "run": {
                        "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits",
                        "cwlVersion": "v1.2",
                        "class": "ExpressionTool",
                        "inputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/num_valid_reads"
                            },
                            {
                                "type": "string",
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/vdj_type"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/num_cores"
                            },
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/num_splits"
                            }
                        ],
                        "expression": "${\n  var ram_per_instance = 192 * 1024;\n  var num_cores = 96;\n  if (inputs.vdj_type == \"BCR\") {\n    ram_per_instance = 144 * 1024;\n    num_cores = 72;\n  }\n  var ram_per_split = 3200;\n  var num_splits_per_instance = parseInt(ram_per_instance / ram_per_split);\n  var num_splits = num_splits_per_instance;\n\n  var num_reads = parseInt(inputs.num_valid_reads);\n  if (num_reads != null) {\n    if (num_reads > 100000000)\n      num_splits = num_splits_per_instance * 2;\n      num_cores = num_cores * 2;\n  }\n\n  return ({\"num_splits\": num_splits, \"num_cores\": num_cores});\n}"
                    }
                }
            ],
            "id": "#VDJ_Preprocess_Reads.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": "mist_vdj_rsec_reads.py",
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_RSEC_Reads.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-valid-reads",
                        "itemSeparator": ","
                    },
                    "id": "#VDJ_RSEC_Reads.cwl/Valid_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--num-splits"
                    },
                    "id": "#VDJ_RSEC_Reads.cwl/num_splits"
                }
            ],
            "outputs": [
                {
                    "type": [
                        {
                            "type": "array",
                            "items": [
                                "null",
                                "File"
                            ]
                        }
                    ],
                    "outputBinding": {
                        "glob": "*RSEC_Reads_Fastq_*.tar.gz"
                    },
                    "id": "#VDJ_RSEC_Reads.cwl/RSEC_Reads_Fastq"
                }
            ],
            "id": "#VDJ_RSEC_Reads.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#VDJ_Settings.cwl/_VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_VGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_Version"
                }
            ],
            "expression": "${\n  var vdjVersion = null;\n  if (!inputs._VDJ_Version) {\n    vdjVersion = null;}\n  else {\n    var _VDJ_Version = inputs._VDJ_Version.toLowerCase();\n    if (_VDJ_Version === \"human\" || _VDJ_Version === \"hs\" || _VDJ_Version === \"human vdj - bcr and tcr\") {\n      vdjVersion = \"human\";\n    } else if (_VDJ_Version === \"humanbcr\" || _VDJ_Version === \"human vdj - bcr only\") {\n      vdjVersion = \"humanBCR\";\n    } else if (_VDJ_Version === \"humantcr\" || _VDJ_Version === \"human vdj - tcr only\") {\n      vdjVersion = \"humanTCR\";\n    } else if (_VDJ_Version === \"mouse\" || _VDJ_Version === \"mm\" || _VDJ_Version === \"mouse vdj - bcr and tcr\") {\n      vdjVersion = \"mouse\";\n    } else if (_VDJ_Version === \"mousebcr\" || _VDJ_Version === \"mouse vdj - bcr only\") {\n      vdjVersion = \"mouseBCR\";\n    } else if (_VDJ_Version === \"mousetcr\" || _VDJ_Version === \"mouse vdj - tcr only\") {\n      vdjVersion = \"mouseTCR\";\n    } else {\n      vdjVersion = inputs._VDJ_Version;\n    }\n  }\n\n  return ({\n  VDJ_Version: vdjVersion,\n  })\n}",
            "id": "#VDJ_Settings.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": "VDJ_Trim_Reads.sh",
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Trim_Reads.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#VDJ_Trim_Reads.cwl/Valid_Reads_Fastq"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "cutadapt.log"
                    },
                    "id": "#VDJ_Trim_Reads.cwl/Trim_Report"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "*vdjtxt.gz"
                    },
                    "id": "#VDJ_Trim_Reads.cwl/Valid_Reads"
                }
            ],
            "id": "#VDJ_Trim_Reads.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:2.0b2",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_print_version.py"
            ],
            "stdout": "output.txt",
            "inputs": [],
            "outputs": [
                {
                    "type": "string",
                    "outputBinding": {
                        "glob": "output.txt",
                        "loadContents": true,
                        "outputEval": "$(self[0].contents)"
                    },
                    "id": "#Version.cwl/version"
                }
            ],
            "id": "#Version.cwl"
        }
    ],
    "cwlVersion": "v1.2",
    "$namespaces": {
        "sbg": "https://sevenbridges.com#",
        "arv": "http://arvados.org/cwl#"
    }
}