#!/usr/bin/env cwl-runner
{
    "$graph": [
        {
            "requirements": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "EnvVarRequirement",
                    "envDef": [
                        {
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                            "envValue": "$(String(runtime.cores))"
                        }
                    ]
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_align_R2.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "id": "#AlignR2.cwl/Extra_Seqs"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--index"
                    },
                    "id": "#AlignR2.cwl/Index"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--r2-fastqs"
                    },
                    "id": "#AlignR2.cwl/R2"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#AlignR2.cwl/Run_Metadata"
                }
            ],
            "id": "#AlignR2.cwl",
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*zip"
                    },
                    "id": "#AlignR2.cwl/Alignments"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AlignR2.cwl/output"
                }
            ]
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                },
                {
                    "ramMin": 2000,
                    "class": "ResourceRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_annotate_R1.py"
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
                    "inputBinding": {
                        "prefix": "--filter-metrics",
                        "itemSeparator": ","
                    },
                    "id": "#AnnotateR1.cwl/Filter_Metrics"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--R1"
                    },
                    "id": "#AnnotateR1.cwl/R1"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#AnnotateR1.cwl/Run_Metadata"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_Annotation_R1.csv.gz"
                    },
                    "id": "#AnnotateR1.cwl/Annotation_R1"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_R1_error_count_table.npy"
                    },
                    "id": "#AnnotateR1.cwl/R1_error_count_table"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_R1_read_count_breakdown.json"
                    },
                    "id": "#AnnotateR1.cwl/R1_read_count_breakdown"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AnnotateR1.cwl/output"
                }
            ],
            "id": "#AnnotateR1.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_annotate_R2.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "id": "#AnnotateR2.cwl/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "id": "#AnnotateR2.cwl/GTF_Annotation"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--R2-zip"
                    },
                    "id": "#AnnotateR2.cwl/R2_zip"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#AnnotateR2.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--transcript-length"
                    },
                    "id": "#AnnotateR2.cwl/Transcript_Length"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*Annotation_R2.csv.gz"
                    },
                    "id": "#AnnotateR2.cwl/Annot_R2"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*-annot.gtf"
                    },
                    "id": "#AnnotateR2.cwl/GTF"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*mapping_R2.BAM"
                    },
                    "id": "#AnnotateR2.cwl/R2_Bam"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_picard_quality_metrics.csv.gz"
                    },
                    "id": "#AnnotateR2.cwl/R2_Quality_Metrics"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AnnotateR2.cwl/output"
                }
            ],
            "id": "#AnnotateR2.cwl"
        },
        {
            "requirements": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entryname": "manifest.json",
                            "entry": "${\n    function getPaths(inputs, attribute) {\n      var fp_arr = []\n      for (var i = 0; i < inputs[attribute].length; i++)\n      {\n          fp_arr.push(inputs[attribute][i].path);\n      }\n      return fp_arr;\n    }\n    var paths = {}\n    paths['annotR1'] = getPaths(inputs, 'R1_Annotation')\n    paths['R1_error_count_table'] = getPaths(inputs, 'R1_error_count_table')\n    paths['R1_read_count_breakdown'] = getPaths(inputs, 'R1_read_count_breakdown')\n    paths['annotR2'] = getPaths(inputs, 'R2_Annotation')\n    paths['r2_quality_metrics_fps'] = getPaths(inputs, 'R2_Quality_Metrics')\n    if(inputs.Filter_Metrics[0] != null){\n        paths['filtering_stat_files'] = getPaths(inputs, 'Filter_Metrics')\n    }\n    var paths_json = JSON.stringify(paths);\n    return paths_json;\n}",
                            "writable": false
                        }
                    ]
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "EnvVarRequirement",
                    "envDef": [
                        {
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                            "envValue": "4"
                        }
                    ]
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_annotate_reads.py"
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
                    "id": "#AnnotateReads.cwl/AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--extra-seqs"
                    },
                    "id": "#AnnotateReads.cwl/Extra_Seqs"
                },
                {
                    "type": {
                        "items": [
                            "null",
                            "File"
                        ],
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/Filter_Metrics"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "id": "#AnnotateReads.cwl/Putative_Cell_Call"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#AnnotateReads.cwl/R1_Annotation"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#AnnotateReads.cwl/R1_error_count_table"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#AnnotateReads.cwl/R1_read_count_breakdown"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#AnnotateReads.cwl/R2_Annotation"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#AnnotateReads.cwl/R2_Quality_Metrics"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#AnnotateReads.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--target-gene-mapping"
                    },
                    "id": "#AnnotateReads.cwl/Target_Gene_Mapping"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Annotation_Read.csv.gz"
                    },
                    "id": "#AnnotateReads.cwl/Annotation_Read"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*read1_error_rate_archive*"
                    },
                    "id": "#AnnotateReads.cwl/Read1_error_rate"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_SeqMetrics.csv.gz"
                    },
                    "id": "#AnnotateReads.cwl/Seq_Metrics"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*Sorted_Valid_Reads.csv.*"
                    },
                    "id": "#AnnotateReads.cwl/Valid_Reads"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).BCR); }"
                    },
                    "id": "#AnnotateReads.cwl/num_valid_ig_reads"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).TCR); }"
                    },
                    "id": "#AnnotateReads.cwl/num_valid_tcr_reads"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AnnotateReads.cwl/output"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_IG_Valid_Reads.fastq.gz"
                    },
                    "id": "#AnnotateReads.cwl/validIgReads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_TCR_Valid_Reads.fastq.gz"
                    },
                    "id": "#AnnotateReads.cwl/validTcrReads"
                }
            ],
            "id": "#AnnotateReads.cwl"
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
            "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    log_file.basename = uuid() + \"-\" + log_file.basename;\n    listing.push(log_file);\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
            "id": "#BundleLogs.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_cell_classifier.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0
                    },
                    "id": "#Cell_Classifier.cwl/molsPerCellMatrix"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*cell_type_experimental.csv"
                    },
                    "id": "#Cell_Classifier.cwl/cellTypePredictions"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#Cell_Classifier.cwl/log"
                }
            ],
            "id": "#Cell_Classifier.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "CheckFastqs does several quality control routines including: (1) ensuring that read pair file names are formatted correctly and contain a read pair mate; (2) disambiguating the \"Subsample Reads\" input and; (3) if not provided, generating a subsampling seed that the downstream instances can use.\n",
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_check_fastqs.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--min-split-size"
                    },
                    "doc": "The minimum size (megabytes) of a file that should get split into chunks of a size designated in NumRecordsPerSplit\n",
                    "id": "#CheckFastqs.cwl/MinChunkSize"
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
                    "id": "#CheckFastqs.cwl/Reads"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "prefix": "--subsample"
                    },
                    "id": "#CheckFastqs.cwl/Subsample"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--subsample-seed"
                    },
                    "id": "#CheckFastqs.cwl/Subsample_Seed"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--subsample-seed"
                    },
                    "id": "#CheckFastqs.cwl/UserInputSubsampleSeed"
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
                                    "name": "#CheckFastqs.cwl/Bead_Version/Library",
                                    "type": "string"
                                },
                                {
                                    "name": "#CheckFastqs.cwl/Bead_Version/bead_version",
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
                    "id": "#CheckFastqs.cwl/Bead_Version"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/filename",
                                    "type": "string"
                                },
                                {
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/readFlag",
                                    "type": "string"
                                },
                                {
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/readPairId",
                                    "type": "string"
                                },
                                {
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/library",
                                    "type": "string"
                                },
                                {
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/beadVersion",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "outputBinding": {
                        "glob": "fastq_read_pairs.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).fastq_read_pairs)\n"
                    },
                    "id": "#CheckFastqs.cwl/FastqReadPairs"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "outputBinding": {
                        "glob": "files_to_skip_split_and_subsample.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).files_to_skip_split_and_subsample)\n"
                    },
                    "id": "#CheckFastqs.cwl/FilesToSkipSplitAndSubsample"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "outputBinding": {
                        "glob": "fastq_read_pairs.json",
                        "loadContents": true,
                        "outputEval": "${\n  var obj = JSON.parse(self[0].contents);\n  var libraries = [];\n  var pairs = obj.fastq_read_pairs\n  for (var i in pairs){\n    if (pairs[i][\"readFlag\"] == \"R1\"){\n      if (libraries.indexOf(pairs[i][\"library\"]) == -1){ \n        libraries.push(pairs[i][\"library\"]);\n      }\n    }\n  }\n  libraries.sort();\n  return(libraries.toString())\n}\n"
                    },
                    "id": "#CheckFastqs.cwl/Libraries"
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
                    "id": "#CheckFastqs.cwl/ReadsList"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "subsampling_info.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).subsampling_seed)\n"
                    },
                    "id": "#CheckFastqs.cwl/SubsampleSeed"
                },
                {
                    "type": "float",
                    "outputBinding": {
                        "glob": "subsampling_info.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).subsampling_ratio)\n"
                    },
                    "id": "#CheckFastqs.cwl/SubsamplingRatio"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#CheckFastqs.cwl/log"
                }
            ],
            "id": "#CheckFastqs.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#CheckReference.cwl/Run_Metadata"
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
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
                    "id": "#DensetoSparseFile.cwl/Cell_Order",
                    "outputBinding": {
                        "glob": "cell_order.json"
                    }
                }
            ],
            "id": "#DensetoSparseFile.cwl"
        },
        {
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": "${return Math.min(Math.max(parseInt(inputs.Total_Molecules.reduce(function(a, b) { return a + b; }, 0) / 4000), 32000), 768000);}"
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L"
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
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--seq-metrics"
                    },
                    "id": "#GetDataTable.cwl/Seq_Metrics"
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
                    "id": "#GetDataTable.cwl/Tag_Names"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "id": "#GetDataTable.cwl/Total_Molecules"
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
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*.png"
                    },
                    "id": "#GetDataTable.cwl/Cell_Label_Filter"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "cell_order.json"
                    },
                    "id": "#GetDataTable.cwl/Cell_Order"
                },
                {
                    "type": "File",
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
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Annotations/*_Annotation_Molecule_Trueno.csv"
                    },
                    "id": "#GetDataTable.cwl/Tag_Annotation"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Trueno/*_Calls.csv"
                    },
                    "id": "#GetDataTable.cwl/Tag_Calls"
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
                        "glob": "Trueno/*csv"
                    },
                    "id": "#GetDataTable.cwl/Trueno_out"
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
                        "glob": "Trueno/*zip"
                    },
                    "id": "#GetDataTable.cwl/Trueno_zip"
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
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [],
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
                    "id": "#InternalSettings.cwl/Subsample_Tags"
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
            "expression": "${\n  var internalInputs = [\n    '_Label_Version',\n    '_Read_Filter_Off',\n    '_Barcode_Num',\n    '_Seq_Run',\n    '_AbSeq_UMI',\n    '_Use_DBEC',\n    '_Extra_Seqs',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n    '_Target_analysis',\n    '_Subsample_Tags',\n    '_VDJ_VGene_Evalue',\n    '_VDJ_JGene_Evalue',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
            "id": "#InternalSettings.cwl"
        },
        {
            "class": "Workflow",
            "label": "BD Rhapsody\u2122 Targeted Analysis Pipeline 1.10.1 Large Input",
            "doc": "The BD Rhapsody\u2122 assays are used to create sequencing libraries from single cell transcriptomes.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files and a reference file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
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
                    "label": "Putative Cell Calling",
                    "doc": "Specify the data to be used for putative cell calling. mRNA is the default selected option. AbSeq (Experimental) is for troubleshooting only.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Putative_Cell_Call/Putative_Cell_Call",
                            "symbols": [
                                "#main/Putative_Cell_Call/Putative_Cell_Call/mRNA",
                                "#main/Putative_Cell_Call/Putative_Cell_Call/AbSeq_Experimental"
                            ]
                        }
                    ],
                    "id": "#main/Putative_Cell_Call"
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
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Reference",
                    "doc": "A fasta file containing the mRNA panel amplicon targets used in the experiment",
                    "id": "#main/Reference"
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
                                "#main/Sample_Tags_Version/Sample_Tags_Version/custom"
                            ]
                        }
                    ],
                    "id": "#main/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Subsample Reads",
                    "doc": "Any number of reads >1 or a fraction between 0 < n < 1 to indicate the percentage of reads to subsample.\n",
                    "id": "#main/Subsample"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Subsample Seed",
                    "doc": "For use when replicating a previous subsampling run only. Obtain the seed generated from the log file for the SplitFastQ node.\n",
                    "id": "#main/Subsample_seed"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "label": "Tag Names",
                    "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Should be alpha numeric, with + - and _ allowed. Any special characters: &, (), [], {}, <>, ?, | will be corrected to underscores. \n",
                    "id": "#main/Tag_Names"
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
                    "label": "Bioproduct Statistics",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Bioproduct_Stats",
                    "id": "#main/Bioproduct_Stats"
                },
                {
                    "label": "Cell Label Filter",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/GetDataTable/Cell_Label_Filter",
                    "id": "#main/Cell_Label_Filter"
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
                    "outputSource": "#main/CellClassifier/cellTypePredictions",
                    "id": "#main/ImmuneCellClassification(Experimental)"
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
                    "label": "vdjDominantContigs",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjDominantContigs",
                    "id": "#main/vdjDominantContigs"
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
                    "label": "vdjUnfilteredContigs",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjUnfilteredContigs",
                    "id": "#main/vdjUnfilteredContigs"
                }
            ],
            "steps": [
                {
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "coresMin": 8,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#AlignR2.cwl",
                    "in": [
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AlignR2/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/Index",
                            "id": "#main/AlignR2/Index"
                        },
                        {
                            "source": "#main/QualityFilterOuter/R2",
                            "id": "#main/AlignR2/R2"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AlignR2/Run_Metadata"
                        }
                    ],
                    "out": [
                        "#main/AlignR2/Alignments",
                        "#main/AlignR2/output"
                    ],
                    "id": "#main/AlignR2"
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
                            "source": "#main/AnnotateReads/Valid_Reads",
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
                    "run": "#AnnotateR1.cwl",
                    "in": [
                        {
                            "source": "#main/QualityFilterOuter/Filter_Metrics",
                            "id": "#main/AnnotateR1/Filter_Metrics"
                        },
                        {
                            "source": "#main/QualityFilterOuter/R1",
                            "id": "#main/AnnotateR1/R1"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateR1/Run_Metadata"
                        }
                    ],
                    "out": [
                        "#main/AnnotateR1/Annotation_R1",
                        "#main/AnnotateR1/R1_error_count_table",
                        "#main/AnnotateR1/R1_read_count_breakdown",
                        "#main/AnnotateR1/output"
                    ],
                    "scatter": [
                        "#main/AnnotateR1/R1"
                    ],
                    "id": "#main/AnnotateR1"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#AnnotateR2.cwl",
                    "in": [
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AnnotateR2/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/AnnotateR2/GTF_Annotation"
                        },
                        {
                            "source": "#main/AlignR2/Alignments",
                            "id": "#main/AnnotateR2/R2_zip"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateR2/Run_Metadata"
                        },
                        {
                            "source": "#main/CheckReference/Transcript_Length",
                            "id": "#main/AnnotateR2/Transcript_Length"
                        }
                    ],
                    "out": [
                        "#main/AnnotateR2/Annot_R2",
                        "#main/AnnotateR2/R2_Bam",
                        "#main/AnnotateR2/GTF",
                        "#main/AnnotateR2/output",
                        "#main/AnnotateR2/R2_Quality_Metrics"
                    ],
                    "scatter": [
                        "#main/AnnotateR2/R2_zip"
                    ],
                    "id": "#main/AnnotateR2"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 32000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#AnnotateReads.cwl",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/AbSeq_UMI",
                            "id": "#main/AnnotateReads/AbSeq_UMI"
                        },
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AnnotateReads/Extra_Seqs"
                        },
                        {
                            "source": "#main/QualityFilterOuter/Filter_Metrics",
                            "id": "#main/AnnotateReads/Filter_Metrics"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/AnnotateReads/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/AnnotateR1/Annotation_R1",
                            "id": "#main/AnnotateReads/R1_Annotation"
                        },
                        {
                            "source": "#main/AnnotateR1/R1_error_count_table",
                            "id": "#main/AnnotateReads/R1_error_count_table"
                        },
                        {
                            "source": "#main/AnnotateR1/R1_read_count_breakdown",
                            "id": "#main/AnnotateReads/R1_read_count_breakdown"
                        },
                        {
                            "source": "#main/AnnotateR2/Annot_R2",
                            "id": "#main/AnnotateReads/R2_Annotation"
                        },
                        {
                            "source": "#main/AnnotateR2/R2_Quality_Metrics",
                            "id": "#main/AnnotateReads/R2_Quality_Metrics"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateReads/Run_Metadata"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AnnotateReads/Target_Gene_Mapping"
                        }
                    ],
                    "out": [
                        "#main/AnnotateReads/Seq_Metrics",
                        "#main/AnnotateReads/Valid_Reads",
                        "#main/AnnotateReads/Read1_error_rate",
                        "#main/AnnotateReads/Annotation_Read",
                        "#main/AnnotateReads/output",
                        "#main/AnnotateReads/validTcrReads",
                        "#main/AnnotateReads/validIgReads",
                        "#main/AnnotateReads/num_valid_tcr_reads",
                        "#main/AnnotateReads/num_valid_ig_reads"
                    ],
                    "id": "#main/AnnotateReads"
                },
                {
                    "run": "#BundleLogs.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/AnnotateReads/output",
                                "#main/AnnotateR1/output",
                                "#main/AnnotateR2/output",
                                "#main/CheckReference/output",
                                "#main/GetDataTable/output",
                                "#main/Metrics/output",
                                "#main/AnnotateMolecules/output",
                                "#main/QualityFilterOuter/output",
                                "#main/CheckFastqs/log",
                                "#main/SplitAndSubsample/log",
                                "#main/Dense_to_Sparse_Datatable/output",
                                "#main/Dense_to_Sparse_Datatable_Unfiltered/output",
                                "#main/CellClassifier/log"
                            ],
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
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#Cell_Classifier.cwl",
                    "in": [
                        {
                            "source": "#main/FindDataTableForCellClassifier/molsPerCellMatrixForCellClassifier",
                            "id": "#main/CellClassifier/molsPerCellMatrix"
                        }
                    ],
                    "out": [
                        "#main/CellClassifier/cellTypePredictions",
                        "#main/CellClassifier/log"
                    ],
                    "id": "#main/CellClassifier"
                },
                {
                    "run": "#CheckFastqs.cwl",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/MinChunkSize",
                            "id": "#main/CheckFastqs/MinChunkSize"
                        },
                        {
                            "source": "#main/Reads",
                            "id": "#main/CheckFastqs/Reads"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Reads",
                            "id": "#main/CheckFastqs/Subsample"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Seed",
                            "id": "#main/CheckFastqs/Subsample_Seed"
                        }
                    ],
                    "out": [
                        "#main/CheckFastqs/SubsampleSeed",
                        "#main/CheckFastqs/SubsamplingRatio",
                        "#main/CheckFastqs/FilesToSkipSplitAndSubsample",
                        "#main/CheckFastqs/FastqReadPairs",
                        "#main/CheckFastqs/Bead_Version",
                        "#main/CheckFastqs/Libraries",
                        "#main/CheckFastqs/ReadsList",
                        "#main/CheckFastqs/log"
                    ],
                    "id": "#main/CheckFastqs"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 1000,
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
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/CheckReference/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/Reference",
                            "id": "#main/CheckReference/Reference"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/CheckReference/Run_Metadata"
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
                    "run": {
                        "cwlVersion": "v1.0",
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "id": "#main/FindDataTableForCellClassifier/3432aa72-4b1b-48ce-8eb0-bc1b2ed45ad2/dataTables"
                            }
                        ],
                        "outputs": [
                            {
                                "type": "File",
                                "id": "#main/FindDataTableForCellClassifier/3432aa72-4b1b-48ce-8eb0-bc1b2ed45ad2/molsPerCellMatrixForCellClassifier"
                            }
                        ],
                        "expression": "${\n  for (var i = 0; i < inputs.dataTables.length; i++) {\n    var dataTable = inputs.dataTables[i];\n    if (dataTable.basename.indexOf(\"_RSEC_MolsPerCell.csv\") >= 0) {\n      return({molsPerCellMatrixForCellClassifier: dataTable});\n    }\n  }\n  return({molsPerCellMatrixForCellClassifier: null});\n}",
                        "id": "#main/FindDataTableForCellClassifier/3432aa72-4b1b-48ce-8eb0-bc1b2ed45ad2"
                    },
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_Datatable/Data_Tables",
                            "id": "#main/FindDataTableForCellClassifier/dataTables"
                        }
                    ],
                    "out": [
                        "#main/FindDataTableForCellClassifier/molsPerCellMatrixForCellClassifier"
                    ],
                    "id": "#main/FindDataTableForCellClassifier"
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
                            "source": "#main/AnnotateReads/Seq_Metrics",
                            "id": "#main/GetDataTable/Seq_Metrics"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Tag_Sample_Names",
                            "id": "#main/GetDataTable/Tag_Names"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Total_Molecules",
                            "id": "#main/GetDataTable/Total_Molecules"
                        }
                    ],
                    "out": [
                        "#main/GetDataTable/Tag_Calls",
                        "#main/GetDataTable/Molecular_Annotation",
                        "#main/GetDataTable/Corrected_Molecular_Annotation",
                        "#main/GetDataTable/Tag_Annotation",
                        "#main/GetDataTable/Annot_Files",
                        "#main/GetDataTable/Cell_Label_Filter",
                        "#main/GetDataTable/Dense_Data_Tables",
                        "#main/GetDataTable/Dense_Data_Tables_Unfiltered",
                        "#main/GetDataTable/Expression_Data",
                        "#main/GetDataTable/Expression_Data_Unfiltered",
                        "#main/GetDataTable/Bioproduct_Stats",
                        "#main/GetDataTable/UMI_Adjusted_CellLabel_Stats",
                        "#main/GetDataTable/Putative_Cells_Origin",
                        "#main/GetDataTable/Protein_Aggregates_Experimental",
                        "#main/GetDataTable/Trueno_out",
                        "#main/GetDataTable/Trueno_zip",
                        "#main/GetDataTable/output",
                        "#main/GetDataTable/Cell_Order",
                        "#main/GetDataTable/Gene_List"
                    ],
                    "id": "#main/GetDataTable"
                },
                {
                    "label": "Internal Settings",
                    "run": "#InternalSettings.cwl",
                    "in": [],
                    "out": [
                        "#main/Internal_Settings/Read_Filter_Off",
                        "#main/Internal_Settings/Barcode_Num",
                        "#main/Internal_Settings/Seq_Run",
                        "#main/Internal_Settings/AbSeq_UMI",
                        "#main/Internal_Settings/Use_DBEC",
                        "#main/Internal_Settings/Extra_Seqs",
                        "#main/Internal_Settings/MinChunkSize",
                        "#main/Internal_Settings/NumRecordsPerSplit",
                        "#main/Internal_Settings/Target_analysis",
                        "#main/Internal_Settings/Subsample_Tags",
                        "#main/Internal_Settings/VDJ_VGene_Evalue",
                        "#main/Internal_Settings/VDJ_JGene_Evalue"
                    ],
                    "id": "#main/Internal_Settings"
                },
                {
                    "run": {
                        "cwlVersion": "v1.0",
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
                                "id": "#main/MergeMultiplex/c48ecf06-344f-42d3-b8ec-f9d2f4d5cdb9/SampleTag_Files"
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
                                "id": "#main/MergeMultiplex/c48ecf06-344f-42d3-b8ec-f9d2f4d5cdb9/Multiplex_out"
                            }
                        ],
                        "expression": "${\n  var fp_array = [];\n  for (var i = 0; i < inputs.SampleTag_Files.length; i++) {\n    var fp = inputs.SampleTag_Files[i];\n    if (fp != null) {\n      fp_array.push(fp);\n    }\n  }\n  return({\"Multiplex_out\": fp_array});\n}",
                        "id": "#main/MergeMultiplex/c48ecf06-344f-42d3-b8ec-f9d2f4d5cdb9"
                    },
                    "in": [
                        {
                            "source": [
                                "#main/GetDataTable/Trueno_out",
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
                            "valueFrom": "Targeted",
                            "id": "#main/Metadata_Settings/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Basic_Algo_Only",
                            "id": "#main/Metadata_Settings/Basic_Algo_Only"
                        },
                        {
                            "source": "#main/CheckFastqs/Bead_Version",
                            "id": "#main/Metadata_Settings/Bead_Version"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                            "id": "#main/Metadata_Settings/Exact_Cell_Count"
                        },
                        {
                            "source": "#main/CheckFastqs/Libraries",
                            "id": "#main/Metadata_Settings/Libraries"
                        },
                        {
                            "valueFrom": "BD Rhapsody Targeted Analysis Pipeline",
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
                            "source": "#main/CheckFastqs/ReadsList",
                            "id": "#main/Metadata_Settings/Reads"
                        },
                        {
                            "source": "#main/Reference",
                            "id": "#main/Metadata_Settings/Reference"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Name",
                            "id": "#main/Metadata_Settings/Run_Name"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Tag_Sample_Names",
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
                            "source": "#main/Subsample_Settings/Subsample_Reads",
                            "id": "#main/Metadata_Settings/Subsample"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Seed",
                            "id": "#main/Metadata_Settings/Subsample_Seed"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/Metadata_Settings/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/Metadata_Settings/Run_Metadata",
                        "#main/Metadata_Settings/Run_Base_Name"
                    ],
                    "id": "#main/Metadata_Settings"
                },
                {
                    "run": "#Metrics.cwl",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Annot_Files",
                            "id": "#main/Metrics/Annot_Files"
                        },
                        {
                            "source": "#main/AnnotateReads/Read1_error_rate",
                            "id": "#main/Metrics/Read1_error_rate"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Metrics/Run_Metadata"
                        },
                        {
                            "source": "#main/GetDataTable/Trueno_zip",
                            "id": "#main/Metrics/Sample_Tag_Archives"
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
                            "source": "#main/VDJ_Compile_Results/vdjMetricsJson",
                            "id": "#main/Metrics/vdjMetricsJson"
                        }
                    ],
                    "out": [
                        "#main/Metrics/Metrics_Summary",
                        "#main/Metrics/Metrics_Archive",
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
                            "source": "#main/Sample_Tags_Version",
                            "id": "#main/Multiplexing_Settings/_Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Tag_Names",
                            "id": "#main/Multiplexing_Settings/_Tag_Sample_Names"
                        }
                    ],
                    "out": [
                        "#main/Multiplexing_Settings/Tag_Sample_Names",
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
                    "run": "#PairReadFiles.cwl",
                    "in": [
                        {
                            "source": "#main/CheckFastqs/FastqReadPairs",
                            "id": "#main/PairReadFiles/FastqReadPairs"
                        },
                        {
                            "source": "#main/SplitAndSubsample/SplitAndSubsampledFastqs",
                            "id": "#main/PairReadFiles/Reads"
                        }
                    ],
                    "out": [
                        "#main/PairReadFiles/ReadPairs"
                    ],
                    "id": "#main/PairReadFiles"
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
                    "run": "#QualityFilterOuter.cwl",
                    "in": [
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/QualityFilterOuter/Run_Metadata"
                        },
                        {
                            "source": "#main/PairReadFiles/ReadPairs",
                            "id": "#main/QualityFilterOuter/Split_Read_Pairs"
                        }
                    ],
                    "out": [
                        "#main/QualityFilterOuter/Filter_Metrics",
                        "#main/QualityFilterOuter/R1",
                        "#main/QualityFilterOuter/R2",
                        "#main/QualityFilterOuter/output"
                    ],
                    "id": "#main/QualityFilterOuter"
                },
                {
                    "run": "#SplitAndSubsample.cwl",
                    "in": [
                        {
                            "source": "#main/Reads",
                            "id": "#main/SplitAndSubsample/Fastqs"
                        },
                        {
                            "source": "#main/CheckFastqs/FilesToSkipSplitAndSubsample",
                            "id": "#main/SplitAndSubsample/FilesToSkipSplitAndSubsample"
                        },
                        {
                            "source": "#main/Internal_Settings/NumRecordsPerSplit",
                            "id": "#main/SplitAndSubsample/NumRecordsPerSplit"
                        },
                        {
                            "source": "#main/CheckFastqs/SubsamplingRatio",
                            "id": "#main/SplitAndSubsample/SubsampleRatio"
                        },
                        {
                            "source": "#main/CheckFastqs/SubsampleSeed",
                            "id": "#main/SplitAndSubsample/SubsampleSeed"
                        }
                    ],
                    "out": [
                        "#main/SplitAndSubsample/SplitAndSubsampledFastqs",
                        "#main/SplitAndSubsample/log"
                    ],
                    "id": "#main/SplitAndSubsample"
                },
                {
                    "run": {
                        "cwlVersion": "v1.0",
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
                                "id": "#main/Start_Time/037fec44-4671-4072-93e3-e7d0bf5ddb73/Start_Time"
                            }
                        ],
                        "expression": "${   \n  var today = new Date();\n  var date = today.toString()\n  return ({Start_Time: date});\n} ",
                        "id": "#main/Start_Time/037fec44-4671-4072-93e3-e7d0bf5ddb73"
                    },
                    "in": [],
                    "out": [
                        "#main/Start_Time/Start_Time"
                    ],
                    "id": "#main/Start_Time"
                },
                {
                    "label": "Subsample Settings",
                    "run": "#SubsampleSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Subsample",
                            "id": "#main/Subsample_Settings/_Subsample_Reads"
                        },
                        {
                            "source": "#main/Subsample_seed",
                            "id": "#main/Subsample_Settings/_Subsample_Seed"
                        }
                    ],
                    "out": [
                        "#main/Subsample_Settings/Subsample_Reads",
                        "#main/Subsample_Settings/Subsample_Seed"
                    ],
                    "id": "#main/Subsample_Settings"
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
                            "source": "#main/AnnotateReads/Seq_Metrics",
                            "id": "#main/VDJ_Compile_Results/Seq_Metrics"
                        },
                        {
                            "source": "#main/CellClassifier/cellTypePredictions",
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
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Compile_Results/vdjVersion"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Compile_Results/vdjCellsDatatable",
                        "#main/VDJ_Compile_Results/vdjCellsDatatableUncorrected",
                        "#main/VDJ_Compile_Results/vdjDominantContigs",
                        "#main/VDJ_Compile_Results/vdjUnfilteredContigs",
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
                            "source": "#main/AnnotateReads/validIgReads",
                            "id": "#main/VDJ_Preprocess_Reads_IG/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#main/AnnotateReads/num_valid_ig_reads",
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
                            "source": "#main/AnnotateReads/validTcrReads",
                            "id": "#main/VDJ_Preprocess_Reads_TCR/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#main/AnnotateReads/num_valid_tcr_reads",
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
                    "id": "#Metadata.cwl/Subsample"
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
                        "float"
                    ],
                    "id": "#Metadata.cwl/Subsample_Tags"
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
                    "id": "#Metadata.cwl/Run_Metadata",
                    "outputBinding": {
                        "glob": "run_metadata.json"
                    }
                }
            ],
            "stdout": "run_metadata.json",
            "arguments": [
                {
                    "prefix": ""
                },
                {
                    "shellQuote": true,
                    "valueFrom": "${\n  var metadata = inputs;\n  var all_bv = {};\n  var customer_bv = \"Original (V1)\";\n  for (var i = 0; i < inputs.Bead_Version.length; i++) {\n      var BeadVer = inputs.Bead_Version[i];\n      var Library = BeadVer[\"Library\"];\n      var bead_version = BeadVer[\"bead_version\"];\n      all_bv[Library] = bead_version  \n      var short_bv =  bead_version.substring(0, 2);\n      if (short_bv == \"V2\"){\n        var customer_bv = \"Enhanced (V2)\";\n      }\n  }\n  metadata[\"Bead_Version\"] = all_bv;\n\n  var pipeline_name = inputs.Pipeline_Name;\n  var assay = inputs.Assay;\n  var version = inputs.Pipeline_Version;\n  var time = inputs.Start_Time;\n  var libraries = inputs.Libraries.split(\",\");\n  var i = 0;\n  var reference_list = []\n  if(inputs.Reference != null){\n      reference_list = reference_list.concat(inputs.Reference);\n  }\n  if(inputs.AbSeq_Reference != null){\n      reference_list = reference_list.concat(inputs.AbSeq_Reference);\n  }\n\n  var supplemental = \"\"\n  if(inputs.Supplemental_Reference != null){\n      supplemental = \"; Supplemental_Reference - \" + inputs.Supplemental_Reference[0][\"basename\"];\n  }\n  var references = [];\n  for (i = 0; i< reference_list.length; i++) {\n      if(reference_list[i] != null){\n          references.push(reference_list[i][\"basename\"]);\n      }\n  }\n  var parameters = [];\n  if(inputs.Sample_Tags_Version != null){\n      var tags = \"Sample Tag Version: \" + inputs.Sample_Tags_Version;\n  } else{ \n      var tags = \"Sample Tag Version: None\";\n  }\n  parameters.push(tags);\n\n  if(inputs.Sample_Tag_Names != null){\n      var tag_names = inputs.Sample_Tag_Names.join(\" ; \")\n      var tag_list = \"Sample Tag Names: \" + tag_names;\n  } else{\n      var tag_list = \"Sample Tag Names: None\";\n  }\n  parameters.push(tag_list);\n \n  if(inputs.VDJ_Version != null){\n      var vdj = \"VDJ Version: \" + inputs.VDJ_Version;\n  } else{ \n      var vdj = \"VDJ Version: None\";\n  }\n  parameters.push(vdj)\n\n  if(inputs.Subsample != null){\n      var subsample = \"Subsample: \" + inputs.Subsample;\n  } else{ \n      var subsample = \"Subsample: None\";\n  }   \n  parameters.push(subsample);\n\n  if(inputs.Putative_Cell_Call == 1){\n      var call = \"Putative Cell Calling Type: AbSeq\";\n  } else{ \n      var call = \"Putative Cell Calling Type: mRNA\";\n  }   \n  parameters.push(call)\n\n  if(inputs.Basic_Algo_Only){\n      var basic = \"Refined Putative Cell Calling: Off\";\n  } else{ \n      var basic = \"Refined Putative Cell Calling: On\";\n  }   \n  parameters.push(basic)\n\n  if(inputs.Exact_Cell_Count != null){\n      var cells = \"Exact Cell Count: \" + inputs.Exact_Cell_Count;\n  } else{ \n      var cells = \"Exact Cell Count: None\";\n  }   \n  parameters.push(cells)\n\n  var name = inputs.Run_Name;\n  if (name == null){\n    var libraries = inputs.Libraries.split(',');\n    name = libraries[0];\n  }        \n\n  var header = [\"####################\"];\n  header.push(\"## \" + pipeline_name + \" Version \" + version);\n  header.push(\"## Analysis Date - \" + time);\n  header.push(\"## Libraries - \" + libraries.join(' | ') + \" - Bead version detected: \" + customer_bv);\n  header.push(\"## References - \" + references.join(' | ') + supplemental);\n  header.push(\"## Parameters - \" + parameters.join(' | '));\n  header.push(\"####################\");\n  metadata[\"Output_Header\"] = header;\n  metadata[\"Run_Base_Name\"] = name;\n  var metadata_json = JSON.stringify(metadata);\n  return metadata_json;\n}\n"
                }
            ],
            "id": "#Metadata.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--read1-error-rate"
                    },
                    "id": "#Metrics.cwl/Read1_error_rate"
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
                        "Any"
                    ],
                    "id": "#MultiplexingSettings.cwl/_Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/_Tag_Sample_Names"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MultiplexingSettings.cwl/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/Tag_Sample_Names"
                }
            ],
            "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  var listTagNames = inputs._Tag_Sample_Names\n  var newTagNames = []\n  for (var num in listTagNames) {\n    var tag = listTagNames[num].replace(/[^A-Za-z0-9-+]/g,\"_\");\n    newTagNames.push(tag);    \n  }  \n  return ({\n  Tag_Sample_Names: newTagNames,\n  Sample_Tags_Version: enumifiedSampleTagsVersion\n  });\n}",
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
            "expression": "${ var name = inputs._Run_Name;\n   if (name != null) {\n     name = name.replace(/[\\W_]+/g,\"-\");}\n   return({'Run_Name' : name });\n }  ",
            "id": "#NameSettings.cwl"
        },
        {
            "doc": "PairReadFiles takes an array of split files and pairs them, such that an R1 file is transferred to the QualityFilter with its corresponding R2 file.\nThe original FASTQ files are paired in CheckFastqs and then split and sub-sampled in SplitAndSubsample. The pairing information is taken from CheckFastqs.\n",
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/filename",
                                    "type": "string"
                                },
                                {
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/readFlag",
                                    "type": "string"
                                },
                                {
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/readPairId",
                                    "type": "string"
                                },
                                {
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/library",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "id": "#PairReadFiles.cwl/FastqReadPairs"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#PairReadFiles.cwl/Reads"
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
                                    "name": "#PairReadFiles.cwl/ReadPairs/R1",
                                    "type": "File"
                                },
                                {
                                    "name": "#PairReadFiles.cwl/ReadPairs/R2",
                                    "type": "File"
                                },
                                {
                                    "name": "#PairReadFiles.cwl/ReadPairs/readPairId",
                                    "type": "int"
                                },
                                {
                                    "name": "#PairReadFiles.cwl/ReadPairs/library",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "id": "#PairReadFiles.cwl/ReadPairs"
                }
            ],
            "expression": "${\n  // use the CheckFastqs read pairing information to create a dictionary\n  // using the original fastq file name without the extension as the key\n  var fastqReadPairs = {}\n  for (var i = 0; i < inputs.FastqReadPairs.length; i++) {\n    var fileDict = inputs.FastqReadPairs[i];\n    var filename = fileDict[\"filename\"];\n\n    if (!fastqReadPairs[filename]) {\n      fastqReadPairs[filename] = {\n        readPairId: null,\n        readFlag: null,\n        library: null,\n      };\n    }\n    else {\n      throw new Error(\"Found non-unique fastq filename '\" + filename + \"' in the FastqReadPairs dictionary from CheckFastqs.\")\n    }\n\n    fastqReadPairs[filename].readPairId = fileDict[\"readPairId\"]\n    fastqReadPairs[filename].readFlag = fileDict[\"readFlag\"]\n    fastqReadPairs[filename].library = fileDict[\"library\"]\n  }\n\n  // now loop through the input read files which could\n  // be the original fastq files if no sub-sampling has\n  // been done, or the sub-sampled fastq files\n  var readPairs = {}\n  for (var i = 0; i < inputs.Reads.length; i++) {\n\n    // Set the fileDict to null\n    var fileDict = null;\n\n    // Get the fastq file\n    var fastqFile = inputs.Reads[i];\n\n    // Remove the .gz from the end of the filename\n    var fileNoGzExt = fastqFile.basename.replace(/.gz$/i, \"\");\n\n    // Remove the next file extension if it exists\n    var fileArrayWithExt = fileNoGzExt.split(\".\");\n    // If an extension exists, splice the array\n    var fileArrayNoExt = null;\n    if (fileArrayWithExt.length > 1) {\n      fileArrayNoExt = fileArrayWithExt.splice(0, fileArrayWithExt.length-1);\n    } else {\n      // No file extension exists, so use the whole array\n      fileArrayNoExt = fileArrayWithExt\n    }\n    var fileRootname = fileArrayNoExt.join(\".\")\n\n    // if the original files were sub-sampled\n    // get the original file and the chunk id\n    if (fileRootname.indexOf(\"-\") != -1) {\n      // Split on the dash to get the name of\n      // the original file and the chunk id\n      // The original file name can also have dashes\n      var chunkFileArray = fileRootname.split(\"-\");\n\n      // Get the original file rootname and chunk id\n      // The rootname without the chunk id and file\n      // extension is the key from CheckFastqs\n      // The chunk id is used later to create a new unique\n      // read pair id for all sub-sampled fastq files\n\n      // The rootname array should contain all elements up to the last dash\n      var fileRootnameArray = chunkFileArray.splice(0, chunkFileArray.length-1);\n      var fileRootnameNoChunkId = fileRootnameArray.join(\"-\");\n\n      // The chunk id is the last element in the array\n      // representing the content after the last dash\n      var orgChunkId = chunkFileArray.pop();\n\n      // if there is no chunk id, use an arbitrary number\n      // the chunk id is unique when the files are sub-sampled\n      // and does not need to be unique when the files are not sub-sampled\n      var chunkId = 9999;\n      if (orgChunkId) {\n        // cast to an integer\n        chunkId = parseInt(orgChunkId);\n      }\n      // double check that we have a chunk id\n      if (chunkId === undefined || chunkId === null) {\n        throw new Error(\"The fastq file sub-sampling id could not be determined!\");\n      }\n\n      // The file rootname without the chunk id and file extension\n      // should match the original file rootname from CheckFastqs\n      // The original file rootname from CheckFastqs is the key for\n      // the dictionary containing the original unique pair id\n      var fileDict = fastqReadPairs[fileRootnameNoChunkId];\n    }\n\n    // If the files are not sub-sampled or the fileDict\n    // is not found, then try to use the original\n    // file rootname without the file extension as the key\n    if (fileDict === undefined || fileDict === null) {\n\n      // if the original files were not sub-sampled,\n      // use the original file rootname and an arbitrary chunk id\n      var chunkId = 9999;\n\n      var fileDict = fastqReadPairs[fileRootname];\n\n      // If the fileDict for this file rootname is not found,\n      // then the filenames are in an unexpected format and\n      // the code to parse the filenames in CheckFastqs,\n      // SplitAndSubsample and here need to match\n      if (fileDict === undefined || fileDict === null) {\n        // Create an error\n        if (fileDict === undefined || fileDict === null) {\n          throw new Error(\"Cannot find the fastq read pair information for '\" + fastqFile.basename + \"'.\");\n        }\n      }\n    }\n\n    // Get the pairing information from CheckFastqs\n    var readPairId = fileDict[\"readPairId\"];\n    var library = fileDict[\"library\"];\n    var flag = fileDict[\"readFlag\"];\n\n    // Add the chunkId to create a new unique read pair id\n    // for each file (sub-sampled or not)\n    var chunkReadPairId = readPairId + \"_\" + chunkId;\n\n    // Create a dictionary for each pair of files\n    if (!readPairs[chunkReadPairId]) {\n      readPairs[chunkReadPairId] = {\n        R1: null,\n        R2: null,\n        library: library,\n        readPairId: null,\n      };\n    }\n    // add in the R1 and R2 files, depending on the flag\n    if (flag === \"R1\") {\n      readPairs[chunkReadPairId].R1 = fastqFile\n    } else if (flag === \"R2\") {\n      readPairs[chunkReadPairId].R2 = fastqFile\n    }\n  }\n  // we are not interested in the read pair ids in readPairs\n  // flatten into an array of objects\n  var readPairsList = [];\n  var i = 1;\n  for (var key in readPairs) {\n    if (readPairs.hasOwnProperty(key)) {\n      var readPair = readPairs[key];\n      readPair.readPairId = i;\n      readPairsList.push(readPair);\n      i++;\n    }\n  }\n  // pass this array to the record array named \"ReadPairs\" on the CWL layer\n  return {ReadPairs: readPairsList}\n}",
            "id": "#PairReadFiles.cwl"
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
            "expression": "${\n  // the basic algorithm flag defaults to false\n  var basicAlgOnlyFlag = false;\n  // the user can set the basic algorithm flag\n  if (inputs._Basic_Algo_Only) {\n    basicAlgOnlyFlag = inputs._Basic_Algo_Only;\n  }\n  // convert the Putative_Cell_Call from a string to an integer\n  var putativeCellCallInt = 0;\n  if (inputs._Putative_Cell_Call) {\n    if (inputs._Putative_Cell_Call === \"mRNA\") {\n      putativeCellCallInt = 0;\n    }\n    else if (inputs._Putative_Cell_Call == \"AbSeq_Experimental\" || inputs._Putative_Cell_Call == \"AbSeq (Experimental)\") {\n      putativeCellCallInt = 1;\n      // for protein-only cell calling, we only have the basic algorithm\n      basicAlgOnlyFlag = true;\n    }\n    else if (inputs._Putative_Cell_Call == \"mRNA_and_AbSeq\") {\n      putativeCellCallInt = 2;\n    }\n  }\n  // check the exact cell count\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      throw(\"Illogical value for exact cell count: \" + inputs._Exact_Cell_Count);\n    }\n  }\n  return ({\n    Putative_Cell_Call: putativeCellCallInt,\n    Exact_Cell_Count: inputs._Exact_Cell_Count,\n    Basic_Algo_Only: basicAlgOnlyFlag,\n  });\n}",
            "id": "#PutativeCellSettings.cwl"
        },
        {
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_quality_filter.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#QualityFilter.cwl/Run_Metadata"
                },
                {
                    "type": {
                        "type": "record",
                        "fields": [
                            {
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/R1",
                                "type": "File",
                                "inputBinding": {
                                    "prefix": "--r1"
                                }
                            },
                            {
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/R2",
                                "type": "File",
                                "inputBinding": {
                                    "prefix": "--r2"
                                }
                            },
                            {
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/readPairId",
                                "type": "int",
                                "inputBinding": {
                                    "prefix": "--read-pair-id"
                                }
                            },
                            {
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/library",
                                "type": "string",
                                "inputBinding": {
                                    "prefix": "--library"
                                }
                            }
                        ]
                    },
                    "id": "#QualityFilter.cwl/Split_Read_Pairs"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*read_quality.csv.gz"
                    },
                    "id": "#QualityFilter.cwl/Filter_Metrics"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_R1*.fastq.gz"
                    },
                    "id": "#QualityFilter.cwl/R1"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_R2*.fastq.gz"
                    },
                    "id": "#QualityFilter.cwl/R2"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#QualityFilter.cwl/output"
                }
            ],
            "id": "#QualityFilter.cwl"
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
                    "type": "File",
                    "id": "#QualityFilterOuter.cwl/Run_Metadata"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/R1",
                                    "type": "File"
                                },
                                {
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/R2",
                                    "type": "File"
                                },
                                {
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/readPairId",
                                    "type": "int"
                                },
                                {
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/library",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "id": "#QualityFilterOuter.cwl/Split_Read_Pairs"
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
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Filter_Metrics",
                    "id": "#QualityFilterOuter.cwl/Filter_Metrics"
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
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R1",
                    "id": "#QualityFilterOuter.cwl/R1"
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
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R2",
                    "id": "#QualityFilterOuter.cwl/R2"
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
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/output",
                    "id": "#QualityFilterOuter.cwl/output"
                }
            ],
            "steps": [
                {
                    "id": "#QualityFilterOuter.cwl/Quality_Filter_Scatter",
                    "run": "#QualityFilter.cwl",
                    "in": [
                        {
                            "source": "#QualityFilterOuter.cwl/Run_Metadata",
                            "id": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Run_Metadata"
                        },
                        {
                            "source": "#QualityFilterOuter.cwl/Split_Read_Pairs",
                            "id": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Split_Read_Pairs"
                        }
                    ],
                    "out": [
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R1",
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R2",
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Filter_Metrics",
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/output"
                    ],
                    "scatter": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Split_Read_Pairs"
                }
            ],
            "id": "#QualityFilterOuter.cwl"
        },
        {
            "doc": "SplitAndSubsample splits, subsamples and formats read files to be deposited in QualityFilter.\n",
            "class": "Workflow",
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#SplitAndSubsample.cwl/Fastqs"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "id": "#SplitAndSubsample.cwl/FilesToSkipSplitAndSubsample"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#SplitAndSubsample.cwl/NumRecordsPerSplit"
                },
                {
                    "type": "float",
                    "id": "#SplitAndSubsample.cwl/SubsampleRatio"
                },
                {
                    "type": "int",
                    "id": "#SplitAndSubsample.cwl/SubsampleSeed"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#SplitAndSubsample.cwl/FlattenOutput/SplitFastqList",
                    "id": "#SplitAndSubsample.cwl/SplitAndSubsampledFastqs"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#SplitAndSubsample.cwl/SplitAndSubsample/log",
                    "id": "#SplitAndSubsample.cwl/log"
                }
            ],
            "steps": [
                {
                    "doc": "After scattering \"SplitAndSubsample\" on a File array, the output of each node is also an array. Thus, we are left with a nestled list. This JS expression flattens this list to deal with the split reads in PairReadFiles.cwl",
                    "in": [
                        {
                            "source": "#SplitAndSubsample.cwl/SplitAndSubsample/SplitAndSubsampledFastqs",
                            "id": "#SplitAndSubsample.cwl/FlattenOutput/nestledSplitFastqList"
                        }
                    ],
                    "out": [
                        "#SplitAndSubsample.cwl/FlattenOutput/SplitFastqList"
                    ],
                    "run": {
                        "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output",
                        "cwlVersion": "v1.0",
                        "class": "ExpressionTool",
                        "inputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": {
                                        "type": "array",
                                        "items": "File"
                                    }
                                },
                                "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output/nestledSplitFastqList"
                            }
                        ],
                        "outputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output/SplitFastqList"
                            }
                        ],
                        "expression": "${\n  return {SplitFastqList: [].concat.apply([], inputs.nestledSplitFastqList)}\n}\n"
                    },
                    "id": "#SplitAndSubsample.cwl/FlattenOutput"
                },
                {
                    "doc": "Allocate one docker/python process per file to do the actual file splitting.",
                    "in": [
                        {
                            "source": "#SplitAndSubsample.cwl/Fastqs",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/Fastq"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/FilesToSkipSplitAndSubsample",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/FilesToSkipSplitAndSubsample"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/NumRecordsPerSplit",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/NumRecordsPerSplit"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/SubsampleRatio",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/SubsampleRatio"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/SubsampleSeed",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/SubsampleSeed"
                        }
                    ],
                    "scatter": [
                        "#SplitAndSubsample.cwl/SplitAndSubsample/Fastq"
                    ],
                    "out": [
                        "#SplitAndSubsample.cwl/SplitAndSubsample/SplitAndSubsampledFastqs",
                        "#SplitAndSubsample.cwl/SplitAndSubsample/log"
                    ],
                    "run": {
                        "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq",
                        "class": "CommandLineTool",
                        "cwlVersion": "v1.0",
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                                "class": "DockerRequirement"
                            }
                        ],
                        "baseCommand": [
                            "mist_split_fastq.py"
                        ],
                        "inputs": [
                            {
                                "type": "File",
                                "inputBinding": {
                                    "prefix": "--fastq-file-path"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/Fastq"
                            },
                            {
                                "type": {
                                    "type": "array",
                                    "items": "string"
                                },
                                "inputBinding": {
                                    "prefix": "--files-to-skip-split-and-subsample",
                                    "itemSeparator": ","
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/FilesToSkipSplitAndSubsample"
                            },
                            {
                                "type": [
                                    "null",
                                    "long"
                                ],
                                "inputBinding": {
                                    "prefix": "--num-records"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/NumRecordsPerSplit"
                            },
                            {
                                "type": "float",
                                "inputBinding": {
                                    "prefix": "--subsample-ratio"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SubsampleRatio"
                            },
                            {
                                "type": "int",
                                "inputBinding": {
                                    "prefix": "--subsample-seed"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SubsampleSeed"
                            }
                        ],
                        "outputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "outputBinding": {
                                    "glob": "*.fastq.gz",
                                    "outputEval": "${ if (self.length === 0) { return [inputs.Fastq]; } else { return self; } }"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SplitAndSubsampledFastqs"
                            },
                            {
                                "type": "File",
                                "outputBinding": {
                                    "glob": "*.log"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/log"
                            }
                        ]
                    },
                    "id": "#SplitAndSubsample.cwl/SplitAndSubsample"
                }
            ],
            "id": "#SplitAndSubsample.cwl"
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
                        "float"
                    ],
                    "id": "#SubsampleSettings.cwl/_Subsample_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SubsampleSettings.cwl/_Subsample_Seed"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#SubsampleSettings.cwl/Subsample_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SubsampleSettings.cwl/Subsample_Seed"
                }
            ],
            "expression": "${\n  var subsamplingOutputs = {\n    Subsample_Reads: inputs._Subsample_Reads,\n    Subsample_Seed: inputs._Subsample_Seed\n  }\n  return subsamplingOutputs;\n}",
            "id": "#SubsampleSettings.cwl"
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
                                "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                                "class": "DockerRequirement"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "cwlVersion": "v1.0",
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
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner/Compressed_File"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "type": "File",
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner/Uncompressed_File"
                            }
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner"
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
                                "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                                "class": "DockerRequirement"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "cwlVersion": "v1.0",
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
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner/Compressed_File"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "type": "File",
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner/Uncompressed_File"
                            }
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner"
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
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
                    "id": "#VDJ_Compile_Results.cwl/vdjVersion"
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
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_Dominant_Contigs.csv.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjDominantContigs"
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
                    "outputBinding": {
                        "glob": "*_VDJ_Unfiltered_Contigs.csv.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjUnfilteredContigs"
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
                        "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR",
                        "cwlVersion": "v1.0",
                        "class": "CommandLineTool",
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.10.1L",
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
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR/theCalls"
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
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR/gatheredCalls"
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
                        "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits",
                        "cwlVersion": "v1.0",
                        "class": "ExpressionTool",
                        "inputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/num_valid_reads"
                            },
                            {
                                "type": "string",
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/vdj_type"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/num_cores"
                            },
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/num_splits"
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
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": "mist_vdj_rsec_reads.py",
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
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": "VDJ_Trim_Reads.sh",
            "inputs": [
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
                    "dockerPull": "bdgenomics/rhapsody:1.10.1L",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_version.py"
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
    "cwlVersion": "v1.0",
    "$namespaces": {
        "sbg": "https://sevenbridges.com#",
        "arv": "http://arvados.org/cwl#"
    }
}