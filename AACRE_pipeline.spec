{
  "all":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "0:05:00",
    "queue": "genD",
    "qos": "test"
  },
  "trim_galore":{
    "name": "{rule}.{wildcards.file}.{pipeline_version}",
    "time": "3:00:00",
    "queue": "genD",
    "qos": "short"
  },
  "concat_reads":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "3:00:00",
    "queue": "genD",
    "qos": "short"
  },
  "Kraken2":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:15:00",
    "queue":"genD",
    "qos": "vshort"
  },
  "Trim_ONT_qcat":{
    "name": "{rule}.{wildcards.ontfile}.{pipeline_version}",
    "time": "3:00:00",
    "queue": "genD",
    "qos": "short"
  },
  "nanoplot":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:30:00",
    "queue": "genD",
    "qos": "vshort",
    "mem" : "500"
  },  
  "Verify_Species":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:05:00",
    "queue": "genD",
    "qos": "test"
  },
  "unicycler":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "22:00:00",
    "qos" : "long",
    "queue": "genD",
    "mem": 1024
  },
  "produce_final_assembly":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:10:00",
    "queue": "genD",
    "qos": "vshort"
  },
  "MLST":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "queue": "genD",
    "time": "0:05:00",
    "qos": "test"
  },
  "centrifuge":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "queue": "genD",
    "time": "3:00:00",
    "qos": "normal"
  },
  "Prokka":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genD",
    "qos": "vshort"
  },
  "RGI":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genD",
    "qos": "vshort"
  },
  "AMRFinder":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genD",
    "qos": "vshort"
  },
  "resfinder":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genD",
    "qos": "vshort"
  },
  "IS_Finder":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genD",
    "qos": "vshort"
  }
}
