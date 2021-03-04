{
  "all":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "0:05:00",
    "queue": "genB,main"
  },
  "Trim_Galore":{
    "name": "{rule}.{wildcards.file}.{pipeline_version}",
    "time": "3:00:00"
  },
  "Concat_Illumina":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "3:00:00"
  },
  "Bracken_Illumina":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:05:00"
  },
  "Trim_ONT_qcat":{
    "name": "{rule}.{wildcards.ontfile}.{pipeline_version}",
    "time": "3:00:00"
  },
  "Concat_ONT":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "3:00:00"
  },
  "NanoStats_ONT":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:10:00",
    "queue": "genB,main"
  },  
  "Bracken_ONT":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:05:00"
  },
  "Verify_Species":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "00:01:00",
    "queue": "genB,main"
  },
  "Assembly":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "22:00:00",
    "qos" : "normal",
    "queue": "main",
  },
  "MLST":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "queue": "genB,main",
    "time": "0:05:00" 
  },
  "centrifuge":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "queue": "genB,main",
    "time": "3:00:00",
    "qos": "normal"
  },
  "Annotation":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genB,main"
  },
  "AMRFinder":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genB,main"
  },
  "resfinder":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genB,main"
  },
  "IS_Finder":{
    "name": "{rule}.{sample}.{pipeline_version}",
    "time": "1:00:00",
    "queue": "genB,main"
  }
}
