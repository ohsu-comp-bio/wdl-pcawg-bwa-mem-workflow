task get_basename {
  File f

  command {
    basename ${f}
  }

  output {
    String base = read_string(stdout())
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task read_header {
  File unalignedBam
  String bamName
  String outputDir

  command {
    samtools view -H ${unalignedBam} | \
    perl -nae 'next unless /^\@RG/; s/\tPI:\t/\t/; s/\tPI:\s*\t/\t/; s/\tPI:\s*\z/\n/; s/\t/\\t/g; print' > "${outputDir}/${bamName}_header.txt"
  }

  output {
    String header = read_string("${outputDir}/${bamName}_header.txt")
    File header_file = "${outputDir}/${bamName}_header.txt"
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task count_reads {
  File unalignedBam
  String bamName
  String outputDir

  command {
    samtools view ${unalignedBam} | \
    wc -l > "${outputDir}/${bamName}_read_count.txt"
  }

  output {
    File counts_file = "${outputDir}/${bamName}_read_count.txt"
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task align {
  File unalignedBam
  String bamHeader
  File reference_gz
  File reference_gz_fai
  File reference_gz_amb
  File reference_gz_ann
  File reference_gz_bwt
  File reference_gz_pac
  File reference_gz_sa
  String bamName
  String outputDir
  Int threads
  Int sortMemMb

  command {
    bamtofastq exlcude=QCFAIL,SECONDARY,SUPPLEMENTARY T=${bamName + ".t"} S=${bamName + ".s"} O=${bamName + ".o"} O2=${bamName + ".o2"} collate=1 tryoq=1 filename=${unalignedBam} | \
    bwa mem -p -t ${threads} -T 0 -R "${bamHeader}" ${reference_gz} - | \
    bamsort blockmb=${sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${bamName + ".sorttmp"} O=${outputDir + "/" + bamName + "_aligned.bam"}
  }

  output {
    File bam_output = "${outputDir}/${bamName}_aligned.bam"
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task bam_stats_qc {
  File bamHeader
  File readCount
  File bam
  String bamName
  String outputDir

  command {
   bam_stats -i ${bam} -o ${outputDir + "/" + bamName + ".bas"} \
   && \
   verify_read_groups.pl --header-file ${bamHeader} --bas-file ${outputDir + "/" + bamName + ".bas"} --input-read-count-file ${readCount}
  }

  output {
    File bam_stats = "${outputDir}/${bamName}.bas"
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task merge {
  Array[File]+ inputBams
  String outputFilePrefix
  String outputDir
  Int threads

  command {
    bammarkduplicates \
    I=${sep=" I=" inputBams} \
    O=${outputDir + "/" + outputFilePrefix + ".bam"} \
    M=${outputDir + "/" + outputFilePrefix + ".metrics"} \
    tmpfile=${outputDir + "/" + outputFilePrefix + ".biormdup"} \
    markthreads=${threads} \
    rewritebam=1 \
    rewritebamlevel=1 \
    index=1 \
    md5=1
  }

  output {
    File merged_bam = "${outputDir}/${outputFilePrefix}.bam"
    File merged_bam_metrics = "${outputDir}/${outputFilePrefix}.metrics"
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task extract_unaligned_reads {
  File inputBam
  File reference_gz
  String outputFilePrefix
  String outputDir
  Int sortMemMb
  Int f

  command {
    samtools view -h -f ${f} ${inputBam} | \
    remove_both_ends_unmapped_reads.pl | \
    bamsort blockmb=${sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${outputFilePrefix + ".sorttmp"} O=${outputDir + "/" + outputFilePrefix + "_unmappedReads" + f + ".bam"}
  }

  output {
    File unmapped_reads = "${outputDir}/${outputFilePrefix}_unmappedReads_f${f}.bam"
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task extract_reads_both_mates_unaligned {
  File inputBam
  String outputFilePrefix
  String outputDir

  command {
    samtools view -h -b -f 12 ${inputBam} > "${outputDir}/${outputFilePrefix}_unmappedReads_f12.bam"
  }

  output {
    File unmapped_reads = "${outputDir}/${outputFilePrefix}_unmappedReads_f12.bam"
  }

  runtime {
    docker: "bwa-workflow"
  }
}

workflow bwa_workflow {
  Array[File]+ unalignedBams
  File reference_gz
  String outputFilePrefix
  String outputDir = "."
  Int sortMemMb
  Int threads

  scatter(bam in unalignedBams) {
    call get_basename {
      input: f=bam
   }

    call read_header {
      input: unalignedBam=bam, 
             bamName=get_basename.base,
             outputDir=outputDir
    }

    call count_reads {
      input: unalignedBam=bam,
             bamName=get_basename.base,
             outputDir=outputDir
    }

    call align {
      input: unalignedBam=bam,
             bamHeader=read_header.header,
             bamName=get_basename.base,
             threads=threads,
             sortMemMb=sortMemMb,
             outputDir=outputDir, 
             reference_gz=reference_gz
    }

    call bam_stats_qc {
      input: bam=align.bam_output,
             bamHeader=read_header.header_file,
             readCount=count_reads.counts_file,
             bamName=get_basename.base,
             outputDir=outputDir
    }
  }

  call merge {
    input: inputBams=align.bam_output,
           threads=threads,
           outputFilePrefix=outputFilePrefix, 
           outputDir=outputDir
  }

  call extract_unaligned_reads as reads_unmapped {
    input: inputBam=merge.merged_bam,
           f=4,
           sortMemMb=sortMemMb,
           outputFilePrefix=outputFilePrefix,
           outputDir=outputDir,
           reference_gz=reference_gz
  }

  call extract_unaligned_reads as reads_mate_unmapped {
    input: inputBam=merge.merged_bam,
           f=8,
           sortMemMb=sortMemMb,
           outputFilePrefix=outputFilePrefix,
           outputDir=outputDir,
           reference_gz=reference_gz
  }

  call extract_reads_both_mates_unaligned {
    input: inputBam=merge.merged_bam, 
           outputFilePrefix=outputFilePrefix, 
           outputDir=outputDir
  }

  call merge as merge_unmapped {
    input: inputBams=[reads_unmapped.unmapped_reads, reads_mate_unmapped.unmapped_reads, extract_reads_both_mates_unaligned.unmapped_reads],
           threads=threads,
           outputFilePrefix=outputFilePrefix,
           outputDir=outputDir
  }
}
