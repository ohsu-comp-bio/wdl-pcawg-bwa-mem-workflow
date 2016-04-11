task prepare_reference_files {
  File reference_gz
  File reference_gz_fai
  File reference_gz_amb
  File reference_gz_ann
  File reference_gz_bwt
  File reference_gz_pac
  File reference_gz_sa
  String referenceDir

  command {
    ln -s ${reference_gz} ${referenceDir + "/" + "genome.fa.gz"}; \
    ln -s ${reference_gz_fai} ${referenceDir + "/" + "genome.fa.gz.fai"}; \
    ln -s ${reference_gz_amb} ${referenceDir + "/" + "genome.fa.gz.64.amb"}; \
    ln -s ${reference_gz_ann} ${referenceDir + "/" + "genome.fa.gz.64.ann"}; \
    ln -s ${reference_gz_bwt} ${referenceDir + "/" + "genome.fa.gz.64.bwt"}; \
    ln -s ${reference_gz_pac} ${referenceDir + "/" + "genome.fa.gz.64.pac"}; \
    ln -s ${reference_gz_sa} ${referenceDir + "/" + "genome.fa.gz.64.sa"}
  }

  output {
    String reference_genome = "${referencDir}/genome.fa.gz"
  }
}

task get_basename {
  File f

  command {
    basename ${f}
  }

  output {
    String base = read_string(stdout())
  }
}

task read_header {
  File unalignedBam
  String bamName
  String outputDir

  command {
    samtools view -H ${unalignedBam} \
    | \
    perl -nae 'next unless /^\\@RG/; s/\\tPI:\\s*\\t/\\t/; s/\\tPI:\\s*\\z/\\n/; s/\\t/\\\\t/g; print' > "${outputDir}/${bamName}_header.txt"
  }

  output {
    File header = "${outputDir}/${bamName}_header.txt"
  }
}

task count_reads {
  File unalignedBam
  String bamName
  String outputDir

  command {
    samtools view ${unalignedBam} | wc -l > "${outputDir}/${bamName}_read_count.txt"
  }

  output {
    File counts_file = "${outputDir}/${bamName}_read_count.txt"
  }
}

task align {
  File unalignedBam
  File bamHeader
  File reference_gz
  String bamName
  String outputDir
  Int threads

  command {
    bamtofastq exlcude=QCFAIL,SECONDARY,SUPPLEMENTARY T=${bamName + ".t"} S=${bamName + ".s"} O=${bamName + ".o"} O2=${bamName + ".o2"} collate=1 tryoq=1 filename=${unalignedBam} \
    | \
    bwa mem -p -t ${threads} -T 0 -R ${sep="" read_lines(bamHeader)} ${reference_gz} - \
    | \
    bamsort inputformat=sam level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${bamName + ".sorttmp"} O=${outputDir + "/" + bamName + ".bam"} 2> ${outputDir + "/" + bamName + "_bamsort_info.txt"}
  }

  output {
    File bam_output = "${outputDir}/${bamName}_aligned.bam"
    File bam_sortinfo = "${outputDir}/${bamName}_bamsort_info.txt"
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
}

task extract_unaligned_reads {
  File inputBam
  String reference_gz
  String outputFilePrefix
  String outputDir
  Int f

  command {
    samtools view -h -f ${f} ${inputBam} \
    | \
    remove_both_ends_unmapped_reads.pl \
    | \
    bamsort inputformat=sam level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${outputFilePrefix + ".sorttmp"} O=${outputDir + "/" + outputFilePrefix + "_unmappedReads" + f + ".bam"}
  }

  output {
    File unmapped_reads = "${outputDir}/${outputFilePrefix}_unmappedReads_f${f}.bam"
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
}


workflow bwa_workflow {
  Array[File]+ unalignedBams
  String outputFilePrefix
  String ouputDir = "/output/"
  Int threads = 8

  call prepare_reference_files

  scatter(bam in unalignedBams) {
    call get_basename {
      input: f=bam
   }

    call read_header {
      input: unalignedBam=bam, 
             bamName=get_basename.base,
             outputDir=outputDir
    }

    call align {
      input: unalignedBam=bam,
             bamHeader=read_header.header,
             bamName=get_basename.base,
             threads=threads,
             outputDir=outputDir, 
             reference_gz=prepare_reference_files.reference_genome
    }

    call count_reads {
      input: unalignedBam=bam,
             bamName=get_basename.base,
             outputDir=outputDir
    }

    call bam_stats_qc {
      input: bam=align.bam_output,
             bamHeader=read_header.header,
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
           outputFilePrefix=outputFilePrefix, 
           outputDir=outputDir,
           reference_gz=prepare_reference_files.reference_genome
  }

  call extract_unaligned_reads as reads_mate_unmapped {
    input: inputBam=merge.merged_bam, 
           f=8, 
           outputFilePrefix=outputFilePrefix, 
           outputDir=outputDir,
           reference_gz=prepare_reference_files.reference_genome
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
