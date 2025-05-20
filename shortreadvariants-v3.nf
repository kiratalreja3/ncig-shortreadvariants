include { fromQuery } from 'plugin/nf-sqldb'

refmode = ["chm13"]

references = ["chm13.XX": params.chm13xx, "chm13.XY": params.chm13xy, "chm13.XX.interval" : params.chm13xxbed, "chm13.XY.interval" : params.chm13xybed, "chm13.XY.par" : params.chm13xypar]


fastqch = Channel
    .fromQuery("SELECT SGID, uniqueID, tubeid, sex, filename1, filename2 FROM illumina where coverage='High_Coverage'", db: 'inputdb')
    .map { row ->
        def ( participantID, uniqueID, tubeid, sex, filename1, filename2) = row
        meta = [RG:tubeid]
        return [participantID, meta, uniqueID,sex, filename1, filename2]
    }



process alignment {

    errorStrategy 'ignore'

    maxForks 50

    input:
    tuple val(sample), val(meta), val (uniqueid), val (sex), val(r1), val(r2)
    each mode

    output:
    tuple val(sample), val(meta), val (mode), val(sex), path("*cram")
    tuple val(sample), path("*cram.crai")

    script:
    def ref = references["${mode}.${sex}"]

    """

    bwa-mem2 mem -t \${PBS_NCPUS} -Y -K 100000000 -R "@RG\\tID:${meta.RG}\\tSM:${sample}\\tPU:Illumina" \
     ${ref} ${r1} ${r2} | samtools sort - --reference ${ref} -T \${PBS_JOBFS} -@ "${task.cpus}" --write-index --output-fmt CRAM \
     -o "${sample}.${meta.RG}.${uniqueid}.${mode}.cram"

    """

    stub:

    """

    touch "${sample}.${meta.RG}.${uniqueid}.${mode}.cram"
    touch "${sample}.${meta.RG}.${uniqueid}.${mode}.cram.crai"

    """

}

process mergecram {

    input:
    tuple val(sample), val(mode), val(sex), val(cramlist)

    output:
    tuple val(sample), val(mode), val(sex), path("*.cram")
    tuple val(sample), path("*cram.crai")

    script:

    def cramfiles = cramlist.join(' ')
    def ref = references["${mode}.${sex}"]

    """

    samtools merge -f --write-index --output-fmt CRAM --threads "${task.cpus}" --reference ${ref} \
     -o "${sample}.${mode}.merged.cram" ${cramfiles}

    """

    stub:

    """

    touch "${sample}.${mode}.merged.cram"
    touch "${sample}.${mode}.merged.cram.crai"
    
    """


}



process publishalignment {

    publishDir "${params.outdir}/${params.cohortID}/participantvariants/${sample}" , mode: 'copy', pattern:"*cram*"

    input:
    tuple val(sample), val (mode), val (sex), val(cram)

    output:
    tuple path ("*cram*"), path ("*.done")

    script:
    def refpath = references["${mode}.${sex}"]

    """
    
    mkdir -p "${params.outdir}/${params.cohortID}/participantvariants/${sample}"
    mkdir -p "${params.outdir}/${params.cohortID}/sql"

    ln ${cram} \${PWD}
    ln "${cram}.crai" \${PWD}

    publishedcram="${params.outdir}/${params.cohortID}/participantvariants/${sample}/\$(basename ${cram})"
    query="INSERT OR REPLACE INTO participantvariants (participantID, cohortID,ref,refpath,alignment) VALUES ('${sample}', '${params.cohortID}', '${mode}', '${refpath}','\${publishedcram}');"

    echo "\${query}" > ${params.outdir}/${params.cohortID}/sql/${sample}.${mode}.${sex}.alignment.sql

    #insertsql.sh -db ${params.outdb} -query "INSERT OR REPLACE INTO participantvariants (participantID, cohortID,ref,refpath,alignment) VALUES ('${sample}', '${params.cohortID}', '${mode}', '${refpath}','\${publishedcram}');"

    touch insertalignment.done

    """

}

process deepvariantdryrun {

    //capture the cram file and crai here and publish it per sample 

    input:
    tuple val(sample), val (mode), val (sex), val(cram)

    output:
    tuple val (sample), val (mode), val (sex), val (cram), env(make_examples_args), env(call_variants_args), env(post_process_args)

    script:

    def ref = references["${mode}.${sex}"]
    def parbed = references["chm13.XY.par"]
    def parbedparameter = (sex == "XX") ? "" : "--par_regions_bed ${parbed}"
    def haploidcontigs = (sex == "XX") ? "" : "--haploid_contigs chrX,chrY"
    def regionsbed = references["${mode}.${sex}.interval"]

    """
     
    run_deepvariant \
    --model_type WGS \
    --output_vcf "${sample}.${mode}.vcf" \
    --output_gvcf "${sample}.${mode}.g.vcf" \
    --ref ${ref} \
    --reads ${cram} \
    ${parbedparameter} ${haploidcontigs} \
    --regions ${regionsbed} \
    --dry_run=true  > commands.txt

    make_examples_args=\$(grep "/opt/deepvariant/bin/make_examples" commands.txt | awk -F'/opt/deepvariant/bin/make_examples' '{print \$2}' | sed 's/--mode calling//g' | sed 's/--ref "[^"]*"//g' | sed 's/--reads "[^"]*"//g' | sed 's/--sample_name "[^"]*"//g' | sed 's/--examples "[^"]*"//g' | sed 's/--gvcf "[^"]*"//g')
    call_variants_args=\$(grep "/opt/deepvariant/bin/call_variants" commands.txt | awk -F'/opt/deepvariant/bin/call_variants' '{print \$2}' | sed 's/--outfile "[^"]*"//g' | sed 's/--examples "[^"]*"//g')
    post_process_args="${haploidcontigs} --par_regions_bed ${parbed}"


    """

    stub:

    """

    make_examples_args="test"
    call_variants_args="test"
    post_process_args="test"

    """
}


process deepvariantmakeexamples {

    input:
    tuple val (sample), val (mode), val (sex), val (cram), val(make_examples_args), val(call_variants_args), val (post_process_args)

    output:
    tuple val (sample),val (mode), val (sex), path('make_examples*.gz{,.example_info.json}'), path('gvcf*.gz{,.example_info.json}'), val(call_variants_args), val(post_process_args)

    script:
    def ref = references["${mode}.${sex}"]


    """

    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer make_examples \
    --mode calling --ref "${ref}" --reads "${cram}" --sample_name "${sample}" --examples "make_examples.tfrecord@${task.cpus}.gz" \
    --gvcf "gvcf.tfrecord@${task.cpus}.gz"  ${make_examples_args}

    """

    stub:

    """

    touch "make_examples.tfrecord@20.gz"
    touch "gvcf.tfrecord@20.gz"
    

    """



}

process deepvariantcallvariants {

    input:
    tuple val(sample), val (mode), val (sex) ,path(make_examples_out), path(gvcf_make_examples_out), val(call_variants_args), val (post_process_args)

    output:
    tuple val(sample), val (mode), val (sex) , path(gvcf_make_examples_out), path('*.gz'), val (post_process_args)

    script:

    def ref = references["${mode}.${sex}"]
    def matcher = make_examples_out[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
    def num_shards = matcher[0][2] as int

    """
    call_variants --outfile "call_variants_output.tfrecord.gz" --examples "make_examples.tfrecord@${num_shards}.gz" ${call_variants_args}
    """

    stub:

    """
    touch call_variants_output.tfrecord.gz
    """

}


process deepvariantpostprocess {

    publishDir "${params.outdir}/${params.cohortID}/participantvariants/${sample}" , mode: 'copy', pattern:"*html*"

    input:
    tuple val(sample), val (mode), val (sex), path(gvcf_make_examples_out), path(call_variants_out), val (post_process_args)

    output:
    tuple val(sample), val (mode), path("${sample}.${mode}.g.vcf.gz")
    tuple val(sample), val (mode), path("${sample}.${mode}.vcf.gz")
    tuple val(sample), val (mode), path("${sample}.${mode}.vcf.gz.tbi")
    tuple val(sample), val (mode), path("${sample}.${mode}.visual_report.html")

    script:
    def matcher = gvcf_make_examples_out[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
    def num_shards = matcher[0][2] as int
    def ref = references["${mode}.${sex}"]

    """
    postprocess_variants --ref "${ref}" --infile "call_variants_output.tfrecord.gz" --outfile "${sample}.${mode}.vcf.gz" \
     --cpus "${task.cpus}" --gvcf_outfile "${sample}.${mode}.g.vcf.gz" --sample_name "${sample}" --nonvariant_site_tfrecord_path "gvcf.tfrecord@${num_shards}.gz" ${post_process_args} > prog.out.\${PBS_JOBID} 2> prog.err.\${PBS_JOBID}

    vcf_stats_report --input_vcf "${sample}.${mode}.vcf.gz" --outfile_base "${sample}.${mode}"
    """

    stub:
    """
    touch "${sample}.${mode}.g.vcf.gz"
    touch "${sample}.${mode}.vcf.gz"
    
    """

}



process publishvariants {

    publishDir "${params.outdir}/${params.cohortID}/participantvariants/${sample}" , mode: 'copy', pattern:"*vcf*"

    input:
    tuple val(sample), val (mode), val (gvcf)
    tuple val(sample), val (mode), val (vcf)
    tuple val(cram), val (done)

    output:
    tuple path ("*g.vcf*"), path ("*.done")
    tuple path ("*vcf.gz*"), path ("*.done")

    script:

    """

    ln ${gvcf}  \${PWD}
    ln "${gvcf}.tbi" \${PWD}

    ln ${vcf}  \${PWD}
    ln "${vcf}.tbi"  \${PWD}

    publishedgvcf="${params.outdir}/${params.cohortID}/participantvariants/${sample}/\$(basename ${gvcf})"
    publishedvcf="${params.outdir}/${params.cohortID}/participantvariants/${sample}/\$(basename ${vcf})"

    query="UPDATE participantvariants SET gvcf = '\${publishedgvcf}', vcf = '\${publishedvcf}' WHERE participantID = '${sample}' AND cohortID = '${params.cohortID}' AND ref = '${mode}'"
    echo "\${query}" > ${params.outdir}/${params.cohortID}/sql/${sample}.${mode}.variants.sql
    
    #insertsql.sh -db ${params.outdb} -query "UPDATE participantvariants SET gvcf = '\${publishedgvcf}', vcf = '\${publishedvcf}' WHERE participantID = '${sample}' AND cohortID = '${params.cohortID}' AND ref = '${mode}'"
    
    touch insertvariants.done

    """

}




process splitgvcf {

    input:
    tuple val(sample), val (mode), val (gvcf)

    output:
    tuple val(sample), val (mode), path("*.g.vcf.gz")

    script:

    """

    bcftools view -h ${gvcf} | grep "^##contig"  | awk -F',' '{print \$1}' | cut -f3 -d "=" | grep -v "chrM" > \${PBS_JOBFS}/chrs.txt
    cat \${PBS_JOBFS}/chrs.txt | parallel -j "\${PBS_NCPUS}" split-reheader-gvcf.sh -input ${gvcf} -output ${sample}.${mode}-{}.g.vcf.gz -chr {}

    """

    stub:

    """

    for i in {1..3} X Y; do
        touch "${sample}.${mode}-chr\${i}.g.vcf.gz"
    done

    """

}


process jointgenotyping {

    errorStrategy 'ignore'

    input:
    tuple val (mode), val(chr), val(gvcfPaths)

    output:
    tuple val (mode), val(chr), path("${chr}.${mode}.vcf.gz")
    tuple val (mode), val(chr), path("${chr}.${mode}.vcf.gz.tbi")

    script:

    """
    cat > \${PBS_JOBFS}/tmp.${chr}_paths.txt <<EOF
    ${gvcfPaths}
    EOF

    cat \${PBS_JOBFS}/tmp.${chr}_paths.txt | tr -d '[],' | tr ' ' '\\n' > \${PBS_JOBFS}/${chr}_paths.txt

    if [ ${chr} = "chrY" ]; then

        while read file; do

            num_records=\$(bcftools stats "\$file" | grep "number of records:" | awk '{print \$NF}')
            if [ "\$num_records" -gt 0 ]; then
                echo "\$file" >> "\${PBS_JOBFS}/${chr}_paths.chrY.txt"
            fi

        done < \${PBS_JOBFS}/${chr}_paths.txt

        rm \${PBS_JOBFS}/${chr}_paths.txt
        mv \${PBS_JOBFS}/${chr}_paths.chrY.txt \${PBS_JOBFS}/${chr}_paths.txt

    fi

    glnexus_cli --list \${PBS_JOBFS}/${chr}_paths.txt \
     --threads "${task.cpus}" --config DeepVariant > \${PBS_JOBFS}/${chr}.bcf

    bcftools view \${PBS_JOBFS}/${chr}.bcf  | bgzip -@ "${task.cpus}" -c > ${chr}.${mode}.vcf.gz

    tabix -p vcf ${chr}.${mode}.vcf.gz

    """
    
    stub:

    """

    touch "${chr}.${mode}.g.vcf.gz"
    touch "${chr}.${mode}.g.vcf.gz.tbi"
  
    """

}

process publishjointgenotyping {

    publishDir "${params.outdir}/${params.cohortID}/cohortvariants" , mode: 'copy', pattern:"*vcf*"

    input:
    tuple val(mode), val (chr), val (vcf)

    output:
    tuple path ("*vcf*"), path ("*.done")

    script:

    """
    
    mkdir -p "${params.outdir}/${params.cohortID}/cohortvariants"

    ln ${vcf}  \${PWD}
    ln "${vcf}.tbi" \${PWD}

    publishedvcf="${params.outdir}/${params.cohortID}/cohortvariants/\$(basename ${vcf})"

    query="INSERT OR REPLACE INTO cohortvariants (cohortID,ref,chr,vcf) VALUES ('${params.cohortID}', '${mode}', '${chr}','\${publishedvcf}');"
    echo "\${query}" > ${params.outdir}/${params.cohortID}/sql/${chr}.${mode}.jointvariants.sql

    #insertsql.sh -db ${params.outdb} -query "INSERT OR REPLACE INTO cohortvariants (cohortID,ref,chr,vcf) VALUES ('${params.cohortID}', '${mode}', '${chr}','\${publishedvcf}');"

    touch insertcohortvariants.done

    """

}



workflow {

    alignmentCh = alignment(fastqch,refmode)
    
    //subsetting alignmentCh to get only cram files
    //and grouping the files by sample, sex and reference
    //to identify which of them need to be merged
    cramCh = alignmentCh[0].groupTuple(by: [0, 1, 2]).map
            { tuple -> tuple[0,2,3,4]}.map 
            { tuple -> [tuple[0], tuple[1], tuple[2].unique()[0], tuple[3]]}

    //splitting the cramCh into two channels
    //based on what needs to be merged and what
    //can be immediately processed 

    toMergeCramCh = cramCh.filter { tuple ->
        def (sample, sex, mode, cramFiles) = tuple
        cramFiles.size() > 1
    }

    noMergeCramCh = cramCh.filter { tuple ->
        def (sample, sex, mode, cramFiles) = tuple
        cramFiles.size() <= 1 
    }
    
    //getting rid of duplicate entires that come 
    //from using groupTuple
    noMergeCramCh = noMergeCramCh.map 
    { tuple -> [tuple[0], tuple[1], tuple[2], tuple[3][0]]}
    
    //Merging the relevant samples
    mergedCramCh = mergecram(toMergeCramCh)
    mergedCramCh = mergedCramCh[0]

    //Merging the merged and non-merged channels
    cramCh = mergedCramCh.mix(noMergeCramCh)

    //Calling variants
    dryrunCh = deepvariantdryrun(cramCh)
    makeexamplesCh =deepvariantmakeexamples(dryrunCh)
    callvariantsCh = deepvariantcallvariants(makeexamplesCh)
    postprocessCh = deepvariantpostprocess(callvariantsCh)
    
    splitgcvfCh = splitgvcf(postprocessCh[0])
    splitgcvfCh = splitgcvfCh.map { tuple ->
            tuple[1,2] 
            }.groupTuple().map { tuple ->
                [tuple[0], tuple[1].flatten()] 
            }.flatMap { ref, paths ->
                    paths.collect { path ->
                    def match = path =~ /.*-(chr\w+)\.g\.vcf\.gz/
                    def chromosome = match ? match[0][1] : "unknown"
                    [ref, chromosome, path]
                }}.groupTuple(by: [0,1])


    genotypeCh = jointgenotyping(splitgcvfCh)
    publishalignmentCh = publishalignment(cramCh)
    publishvariants(postprocessCh[0],postprocessCh[1],publishalignmentCh)
    publishjointgenotyping(genotypeCh[0])

}


