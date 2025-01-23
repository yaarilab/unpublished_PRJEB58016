$HOSTNAME = ""
params.outdir = 'results'  

evaluate(new File("${params.projectDir}/nextflow_header.config"))
params.metadata.metadata = "${params.projectDir}/tools.json"


if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 

Channel.value(params.mate).into{g_1_mate_g_63;g_1_mate_g_78;g_1_mate_g_71;g_1_mate_g_85;g_1_mate_g_79;g_1_mate_g38_11;g_1_mate_g38_9;g_1_mate_g38_12;g_1_mate_g80_14;g_1_mate_g82_9;g_1_mate_g28_15;g_1_mate_g28_19;g_1_mate_g28_12;g_1_mate_g73_15;g_1_mate_g73_19;g_1_mate_g73_12}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_64_reads_g_63}
 } else {  
	g_64_reads_g_63 = Channel.empty()
 }



process unizp {

input:
 set val(name),file(reads) from g_64_reads_g_63
 val mate from g_1_mate_g_63

output:
 set val(name),file("*.fastq")  into g_63_reads0_g38_11

script:

if(mate=="pair"){
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	
	"""
	case "$R1" in
	*.gz | *.tgz ) 
	        gunzip -c $R1 > R1.fastq
	        ;;
	*)
	        cp $R1 ./R1.fastq
	        echo "$R1 not gzipped"
	        ;;
	esac
	
	case "$R2" in
	*.gz | *.tgz ) 
	        gunzip -c $R2 > R2.fastq
	        ;;
	*)
	        cp $R2 ./R2.fastq
	        echo "$R2 not gzipped"
	        ;;
	esac
	"""
}else{
	"""
	case "$reads" in
	*.gz | *.tgz ) 
	        gunzip -c $reads > R1.fastq
	        ;;
	*)
	        cp $reads ./R1.fastq
	        echo "$reads not gzipped"
	        ;;
	esac
	"""
}
}


process Mask_Primer_MaskPrimers {

input:
 val mate from g_1_mate_g38_11
 set val(name),file(reads) from g_63_reads0_g38_11

output:
 set val(name), file("*_primers-pass.fast*") optional true  into g38_11_reads0_g_71
 set val(name), file("*_primers-fail.fast*") optional true  into g38_11_reads_failed11
 set val(name), file("MP_*")  into g38_11_logFile2_g38_9
 set val(name),file("out*")  into g38_11_logFile3_g61_0

script:
method = params.Mask_Primer_MaskPrimers.method
barcode_field = params.Mask_Primer_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_MaskPrimers.primer_field
barcode = params.Mask_Primer_MaskPrimers.barcode
revpr = params.Mask_Primer_MaskPrimers.revpr
mode = params.Mask_Primer_MaskPrimers.mode
failed = params.Mask_Primer_MaskPrimers.failed
fasta = params.Mask_Primer_MaskPrimers.fasta
nproc = params.Mask_Primer_MaskPrimers.nproc
maxerror = params.Mask_Primer_MaskPrimers.maxerror
umi_length = params.Mask_Primer_MaskPrimers.umi_length
start = params.Mask_Primer_MaskPrimers.start
extract_length = params.Mask_Primer_MaskPrimers.extract_length
maxlen = params.Mask_Primer_MaskPrimers.maxlen
skiprc = params.Mask_Primer_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""
fasta = (fasta=="true") ? "--fasta" : ""
def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk} ${pf} ${bf}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log & \
	wait
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${R1}_MP.log
	"""
}

}


process Cluster_UMI_p11 {

input:
 set val(name),file(reads) from g38_11_reads0_g_71
 val mate from g_1_mate_g_71

output:
 set val(name),file("*_umi-pass.fastq") optional true  into g_71_reads0_g80_14

script:
umi_field = params.Cluster_UMI_p11.umi_field
umi_cluster_script = params.Cluster_UMI_p11.umi_cluster_script

readArray = reads.toString().split(' ')	

R1 = readArray[0]
R2 = readArray[1]



filename_R1=file(R1).getSimpleName()
filename_R2=file(R2).getSimpleName()



"""
#!/bin/bash

# pair awk between R1 and R2, so isotype information will be transfered
awk 'NR%4==1 {split(\$0,a,"|"); split (a[4],isoType,"="); print(a[1]"|"a[2]"|"a[3]"|ISOTYPE="isoType[2]);} NR%4!=1 {print \$0}' ${R2} > tmp_R2.fastq
awk 'NR==FNR && NR%4==1 {split(\$0,a,"|"); split(a[1],b,"/"); id=b[1]; z[id]=1 ; split(a[4],c,"ISOTYPE="); bc[id]=c[2];} NR!=FNR && FNR%4==1 {split(\$0,a,"|"); split(a[1],b,"/"); id=b[1]; if (z[id]==1) {header=a[1]"|"a[2]"|"a[3]"|ISOTYPE="bc[id]; print (header);}} NR!=FNR && FNR%4!=1 {if (z[id]==1) print}' tmp_R2.fastq ${R1} > tmp_R1.fastq


cat tmp_R1.fastq tmp_R2.fastq | awk '{if (NR%4==1) {split(\$0,a,"${umi_field}="); print a[2]}}' > pair_pass.umi

python3 ${umi_cluster_script} -i pair_pass.umi -o pair_pass.umi.convert
awk -F'\t' 'NR==FNR && NR>1 {umi[\$1]=\$2;} NR!=FNR {if(/UMI=/){split(\$0,a,"${umi_field}=");split(a[2],b,"|");sub(/[/]/, "", b[1]);sub(/[/12]/, "", b[1]); print a[1]"UMI="b[1]"|"b[2] ;} else {print}}' pair_pass.umi.convert tmp_R1.fastq > 'tmp_R1_umi-pass.fastq'
awk -F'\t' 'NR==FNR && NR>1 {umi[\$1]=\$2;} NR!=FNR {if(/UMI=/){split(\$0,a,"${umi_field}=");split(a[2],b,"|");sub(/[/]/, "", b[1]);sub(/[/12]/, "", b[1]); print a[1]"UMI="b[1]"|"b[2] ;} else {print}}' pair_pass.umi.convert tmp_R2.fastq > 'tmp_R2_umi-pass.fastq'

#awk -F'\t' 'NR==FNR && NR>1 {umi[\$1]=\$2;} NR!=FNR {if(/UMI=/){split(\$0,a,"${umi_field}=");split(a[2],b,"|"); print a[1]"UMI="b[1];} else {print}}' pair_pass.umi.convert ${R1} > ${R1}'_umi-pass.fastq'
#awk -F'\t' 'NR==FNR && NR>1 {umi[\$1]=\$2;} NR!=FNR {if(/UMI=/){split(\$0,a,"${umi_field}=");split(a[2],b,"|"); print a[1]"UMI="b[1]"|"b[2] ;} else {print}}' pair_pass.umi.convert ${R2} > ${R2}'_umi-pass.fastq'


"""


}


process Cluster_Sets_cluster_sets {

input:
 set val(name),file(reads) from g_71_reads0_g80_14
 val mate from g_1_mate_g80_14

output:
 set val(name),file("*_cluster-pass.fastq")  into g80_14_reads0_g_81
 set val(name),file("*_cluster-fail.fastq") optional true  into g80_14_reads_failed11

script:
method = params.Cluster_Sets_cluster_sets.method
failed = params.Cluster_Sets_cluster_sets.failed
nproc = params.Cluster_Sets_cluster_sets.nproc
cluster_field = params.Cluster_Sets_cluster_sets.cluster_field
ident = params.Cluster_Sets_cluster_sets.ident
length = params.Cluster_Sets_cluster_sets.length
prefix = params.Cluster_Sets_cluster_sets.prefix
cluster_tool = params.Cluster_Sets_cluster_sets.cluster_tool
cluster_exec = params.Cluster_Sets_cluster_sets.cluster_exec
usearch_version = params.Cluster_Sets_cluster_sets.usearch_version
set_field = params.Cluster_Sets_cluster_sets.set_field
start = params.Cluster_Sets_cluster_sets.start
end = params.Cluster_Sets_cluster_sets.end
barcode_field = params.Cluster_Sets_cluster_sets.barcode_field
//* @style @condition:{method="set",set_field,start,end},{method="all",start,end},{method="barcode",barcode_field} @array:{method,failed,cluster_field,ident,length,prefix,cluster_tool,cluster_exec,set_field,start,end,barcode_field}  @multicolumn:{method,failed,nproc,cluster_field,ident,length,prefix,cluster_tool,cluster_exec},{set_field,start,end,barcode_field}

method = (method.size==2) ? method : [method[0],method[0]]
failed = (failed.size==2) ? failed : [failed[0],failed[0]]
cluster_field = (cluster_field.size==2) ? cluster_field : [cluster_field[0],cluster_field[0]]
ident = (ident.size==2) ? ident : [ident[0],ident[0]]
length = (length.size==2) ? length : [length[0],length[0]]
prefix = (prefix.size==2) ? prefix : [prefix[0],prefix[0]]
set_field = (set_field.size==2) ? set_field : [set_field[0],set_field[0]]
start = (start.size==2) ? start : [start[0],start[0]]
end = (end.size==2) ? end : [end[0],end[0]]
barcode_field = (barcode_field.size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]

def args_values = [];
[method, failed, cluster_field, ident, length, prefix, set_field, start, end, barcode_field].transpose().each { m, f, cf, i, l, p, sf, s, e, bf -> {
    f = (f=="true") ? "--failed" : ""
    p = (p=="") ? "" : "--prefix ${p}" 
    //ce = (ce=="") ? "" : "--exec ${ce}" 
    sf = (m=="set") ? "-f ${sf}" : ""
    s = (m=="barcode") ? "" : "--start ${s}" 
    e = (m=="barcode") ? "" : (e=="") ? "" : "--end ${e}" 
    bf = (m=="barcode") ? "-f ${bf}" : ""
    args_values.add("${m} ${f} -k ${cf} --ident ${i} --length ${l} ${p} ${sf} ${s} ${e} ${bf}")
}}


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	args_1 = args_values[0]
	args_2 = args_values[1]
	
	"""
	
	if  [ "${cluster_tool}" == "usearch"] && ["${cluster_exec}"==""]; then
		wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
		gunzip usearch${usearch_version}_i86linux32.gz
		chmod +x usearch${usearch_version}_i86linux32
		mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
		ce=" --exec /usr/local/bin/usearch2"
		echo "usearch"
	else
		echo "no usearch"
		if [ "${cluster_exec}" != "" ]; then
			ce=" --exec ${cluster_exec}"
			ct=" --cluster ${cluster_tool}"
		else
			ce=""
			ct=" --cluster ${cluster_tool}"
		fi
	fi
	
	ClusterSets.py ${args_1} -s $R1  --nproc ${nproc} \$ce 
	ClusterSets.py ${args_2} -s $R2  --nproc ${nproc} \$ce 
	"""
}else{
	args_1 = args_values[0]
	"""
	if  [ "${cluster_tool}" == "usearch" ] && ["${cluster_exec}"==""]; then
		wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
		gunzip usearch${usearch_version}_i86linux32.gz
		chmod +x usearch${usearch_version}_i86linux32
		mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
		ce=" --exec /usr/local/bin/usearch2"
		ct=" --cluster userach"
	else
		if [ "${cluster_exec}" != "" ]; then
			ce=" --exec ${cluster_exec}"
			ct=" --cluster ${cluster_tool}"
		else
			ce=""
			ct=" --cluster ${cluster_tool}"
		fi
	fi
	
	ClusterSets.py ${args_1} -s $reads --nproc ${nproc} \$ce \$ct
	"""
}


}


process take_the_biggest_cluster_of_each_UMI_p11 {

input:
 set val(name),file(reads) from g80_14_reads0_g_81

output:
 set val(name),file("*_select-pass.fastq")  into g_81_reads0_g_85

script:

readArray = reads.toString().split(' ')	

R1 = readArray[0]
R2 = readArray[1]

"""
#!/bin/bash

TMP_SUFFIX=".UMI_TO_CLUSTER.txt"

awk 'BEGIN{LARGEST_CLID="";LARGEST_SIZE=0;print "UMI\tLARGEST_CLID\tCLID_SIZE";} 
	 NR%4==1 {split(\$0,a,"|"); split(a[3],b,"UMI="); UMI=b[2]; split(a[5],c,"CLID="); CLID=c[2];  if(UMI==PREV_UMI || NR==1) {ARR[UMI"_"CLID]+=1; if(ARR[UMI"_"CLID]>LARGEST_SIZE) {LARGEST_SIZE=ARR[UMI"_"CLID]; LARGEST_CLID=CLID} ; if (NR==1) {PREV_UMI=UMI}; next } else { print PREV_UMI"\t"LARGEST_CLID"\t"LARGEST_SIZE; ARR[UMI"_"CLID]+=1; LARGEST_CLID=CLID;LARGEST_SIZE=ARR[UMI"_"CLID]; PREV_UMI=UMI; next }} 
	 END { print PREV_UMI"\t"LARGEST_CLID"\t"LARGEST_SIZE;} ' \
	 ${R1}  > R1\$TMP_SUFFIX
				
awk 'NR==FNR && NR>1 {UMI=\$1;CLID=\$2;ID[UMI]=CLID} 
	 NR!=FNR && FNR%4==1 {split(\$0,a,"|"); split(a[3],b,"UMI="); UMI=b[2]; split(a[5],c,"CLID="); CLID=c[2];  if(CLID==ID[UMI]){flag=1; print}else{flag=0}}
	 NR!=FNR && FNR%4!=1 && flag==1 {print}' \
	 R1\$TMP_SUFFIX ${R1} > R1_select-pass.fastq
		 
# combine cluster and UMI fields
sed -i 's/UMI=\\([^|]*\\)|ISOTYPE=\\([^|]*\\)|CLID=\\([^|]*\\)/UMI=\\1-\\3|ISOTYPE=\\2/' R1_select-pass.fastq

## pair awk between M1S and Z, so cluster-UMI information will be transfered
awk 'NR==FNR && NR%4==1 {split(\$0,a,"|"); id=a[1]; split(a[3],b,"UMI="); umi=b[2]; z[id]=umi;} 
	 NR!=FNR && FNR%4==1 {split(\$0,a,"|"); id=a[1]; header=a[1]"|"a[2]"|"a[4]"|UMI="z[id]; print (header);} 
	 NR!=FNR && FNR%4!=1 {print}' \
	 R1_select-pass.fastq ${R2} > R2_select-pass.fastq

"""
}


process align_sets {

input:
 set val(name),file(reads) from g_81_reads0_g_85
 val mate from g_1_mate_g_85

output:
 set val(name),file("*_align-pass.fastq")  into g_85_reads0_g82_9
 set val(name), file("AS_*")  into g_85_logFile11
 set val(name),file("*_align-fail.fastq") optional true  into g_85_reads_failed22
 set val(name), file("out*") optional true  into g_85_logFile33

script:
method = params.align_sets.method
bf = params.align_sets.bf
div = params.align_sets.div
failed = params.align_sets.failed
nproc = params.align_sets.nproc
alignset_script = params.align_sets.alignset_script
alignset_script_batch_size = params.align_sets.alignset_script_batch_size
muscle_exec = params.align_sets.muscle_exec
muscle_version = params.align_sets.muscle_version
offset_table = params.align_sets.offset_table
pf = params.align_sets.pf
mode = params.align_sets.mode

primer_file = params.align_sets.primer_file
reverse = params.align_sets.reverse

//* @style @condition:{method="muscle",muscle_exec,muscle_version}, {method="offset",offset_table,pf,mode}, {method="table",muscle_exec,primer_file,reverse} @multicolumn:{method,bf,div,nproc},{offset,pf,mode}, {primer_file,reverse}


readArray = reads.toString().split(' ')	

reverse_arg = (reverse=="false") ? "" : "--reverse"
div_arg = (div=="false") ? "" : "--div"
failed_arg = (failed=="true") ? "--failed" : "" 
bf = (bf=="") ? "" : "--bf ${bf}"

primer_file_argv = ""

if(method=="offset"){
	pf = "--pf ${pf}"
	mode = "--mode ${mode}"
	offset_table_argv = "-d ${offset_table}"
	muscle_exec_argv = ""
}else{
	pf = ""
	mode = ""
	offset_table_argv = ""
	
	if(method=="table"){
		primer_file_argv = "-p ${primer_file}"
	}
}

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	
	"""
	if [ "${method}" == "muscle" ]; then
		if  [ "${muscle_version}" != "" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/muscle/downloads${muscle_version}/muscle${muscle_version}_i86linux64.tar.gz
			tar -xvzf muscle${muscle_version}_i86linux64.tar.gz
			chmod +x muscle${muscle_version}_i86linux64
			mv muscle${muscle_version}_i86linux64 /usr/local/bin/muscle2
			muscle_exec_argv="--exec /usr/local/bin/muscle2"
		else
			muscle_exec_argv="--exec ${muscle_exec}"
		fi
	else
		muscle_exec_argv=""
	fi
	
	AlignSets.py ${method} -s ${R1} ${bf} \$muscle_exec_argv ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R1_${name}.log --nproc ${nproc} | tee -a out_${R1}_AS.log
	#AlignSets.py ${method} -s ${R2} ${bf} \$muscle_exec_argv ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R2_${name}.log --nproc ${nproc} | tee -a out_${R1}_AS.log
	# adapt for P11
	chmod +x ${alignset_script}
	${alignset_script} ${R2} ${alignset_script_batch_size}
	"""
	
}else{
	R1 = readArray[0]
	"""
	if [ "${method}" == "muscle" ]; then
		if  [ "${muscle_version}" != "" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/muscle/downloads${muscle_version}/muscle${muscle_version}_i86linux64.tar.gz
			tar -xvzf muscle${muscle_version}_i86linux64.tar.gz
			chmod +x muscle${muscle_version}_i86linux64
			mv muscle${muscle_version}_i86linux64 /usr/local/bin/muscle2
			muscle_exec_argv="--exec /usr/local/bin/muscle2"
		else
			muscle_exec_argv="--exec ${muscle_exec}"
		fi
	else
		muscle_exec_argv=""
	fi
	
	AlignSets.py ${method} -s ${R1} ${bf} \$muscle_exec_argv ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R1_${name}.log --nproc ${nproc} | tee -a out_${R1}_AS.log
	"""
}

}


process Pair_Sequence_per_consensus_pair_seq {

input:
 set val(name),file(reads) from g_85_reads0_g82_9
 val mate from g_1_mate_g82_9

output:
 set val(name),file("*_pair-pass.fastq")  into g82_9_reads0_g_78
 set val(name),file("out*")  into g82_9_logFile11

script:
coord = params.Pair_Sequence_per_consensus_pair_seq.coord
act = params.Pair_Sequence_per_consensus_pair_seq.act
copy_fields_1 = params.Pair_Sequence_per_consensus_pair_seq.copy_fields_1
copy_fields_2 = params.Pair_Sequence_per_consensus_pair_seq.copy_fields_2
failed = params.Pair_Sequence_per_consensus_pair_seq.failed
nproc = params.Pair_Sequence_per_consensus_pair_seq.nproc

if(mate=="pair"){
	
	act = (act=="none") ? "" : "--act ${act}"
	failed = (failed=="true") ? "--failed" : "" 
	copy_fields_1 = (copy_fields_1=="") ? "" : "--1f ${copy_fields_1}" 
	copy_fields_2 = (copy_fields_2=="") ? "" : "--2f ${copy_fields_2}"
	
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	PairSeq.py -1 ${R1} -2 ${R2} ${copy_fields_1} ${copy_fields_2} --coord ${coord} ${act} ${failed} >> out_${R1}_PS.log
	"""
}else{
	
	"""
	echo -e 'PairSeq works only on pair-end reads.'
	"""
}


}

boolean isCollectionOrArray_bc(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep){
	def args_values;
    if(isCollectionOrArray_bc(barcode_field) || isCollectionOrArray_bc(primer_field) || isCollectionOrArray_bc(copy_field) || isCollectionOrArray_bc(mincount) || isCollectionOrArray_bc(minqual) || isCollectionOrArray_bc(minfreq) || isCollectionOrArray_bc(maxerror) || isCollectionOrArray_bc(prcons) || isCollectionOrArray_bc(maxgap) || isCollectionOrArray_bc(maxdiv) || isCollectionOrArray_bc(dep)){
    	primer_field = (isCollectionOrArray_bc(primer_field)) ? primer_field : [primer_field,primer_field]
    	act = (isCollectionOrArray_bc(act)) ? act : [act,act]
    	copy_field = (isCollectionOrArray_bc(copy_field)) ? copy_field : [copy_field,copy_field]
    	mincount = (isCollectionOrArray_bc(mincount)) ? mincount : [mincount,mincount]
    	minqual = (isCollectionOrArray_bc(minqual)) ? minqual : [minqual,minqual]
    	minfreq = (isCollectionOrArray_bc(minfreq)) ? minfreq : [minfreq,minfreq]
    	maxerror = (isCollectionOrArray_bc(maxerror)) ? maxerror : [maxerror,maxerror]
    	prcons = (isCollectionOrArray_bc(prcons)) ? prcons : [prcons,prcons]
    	maxgap = (isCollectionOrArray_bc(maxgap)) ? maxgap : [maxgap,maxgap]
    	maxdiv = (isCollectionOrArray_bc(maxdiv)) ? maxdiv : [maxdiv,maxdiv]
    	dep = (isCollectionOrArray_bc(dep)) ? dep : [dep,dep]
    	args_values = []
        [barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep].transpose().each { bf,pf,a,cf,mc,mq,mf,mr,pc,mg,md,d -> {
            bf = (bf=="") ? "" : "--bf ${bf}"
            pf = (pf=="") ? "" : "--pf ${pf}" 
            a = (a=="none") ? "" : "--act ${a}" 
            cf = (cf=="") ? "" : "--cf ${cf}" 
            mr = (mr=="none") ? "" : "--maxerror ${mr}" 
            pc = (pc=="none") ? "" : "--prcons ${pc}" 
            mg = (mg=="none") ? "" : "--maxgap ${mg}" 
            md = (md=="none") ? "" : "--maxdiv ${md}" 
            mc = (mc=="none") ? "" : "--n ${mc}" 
            d = (d=="true") ? "--dep" : "" 
            args_values.add("${bf} ${pf} ${a} ${cf} ${mc} -q ${mq} --freq ${mf} ${mr} ${pc} ${mg} ${md} ${d}")
        }}
    }else{
        barcode_field = (barcode_field=="") ? "" : "--bf ${barcode_field}"
        primer_field = (primer_field=="") ? "" : "--pf ${primer_field}" 
        act = (act=="none") ? "" : "--act ${act}" 
        copy_field = (copy_field=="") ? "" : "--cf ${copy_field}" 
        maxerror = (maxerror=="none") ? "" : "--maxerror ${maxerror}" 
        prcons = (prcons=="none") ? "" : "--prcons ${prcons}" 
        maxgap = (maxgap=="none") ? "" : "--maxgap ${maxgap}" 
        maxdiv = (maxdiv=="none") ? "" : "--maxdiv ${maxdiv}" 
        dep = (dep=="true") ? "--dep" : "" 
        args_values = "${barcode_field} ${primer_field} ${act} ${copy_field} -n ${mincount} -q ${minqual} --freq ${minfreq} ${maxerror} ${prcons} ${maxgap} ${maxdiv} ${dep}"
    }
    return args_values
}


process build_consensus {

input:
 set val(name),file(reads) from g82_9_reads0_g_78
 val mate from g_1_mate_g_78

output:
 set val(name),file("*_consensus-pass.fastq")  into g_78_reads0_g_79
 set val(name),file("BC*")  into g_78_logFile11
 set val(name),file("out*")  into g_78_logFile22

script:
failed = params.build_consensus.failed
nproc = params.build_consensus.nproc
barcode_field = params.build_consensus.barcode_field
primer_field = params.build_consensus.primer_field
act = params.build_consensus.act
copy_field = params.build_consensus.copy_field
mincount = params.build_consensus.mincount
minqual = params.build_consensus.minqual
minfreq = params.build_consensus.minfreq
maxerror = params.build_consensus.maxerror
prcons = params.build_consensus.prcons
maxgap = params.build_consensus.maxgap
maxdiv = params.build_consensus.maxdiv
dep = params.build_consensus.dep
//* @style @condition:{act="none",},{act="min",copy_field},{act="max",copy_field},{act="sum",copy_field},{act="set",copy_field},{act="majority",copy_field} @array:{barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep} @multicolumn:{failed,nproc},{barcode_field,primer_field,act,copy_field}, {mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep}

args_values_bc = args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep)

// args 
if(isCollectionOrArray_bc(args_values_bc)){
	args_1 = args_values_bc[0]
	args_2 = args_values_bc[1]
}else{
	args_1 = args_values_bc
	args_2 = args_values_bc
}

failed = (failed=="true") ? "--failed" : "" 


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	
	"""
	BuildConsensus.py --version
	BuildConsensus.py -s $R1 ${args_1} --log BC_${name}_R1.log ${failed} --nproc ${nproc} 2>&1 | tee -a out_${R1}_BC.log
	BuildConsensus.py -s $R2 ${args_2} --log BC_${name}_R2.log ${failed} --nproc ${nproc} 2>&1 | tee -a out_${R1}_BC.log
	"""
}else{
	"""
	BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc} 2>&1 | tee -a out_${R1}_BC.log
	"""
}


}


process PairAwk_P11 {

input:
 set val(name), file(reads) from g_78_reads0_g_79
 val mate from g_1_mate_g_79

output:
 set val(name), file("*pair-pass.fastq")  into g_79_reads0_g28_12

script:

if(mate=="pair"){
	readArray = reads.toString().split(' ')	
	
	R1 = readArray[0].toString()
	R2 = readArray[1].toString()
	

	
	"""
	BEGINING1=\$(echo $R1|awk '{split(\$0,a,".fa");print a[1];}')
BEGINING2=\$(echo $R2|awk '{split(\$0,a,".fa");print a[1];}')
awk -v out1="\${BEGINING1}_pair-pass.fastq" -v out2="\${BEGINING2}_pair-pass.fastq" 'NR==FNR{
  if(NR%4==1){
    split(\$0,a,"|");
    NAME=a[1];
    split(a[2],b,"=");
    split(a[3],c,"=");
    split(a[4],d,"=");
    split(a[5],e,"=");
    split(a[6],f,"=");
    CONSCOUNT[a[1]]=b[2];
    PRCONS[a[1]]=c[2];
    PRFREQ[a[1]]=d[2];
    BARCODE[a[1]]=e[2];
    BARCODE_COUNT[a[1]]=f[2];
  }
  if(NR%4==2)SEQ[NAME]=\$0;
  if(NR%4==0)QUAL[NAME]=\$0;
  next;
}
NR%4==1{
  flag=0;
  split(\$0,a,"|");
  if(a[1] in SEQ)flag=1;
    split(a[2],b,"=");
    split(a[3],c,"=");
    split(a[4],d,"=");
    split(a[5],e,"=");
    split(a[6],f,"=");
}
flag==1{
  if(NR%4==1){
  
   split(c[2], isotypes, ",")
   delete uniq_isotypes
   for (i in isotypes) {
        uniq_isotypes[isotypes[i]] = 1
    }
   collapsed_isotypes = ""
   sep = ""
   for (key in uniq_isotypes) {
        collapsed_isotypes = collapsed_isotypes sep key
        sep = ","
   }
   print a[1] "|CONSCOUNT=" CONSCOUNT[a[1]] "|PRCONS=" PRCONS[a[1]] "|PRFREQ=" PRFREQ[a[1]] "|ISOTYPE=" collapsed_isotypes  > out1;
   print a[1] "|CONSCOUNT=" CONSCOUNT[a[1]] "|PRCONS=" PRCONS[a[1]] "|PRFREQ=" PRFREQ[a[1]] "|ISOTYPE=" collapsed_isotypes  > out2;
    next;
  }
  if(NR%4==2){
    print SEQ[a[1]] > out1;
    print \$0 > out2;
    next;
  }
  if(NR%4==3){
    print "+" > out1;
    print \$0 > out2;
    next;
  }
  if(NR%4==0){
    print QUAL[a[1]] > out1;
    print \$0 > out2;
    next;
  }
} ' $R1 $R2
	
	"""
}else{
	
	"""
	echo -e 'PairAwk works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_align_assemble_pairs {

input:
 set val(name),file(reads) from g_79_reads0_g28_12
 val mate from g_1_mate_g28_12

output:
 set val(name),file("*_assemble-pass.f*")  into g28_12_reads0_g_84
 set val(name),file("AP_*")  into g28_12_logFile1_g28_15
 set val(name),file("*_assemble-fail.f*") optional true  into g28_12_reads_failed2_g73_12
 set val(name),file("out*")  into g28_12_logFile33

script:
method = params.Assemble_pairs_align_assemble_pairs.method
coord = params.Assemble_pairs_align_assemble_pairs.coord
rc = params.Assemble_pairs_align_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_align_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_align_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_align_assemble_pairs.failed
fasta = params.Assemble_pairs_align_assemble_pairs.fasta
nproc = params.Assemble_pairs_align_assemble_pairs.nproc
alpha = params.Assemble_pairs_align_assemble_pairs.alpha
maxerror = params.Assemble_pairs_align_assemble_pairs.maxerror
minlen = params.Assemble_pairs_align_assemble_pairs.minlen
maxlen = params.Assemble_pairs_align_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_align_assemble_pairs.scanrev
minident = params.Assemble_pairs_align_assemble_pairs.minident
evalue = params.Assemble_pairs_align_assemble_pairs.evalue
maxhits = params.Assemble_pairs_align_assemble_pairs.maxhits
fill = params.Assemble_pairs_align_assemble_pairs.fill
aligner = params.Assemble_pairs_align_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_align_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_align_assemble_pairs.// dbexec
gap = params.Assemble_pairs_align_assemble_pairs.gap
usearch_version = params.Assemble_pairs_align_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_align_assemble_pairs.assemble_reference
head_seqeunce_file = params.Assemble_pairs_align_assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains(""+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Assemble_pairs_align_parse_log_AP {

input:
 set val(name),file(log_file) from g28_12_logFile1_g28_15
 val mate from g_1_mate_g28_15

output:
 file "*table.tab"  into g28_15_logFile0_g28_25, g28_15_logFile0_g28_19

script:
field_to_parse = params.Assemble_pairs_align_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_align_report_assemble_pairs {

input:
 file log_files from g28_15_logFile0_g28_19
 val matee from g_1_mate_g28_19

output:
 file "*.rmd"  into g28_19_rMarkdown0_g28_25



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	name = assemble-"_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">AP_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_align_presto_render_rmarkdown {

input:
 file rmk from g28_19_rMarkdown0_g28_25
 file log_file from g28_15_logFile0_g28_25

output:
 file "*.html" optional true  into g28_25_outputFileHTML00
 file "*csv" optional true  into g28_25_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_parse_log_MP {

input:
 val mate from g_1_mate_g38_9
 set val(name), file(log_file) from g38_11_logFile2_g38_9

output:
 file "*table.tab"  into g38_9_logFile0_g38_12, g38_9_logFile0_g38_19

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_try_report_maskprimer {

input:
 file primers from g38_9_logFile0_g38_12
 val mate from g_1_mate_g38_12

output:
 file "*.rmd"  into g38_12_rMarkdown0_g38_19


shell:

if(mate=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray[0]
	primers_2 = readArray[1]
	name = primers_1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	
	primer_log1_error <- any(is.na(primer_log_1[['ERROR']]))
	primer_log2_error<- any(is.na(primer_log_2[['ERROR']]))
	
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	# check the error column exists 
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray[0]
	name = primers - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_presto_render_rmarkdown {

input:
 file rmk from g38_12_rMarkdown0_g38_19
 file log_file from g38_9_logFile0_g38_19

output:
 file "*.html" optional true  into g38_19_outputFileHTML00
 file "*csv" optional true  into g38_19_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process make_report_pipeline_cat_all_file {

input:
 set val(name), file(log_file) from g38_11_logFile3_g61_0

output:
 set val(name), file("all_out_file.log")  into g61_0_logFile0_g61_2, g61_0_logFile0_g61_10

script:
readArray = log_file.toString()

"""

echo $readArray
cat out* >> all_out_file.log
"""

}


process make_report_pipeline_report_pipeline {

input:
 set val(name), file(log_files) from g61_0_logFile0_g61_2

output:
 file "*.rmd"  into g61_2_rMarkdown0_g61_10


shell:

readArray = log_files.toString().split(' ')
R1 = readArray[0]

'''
#!/usr/bin/env perl


my $script = <<'EOF';


```{r, message=FALSE, echo=FALSE, results="hide"}
# Setup
library(prestor)
library(knitr)
library(captioner)

plot_titles <- c("Read 1", "Read 2")
if (!exists("tables")) { tables <- captioner(prefix="Table") }
if (!exists("figures")) { figures <- captioner(prefix="Figure") }
tables("count", 
       "The count of reads that passed and failed each processing step.")
figures("steps", 
        paste("The number of reads or read sets retained at each processing step. 
               Shown as raw counts (top) and percentages of input from the previous 
               step (bottom). Steps having more than one column display individual values for", 
              plot_titles[1], "(first column) and", plot_titles[2], "(second column)."))
```

```{r, echo=FALSE}
console_log <- loadConsoleLog(file.path(".","!{R1}"))
```

# Summary of Processing Steps

```{r, echo=FALSE}
count_df <- plotConsoleLog(console_log, sizing="figure")

df<-count_df[,c("task", "pass", "fail")]

write.csv(df,"pipeline_statistics.csv") 
```

`r figures("steps")`

```{r, echo=FALSE}
kable(count_df[c("step", "task", "total", "pass", "fail")],
      col.names=c("Step", "Task", "Input", "Passed", "Failed"),
      digits=3)
```

`r tables("count")`


EOF
	
open OUT, ">pipeline_statistic_!{name}.rmd";
print OUT $script;
close OUT;

'''
}


process make_report_pipeline_presto_render_rmarkdown {

input:
 file rmk from g61_2_rMarkdown0_g61_10
 file log_file from g61_0_logFile0_g61_10

output:
 file "*.html" optional true  into g61_10_outputFileHTML00
 file "*csv" optional true  into g61_10_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Assemble_pairs_reference_assemble_pairs {

input:
 set val(name),file(reads) from g28_12_reads_failed2_g73_12
 val mate from g_1_mate_g73_12

output:
 set val(name),file("*_assemble-pass.f*")  into g73_12_reads0_g_84
 set val(name),file("AP_*")  into g73_12_logFile1_g73_15
 set val(name),file("*_assemble-fail.f*") optional true  into g73_12_reads_failed22
 set val(name),file("out*")  into g73_12_logFile33

script:
method = params.Assemble_pairs_reference_assemble_pairs.method
coord = params.Assemble_pairs_reference_assemble_pairs.coord
rc = params.Assemble_pairs_reference_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_reference_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_reference_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_reference_assemble_pairs.failed
fasta = params.Assemble_pairs_reference_assemble_pairs.fasta
nproc = params.Assemble_pairs_reference_assemble_pairs.nproc
alpha = params.Assemble_pairs_reference_assemble_pairs.alpha
maxerror = params.Assemble_pairs_reference_assemble_pairs.maxerror
minlen = params.Assemble_pairs_reference_assemble_pairs.minlen
maxlen = params.Assemble_pairs_reference_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_reference_assemble_pairs.scanrev
minident = params.Assemble_pairs_reference_assemble_pairs.minident
evalue = params.Assemble_pairs_reference_assemble_pairs.evalue
maxhits = params.Assemble_pairs_reference_assemble_pairs.maxhits
fill = params.Assemble_pairs_reference_assemble_pairs.fill
aligner = params.Assemble_pairs_reference_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_reference_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_reference_assemble_pairs.// dbexec
gap = params.Assemble_pairs_reference_assemble_pairs.gap
usearch_version = params.Assemble_pairs_reference_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_reference_assemble_pairs.assemble_reference
head_seqeunce_file = params.Assemble_pairs_reference_assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains(""+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process combine_fasta_p11 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /R1.f.*$/) "reads_before_split/$filename"}
input:
 set val(name),file(reads1) from g28_12_reads0_g_84
 set val(name),file(reads2) from g73_12_reads0_g_84

output:
 set val(name),file("R1.f*")  into g_84_fastaFile0_g_86

"""
#shell example: 

#!/bin/sh 

cat ${reads1} ${reads2} > R1.fasta
"""
}


process split_seq {

input:
 set val(name),file(reads) from g_84_fastaFile0_g_86

output:
 set val(name), file("*_atleast-*.fast*")  into g_86_fastaFile0_g_72
 set val(name),file("out*") optional true  into g_86_logFile11

script:
field = params.split_seq.field
num = params.split_seq.num
fasta = params.split_seq.fasta

readArray = reads.toString()

if(num!=0){
	num = " --num ${num}"
}else{
	num = ""
}

fasta = (fasta=="false") ? "" : "--fasta"

"""
SplitSeq.py group -s ${readArray} -f ${field} ${num} ${fasta} >> out_${readArray}_SS.log
"""

}


process split_constant {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /light\/.*.fasta$/) "reads/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /heavy\/.*.fasta$/) "reads/$filename"}
input:
 set val(name),file(reads) from g_86_fastaFile0_g_72

output:
 set name, file("light/*.fasta") optional true  into g_72_fastaFile00
 set name, file("heavy/*.fasta") optional true  into g_72_fastaFile11

script:
split_col = params.split_constant.split_col

"""
#!/bin/sh 
mkdir heavy
mkdir light
#awk '/^>/{f=""; split(\$0,b,"${split_col}="); if(substr(b[2],2,3)=="IGK"){f="light/${name}_IGK.fasta"} else {if(substr(b[2],2,3)=="IGL"){f="light/${name}_IGL.fasta"} else {f="heavy/${name}.fasta"}}; print \$0 > f ; next } {print \$0 > f} ' ${reads}
awk '/^>/{f=""; split(\$0,b,"${split_col}="); if(substr(b[2],2,3)=="IGK"){f="light/${name}.fasta"} else {if(substr(b[2],2,3)=="IGL"){f="light/${name}.fasta"} else {f="heavy/${name}.fasta"}}; print \$0 > f ; next } {print \$0 > f} ' ${reads}
"""

}


process Assemble_pairs_reference_parse_log_AP {

input:
 set val(name),file(log_file) from g73_12_logFile1_g73_15
 val mate from g_1_mate_g73_15

output:
 file "*table.tab"  into g73_15_logFile0_g73_25, g73_15_logFile0_g73_19

script:
field_to_parse = params.Assemble_pairs_reference_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_reference_report_assemble_pairs {

input:
 file log_files from g73_15_logFile0_g73_19
 val matee from g_1_mate_g73_19

output:
 file "*.rmd"  into g73_19_rMarkdown0_g73_25



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	name = assemble-"_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">AP_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_reference_presto_render_rmarkdown {

input:
 file rmk from g73_19_rMarkdown0_g73_25
 file log_file from g73_15_logFile0_g73_25

output:
 file "*.html" optional true  into g73_25_outputFileHTML00
 file "*csv" optional true  into g73_25_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process metadata {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.json$/) "metadata/$filename"}

output:
 file "*.json"  into g_69_jsonFile00

script:
metadata = params.metadata.metadata
"""
#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

data <- read_json("${metadata}") 

versions <- lapply(1:length(data), function(i){
	
	docker <- data[i]
	tool <- names(data)[i]
	
	if(grepl("Custom", docker)){
		ver <- "0.0"
	}else{
		ver <- system(paste0(tool," --version"), intern = TRUE)
		ver <- gsub(paste0(tool,": "), "", ver)
	}
	ver
	
})

names(versions) <- names(data)

json_data <- list(
  sample = list(
    data_processing = list(
      preprocessing = list(
        software_versions = versions
	   )
	 )
  )
)

# Convert to JSON string without enclosing scalar values in arrays
json_string <- toJSON(json_data, pretty = TRUE, auto_unbox = TRUE)
print(json_string)
# Write the JSON string to a file
writeLines(json_string, "pre_processed_metadata.json")
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
