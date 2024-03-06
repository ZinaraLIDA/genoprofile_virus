rule joinFile :
	input:
		visualizeResultsDone = expand("{outputPath}/visualizeResults/{prefix}_chr_{index}.done", outputPath=outputPath, prefix=prefix, index=index)

	params :
		allVirus = expand("{outputPath}/filterSummary/{prefix}_chr_{index}_markdup_all_virus.tsv", outputPath=outputPath, prefix=prefix, index=index)
	
	log :
		f"{outputPath}/pipelineLog/joinedFiles/joinedFiles.log"
	
	output :
		f"{outputPath}/joinedFiles/AllFilesJoined.done"
	
	shell :
		"""
		BASEDIR=$(dirname {output})
		for table in {params.allVirus}
			do
				if [ -f $table ]
					then
					prefixName=$(echo $(basename $table) | sed "s/\..*//")
					sample=$(echo $(basename $table) | sed "s/_chr.*//")
					sed "1s/$/\tsample/" $table > $BASEDIR/$prefixName.modified.tsv
					sed -i -e "2,\$s/\$/\\t$sample/" $BASEDIR/$prefixName.modified.tsv
					if [ ! -f $BASEDIR/all_sample_all_virus.tsv ]
						then
						cat $BASEDIR/$prefixName.modified.tsv > $BASEDIR/all_sample_all_virus.tsv
					else
						tail -n +2 $BASEDIR/$prefixName.modified.tsv >> $BASEDIR/all_sample_all_virus.tsv
					fi
				fi
			done >> {log} 2>&1
		touch {output}
		"""