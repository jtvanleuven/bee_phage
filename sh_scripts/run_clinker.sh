#!/bin/bash

source /home/james/miniconda3/etc/profile.d/conda.sh


for d in *.fasta; do
        name=$(echo $d | sed "s/.fasta//g")
        echo "Processing cluster "$name".............................................................................................."
        cd $name
        for f in *.fasta; do
			phage=$(echo $f | sed "s/.fasta//g")
			echo "Running prokka on "$phage".............................................................................................."
			fasta=$(echo "${phage}.fasta")
			cmdstring1=$(echo "prokka $fasta --outdir prokka$phage --prefix $phage --force")
			conda activate prokka
			eval $cmdstring1
		done
	done
		
		mkdir gbks
		cp prokka*/*.gbk gbks/
		cmdstring2=$(echo "clinker gbks/*.gbk -p $name.html -o $name.clinkerout -f")
		conda activate clinker
		eval $cmdstring2
		cd ..  
