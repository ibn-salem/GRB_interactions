#!/bin/bash

#=======================================================================
#
# This script is suppost do document downlaods of public data in the data 
# folder. It should be excuted from within the /data directory. 
#
#=======================================================================

# set some variables here:
BIN=../bin
mkdir -p ${BIN}

#=======================================================================
# General genome assembly based data and tools from UCSC:
#=======================================================================

#-----------------------------------------------------------------------
# UCSC liftover data and tool
#-----------------------------------------------------------------------

# UCSC liftover chains
mkdir -p UCSC
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

gunzip UCSC/*.gz

# download liftOver tool from UCSC:
wget -P ${BIN} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x ${BIN}/liftOver

#=======================================================================
# Convert GRBs and CNES from hg38 to hg19
#=======================================================================

BED_FILES="
GRB_Interactions/grbsHg38.bed
GRB_Interactions/human_dog_cnes_hg38.bed
GRB_Interactions/human_mouse_cnes_hg38.bed
"

for BED in $BED_FILES ; do

  # convert mouse domains from mm9 to mm10 assembly
  ${BIN}/liftOver \
  	${BED}  \
  	UCSC/hg38ToHg19.over.chain \
  	${BED}.hg19.bed \
  	${BED}_unmapped.bed

done

#=======================================================================
# Capture Hi-C data from Mifsud2015
#=======================================================================
mkdir -p Mifsud2015
wget -P Mifsud2015 http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip
unzip Mifsud2015/E-MTAB-2323.additional.1.zip -d Mifsud2015

#=======================================================================
# Vista Enhancers
#=======================================================================
mkdir -p VistaEnhancers
wget -O VistaEnhancers/vista.fasta "https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page_size=100;show=1;search.result=yes;form=search;search.form=no;action=search;search.sequence=1"

cat VistaEnhancers/vista.fasta \
  |sed 's/<pre>//g' \
  |sed 's/<\/pre>//g' \
  |grep ">" \
  > VistaEnhancers/vista_formated.header.txt
