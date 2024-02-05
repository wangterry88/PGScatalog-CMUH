#!/bin/bash

# Prepare require folder
echo "==========================================================================================="
echo ""
echo " Step 0: Preparing Required folders and download PGS information files ....."
echo ""
echo "==========================================================================================="
mkdir ./Result/${PROJECT}
mkdir ./Result/${PROJECT}/0.PGS-RawData
mkdir ./Result/${PROJECT}/1.Process
mkdir ./Result/${PROJECT}/2.Done
mkdir ./Result/${PROJECT}/3.PRSresult
mkdir ./Result/${PROJECT}/4.PRSresult-all
mkdir ./Result/${PROJECT}/4.PRSresult-all/done/
mkdir ./Result/${PROJECT}/5.PRSresult-merge/
mkdir ./Result/${PROJECT}/Intermediate-files

for PGSlist in $(cat ${LIST} | sed 's/\r//g') ; do  mkdir ./Result/${PROJECT}/3.PRSresult/${PGSlist}  ; done

wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata.xlsx --no-check-certificate -P ./Result/${PROJECT}/

echo "==========================================================================================="
echo ""
echo " Step 1: Preparing Required Packages ....."
echo ""
echo "==========================================================================================="

# 0. Prepare require packages

pip install csvkit

echo "==========================================================================================="
echo ""
echo " Step 2: Downloading PGS Scoring File form PGS Catalog ....."
echo ""
echo "==========================================================================================="

# 1. Download the data

for PGSlist in $(cat ${LIST} | sed 's/\r//g') ; do wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${PGSlist}/ScoringFiles/Harmonized/${PGSlist}_hmPOS_GRCh38.txt.gz --no-check-certificate -P ./Result/${PROJECT}/0.PGS-RawData ; done

echo "==========================================================================================="
echo ""
echo " Step 3: Processing Scoring File ....."
echo ""
echo "==========================================================================================="

# 2. Unzip tha data

cd ./Result/${PROJECT}/0.PGS-RawData

gzip -dv ./*.gz

# 3. Get the file name to list

ls |cut -d "_" -f 1 > ./file_list.txt

# 4. Get file name and add column + change headers

cd ../../../

for PGSlist in $(cat ${LIST} | sed 's/\r//g') ; do awk '{print}' ./Result/${PROJECT}/0.PGS-RawData/${PGSlist}_hmPOS_GRCh38.txt|grep -v '^#'|csvcut -t -c hm_chr,hm_pos,effect_allele,effect_weight|csvformat -T | awk '$(NF+1) = 1'|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'|sed -e '1s/hm_chr/CHR/' -e '1s/hm_pos/BP/' -e '1s/effect_weight/Effect/' -e '1s/1/P/'| sed -e '1s/effect_allele/A1/'|awk '{print $1"\t"$2"\t""chr"$1"_"$2"\t"$3"\t"$4"\t"$5}'|sed -e '1s/chrCHR_BP/SNP/' > ./Result/${PROJECT}/1.Process/${PGSlist}.txt ; done

# 5. Process the output file

## Get the correct columns and with header

for PGSlist in $(cat ${LIST} | sed 's/\r//g') ; do awk '$1<23' ./Result/${PROJECT}/1.Process/${PGSlist}.txt|sed '1iCHR\tBP\tSNP\tA1\tEffect\tP'> ./Result/${PROJECT}/2.Done/${PGSlist}.txt ; done

# 6. Run PRS (All)

echo "==========================================================================================="
echo ""
echo " Step 4: Calculating PRS with scoring files ....."
echo ""
echo "==========================================================================================="

for PGSlist in $(cat ${LIST} | sed 's/\r//g') ; do Rscript ./software/PRSice-2/PRSice.R --prsice ./software/PRSice-2/PRSice_linux --base ./Result/${PROJECT}/2.Done/${PGSlist}.txt --target ${BFILE} --bar-levels 1 --stat Effect --beta --A1 A1 --pvalue P --score std --thread 64 --print-snp --seed 123456789 --no-default --chr CHR --bp BP --chr-id C_L --no-regress --no-clump --fastscore --out ./Result/${PROJECT}/3.PRSresult/${PGSlist}/PGSresult_${PGSlist} ; Rscript ./software/PRSice-2/PRSice.R --prsice ./software/PRSice-2/PRSice_linux --base ./Result/${PROJECT}/2.Done/${PGSlist}.txt --target ${BFILE} --bar-levels 1 --stat Effect --beta --A1 A1 --pvalue P --score std --thread 64 --print-snp --seed 123456789 --no-default --chr CHR --bp BP --chr-id C_L --no-regress --no-clump --fastscore --out ./Result/${PROJECT}/3.PRSresult/${PGSlist}/PGSresult_${PGSlist} --extract ./Result/${PROJECT}/3.PRSresult/${PGSlist}/PGSresult_${PGSlist}.valid ; done > /dev/null 2>&1

echo ""
echo "==========================================================================================="
echo ""
echo " Please ignore Error messages from PRSice-2, this will NOT interrupt the Program ....."
echo ""
echo "==========================================================================================="

# 7. Get Full PRS calculation (ABC step have to work together, if no work together PGS"$i".done will be broken)
echo "==========================================================================================="
echo ""
echo " Step 5: Processing calculated PGS files ....."
echo ""
echo "==========================================================================================="

### 7A.

for PGSlist in $(cat ${LIST} | sed 's/\r//g') ; do cp ./Result/${PROJECT}/3.PRSresult/${PGSlist}/PGSresult_${PGSlist}.all_score ./Result/${PROJECT}/4.PRSresult-all/${PGSlist}.txt ; done

### 7B.

cd ./Result/${PROJECT}/4.PRSresult-all

### 7C.

for f in *.txt ; do sed -i "1iFID\tIID\t${f%%.*}" "$f" ; sed '2d' "$f" > ./done/"$f".done ; done

# 8. Make a big table of PRS

cd ../../../

# First, Get file First IID,PRS

awk '{print $2}' ./Result/${PROJECT}/4.PRSresult-all/done/$(cat ${LIST} | head -n 1 | sed 's/\r//g').txt.done > ./Result/${PROJECT}/5.PRSresult-merge/SampleID.txt

# Second, Get file Second to N PRS

for PGSlist in $(cat ${LIST} | sed 's/\r//g') ; do awk '{print $3}' ./Result/${PROJECT}/4.PRSresult-all/done/${PGSlist}.txt.done > ./Result/${PROJECT}/5.PRSresult-merge/${PGSlist}.tomerge.txt ; done

# Third merge the two files

cd ./Result/${PROJECT}/5.PRSresult-merge

file_list=$(for S in *.tomerge.txt ; do echo "$S"; done)

paste --delimiters='\t' SampleID.txt ${file_list} > ../${PROJECT}.PGS-Analysis-done.txt

echo "==========================================================================================="
echo ""
echo " Step 6: PGS Catalog Analysis Done ........"
echo ""
echo "==========================================================================================="
echo "==========================================================================================="
echo ""
echo " Step 7: Moving intermediate files ........"
echo ""
echo "==========================================================================================="
cd ../../../
mv ./Result/${PROJECT}/0.PGS-RawData/ ./Result/${PROJECT}/Intermediate-files
mv ./Result/${PROJECT}/1.Process/ ./Result/${PROJECT}/Intermediate-files
mv ./Result/${PROJECT}/2.Done/ ./Result/${PROJECT}/Intermediate-files
mv ./Result/${PROJECT}/3.PRSresult/ ./Result/${PROJECT}/Intermediate-files
mv ./Result/${PROJECT}/4.PRSresult-all/ ./Result/${PROJECT}/Intermediate-files
mv ./Result/${PROJECT}/5.PRSresult-merge/ ./Result/${PROJECT}/Intermediate-files
echo "==========================================================================================="
echo ""
echo " Moving intermediate files done ........"
echo ""
echo "==========================================================================================="
echo "==========================================================================================="
echo ""
echo " PGS Catalog Analysis Result is in: ./Result/${PROJECT}/${PROJECT}.PGS-Analysis-done.txt"
echo ""
echo "==========================================================================================="
echo ""
echo "==========================================================================================="
echo ""
echo " Do you want to print plot of those PGS result ? "
echo ""
echo "[1] Yes  [2] No  [Other] Exit the Program"
echo ""
read PLOTMODE
echo "==========================================================================================="

if [ $PLOTMODE -eq 1 ]; then
    echo ""
    echo "==========================================================================================="
    echo " Your Select is [1] : Now start the plotting program ...."
    echo "==========================================================================================="
    mkdir ./Result/${PROJECT}/PGSplot
    Rscript ./software/PGS-Catalog/Plot.R ${PROJECT}
    echo "==========================================================================================="
    echo " The plotting program is done ...."
    echo ""
    echo " PGS Catalog Analysis Result and plot is in: ./Result/${PROJECT}/"
    echo ""
    echo " Exit the program .... "
    echo ""
    echo " Bye! Bye! "
    echo "==========================================================================================="
elif [ $PLOTMODE -eq 2 ]; then
    echo "==========================================================================================="
    echo " Your Select is [2] : Now exit the program .... "
    echo ""
    echo " PGS Catalog Analysis Result is in: ./Result/${PROJECT}/"
    echo ""
    echo " Bye! Bye! "
    echo "==========================================================================================="
else 
    echo "==========================================================================================="
    echo " Your Select is [Other] : Now exit the program .... "
    echo ""
    echo " PGS Catalog Analysis Result is in: ./Result/${PROJECT}/"
    echo ""
    echo " Bye! Bye! "
    echo "==========================================================================================="
fi