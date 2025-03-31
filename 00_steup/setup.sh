ssh aagrawal@log2.wynton.ucsf.edu

cd /gladstone/jain/boinformatics-collaboration/

mkdir sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025

mkdir sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/assets
mkdir sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/data
mkdir sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/results
mkdir sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/scripts
mkdir sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp

vi sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/METADATA.txt
# fill out all the details

exit

# start the data transfer from 3-=day-transfer area to jain lab folder
sft list-servers

# login to the Gladstone data transfer node
sft ssh gsdt02

# start tmux session
tmux

# copy the data 
cd /gladstone/30-day-transfer-area/Bioinformatics\ Core/jain-lab-rna-flex-2025/
cp -r Jain_Flex_RNA_2025/ /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/data/
