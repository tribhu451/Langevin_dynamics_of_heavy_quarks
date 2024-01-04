#!/bin/bash


check_directory_existence()
{
    local DIR="$1"
    if [ -d "$DIR"  ]
    then
        echo "[ $DIR ]  exist ..."
    else
        echo "[ $DIR ]  does not exist."
        exit 1
    fi
}

check_file_existence()
{
    local FILE="$1"
    if [ -f "$FILE"  ]
    then
        echo "[ $FILE ]  exist ..."
    else
        echo "[ $FILE ]  does not exist."
        exit 1
    fi
}

file_existence_flag()
{
    NUM=-999
    local FILE="$1"
    if [ -f "$FILE"  ]
    then
        NUM=1
    else
        NUM=0
    fi

    echo $NUM
}



#====================================================================================
#=================                      MAIN                 ========================
#====================================================================================




#================= Generating input file and RUN directories ========================
echo -e "Generating the input files and creating the root RUN directories ... "
python3 generate_folder_for_diff_param_run.py
echo -e "All input files moved into respective directories ... "
#====================================================================================


#=================  All source codes copied to root RUN directories =================
NUM_OF_FOLDERS=$(ls -dq *RUN* | wc -l)
echo -e "Total runs to be given : $NUM_OF_FOLDERS"

check_directory_existence "LDHQ/"
check_directory_existence "pythia8309/"
check_file_existence "run_langevin.sh"

# Place the source folder and input parameter files to desired place
for((i=0; i<$NUM_OF_FOLDERS; i++)); do
cp -r LDHQ/ RUN$i
cp -r pythia8309/ RUN$i
cp -r run_langevin.sh RUN$i
mv RUN$i/input_parameters RUN$i/LDHQ
done



#============================== LOOP OVER RUNS ======================================
for((i=0; i<$NUM_OF_FOLDERS; i++)); do
cd RUN$i

STARTFILE=0
ENDFILE=50

# Generate 50 pythia events
for((j=$STARTFILE; j<$ENDFILE; j++)); do
cp -r pythia8309/ PYTHIA$j
cd PYTHIA$j
cd examples/
make main01
nohup ./main01 > run.log &
#JID=$!
#wait $JID
cd ../
cd ../
done


sleep 180












for((j=$STARTFILE; j<$ENDFILE; j++)); do
cp -r LDHQ/ FILE$j
mv PYTHIA$j FILE$j/pythia8309
check_file_existence "run_langevin.sh"
cp -r run_langevin.sh FILE$j
cd FILE$j
nohup ./run_langevin.sh > run_langevin.log  &
cd ..
done



#======= Wait till all langevin run finish ====================
count=0
while [ $count -lt $ENDFILE ]
do
for((jj=$STARTFILE; jj<$ENDFILE; jj++)); do
NUM=$(file_existence_flag "FILE$jj/run_finished_flag.dat")
count=$(($count+$NUM))
done
if [ $count -lt $ENDFILE ]
then
        echo -e "Langevin running ... "
        date
        count=0
        sleep 60
else
        echo -e "Langevin run finished ... "
        date
fi
done
#===============================================================


#========= conversion to binary files and File Renaming =======
check_file_existence "file_renaming_and_event_count_info.txt"
filename='file_renaming_and_event_count_info.txt'
LineNo=0
while read line; do

if [ $LineNo -eq 0 ]
then
ReNameTo=$line
fi

LineNo=$((LineNo+1))
done < $filename


#========= Collect all binary files =============
for((j=$STARTFILE; j<$ENDFILE; j++)); do
mv FILE$j/heavy_hadron_list.bin heavy_hadron_list_set_$j.bin
mv FILE$j/heavy_quark_list.bin heavy_quark_list_set_$j.bin
mv FILE$j/init_quark_list.bin  init_quark_list_$j.bin
done

rm -rf FILE*
rm -rf pythia8*
rm -rf LDHQ

cd ..
# out of RUN folder now
mv RUN$i $ReNameTo
done













