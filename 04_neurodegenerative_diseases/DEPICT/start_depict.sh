{
for i in bip bip-sh mdd mdd-sh scz scz-sh happiness happiness-sh sh
do
echo Starting calculation for trait ${i}
python /home/common/projects/depict_software/DEPICT_v194/src/python/depict.py ${i}.cfg
done
} 2>&1 | tee /mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases/DEPICT/log_pgc.txt
