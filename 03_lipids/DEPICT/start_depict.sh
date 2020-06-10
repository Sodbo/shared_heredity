{
for i in 1287003 1287004 SH 1287001-SH 1287003-SH 1287004-SH
do
echo Starting calculation for trait $i
python /home/common/projects/depict_software/DEPICT_v194/src/python/depict.py ID_${i}.cfg
done
} 2>&1 | tee /mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/DEPICT/logs.txt
