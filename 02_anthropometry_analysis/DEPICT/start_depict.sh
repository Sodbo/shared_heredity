{
for i in 191 192 193 194 199 SH 202 203 204 205 206
do
python /home/common/projects/depict_software/DEPICT_v194/src/python/depict.py ID_${i}.cfg
done
} 2>&1 | tee /mnt/polyomica/projects/shared_heredity/data/01_anthropometry_results/five_traits/DEPICT/logs_ID_191-199.txt
