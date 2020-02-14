#The script for grep clumped shared heredity SNPs

path='/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/GWAS/scaled_filtered/'
for rs in  rs1077514  rs11802413  rs12748152  rs17513135  rs11591147  rs12066643  rs10789119  rs6603981  rs12740374  rs267733  rs2642438  rs1321257  rs558971  rs1107851  rs1367117  rs2194562  rs1260326  rs4401177  rs6544713  rs2710642  rs17508045  rs2030746  rs6730157  rs4988235  rs13389219  rs2287623  rs11694172  rs1174604  rs1250229  rs2972146  rs11563251  rs7616006  rs7640978  rs13315871  rs17345563  rs645040  rs10513688  rs6831256  rs2869433  rs442177  rs13133548  rs9686661  rs7703051  rs7727150  rs4530754  rs6882076  rs2294261  rs2274089  rs1800562  rs3132625  rs2247056  rs3130679  rs9391858  rs2814982  rs3800406  rs2758886  rs998584  rs17789218  rs868943  rs719726  rs9376090  rs634869  rs11753995  rs6935921  rs1997243  rs6968554  rs12670798  rs4722551  rs2073547  rs2909969  rs714052  rs38855  rs287621  rs2976940  rs9987289  rs6995541  rs1062219  rs4921914  rs12679834  rs10102164  rs4738684  rs2737252  rs2954022  rs7832643  rs3780181  rs581080  rs7033354  rs10757056  rs1883025  rs579459  rs1781930  rs10904908  rs970548  rs7897379  rs2068888  rs1129555  rs10832962  rs7932354  rs10501321  rs4752805  rs1535  rs10790162  rs12282721  rs11603023  rs7117842  rs11220462  rs11613352  rs10861661  rs12321904  rs3184504  rs17630235  rs1169288  rs10773003  rs11057408  rs9534262  rs17061870  rs1341267  rs8017377  rs2412710  rs492571  rs10468017  rs2652834  rs3198697  rs749671  rs9930333  rs247616  rs16942887  rs17343777  rs2000999  rs8044476  rs2925979  rs8069974  rs314253  rs4791641  rs9972882  rs8077889  rs7225700  rs1801689  rs12602912  rs2886232  rs9894524  rs2156552  rs7248104  rs4804311  rs10403668  rs6511720  rs4808802  rs10401969  rs4808993  rs731839  rs1688030  rs1594895  rs7254892  rs7255743  rs516246  rs1132990  rs103294  rs364585  rs2328223  rs7264396  rs6016381  rs6065311  rs1800961  rs4810479  rs1211644  rs181360  rs5763662  rs138764  rs3761445  rs4253772   
do
	grep "$rs," $path'ID_1287001.csv' >> $path'ID_1287001_SH_clump.txt'
	grep "$rs," $path'ID_1287003.csv' >> $path'ID_1287003_SH_clump.txt'
	grep "$rs," $path'ID_1287004.csv' >> $path'ID_1287004_SH_clump.txt'
	grep -P "$rs\t" '/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/linear_combination/SH_Lipids_GWAS.txt' >> $path'SH_clump.txt'
done
