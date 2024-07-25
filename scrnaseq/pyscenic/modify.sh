for i in `ls /camp/home/hungm/template_pips/pyscenic | grep 'set'`
do
   set_ipath='iteration_path = "/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/data/functional/pyscenic/iterations/"'
   sed -i "1s|.*|${set_ipath}|" /camp/home/hungm/template_pips/pyscenic/${i}/pyscenic.py
   cd /camp/home/hungm/template_pips/pyscenic/${i}/
   sbatch /camp/home/hungm/template_pips/pyscenic/${i}/pyscenic.sh
done

