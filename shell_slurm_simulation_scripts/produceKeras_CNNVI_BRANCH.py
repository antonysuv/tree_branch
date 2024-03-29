#!/bin/bash
batch=keras_cnnvi.sh

if [ -e $batch ]
then
    continue
fi
echo '#!/bin/bash' > $batch
echo '#SBATCH --job-name=cnnvi' >> $batch
echo '#SBATCH --ntasks=1' >> $batch
echo '#SBATCH --cpus-per-task=1' >> $batch
echo '#SBATCH --mem=60G' >> $batch
echo '#SBATCH --time=100:00:00' >> $batch
echo '#SBATCH --partition=volta-gpu' >> $batch
echo '#SBATCH --output=run-cnnvi-%j.log' >> $batch
echo '#SBATCH --gres=gpu:1' >> $batch
echo '#SBATCH --qos=gpu_access' >> $batch
echo 'source activate tf-gpu' >> $batch
echo 'python3.9 ~/scripts/keras_CNNVI_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl TRAIN/train_simulation.Y.txt --tel TEST/test_simulation.Y.txt' >> $batch
sbatch $batch
