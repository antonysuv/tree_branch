#!/bin/bash
batch=keras_cnn_mphate.sh

if [ -e $batch ]
then
    continue
fi
echo '#!/bin/bash' > $batch
echo '#SBATCH --job-name=cnn_mphate' >> $batch
echo '#SBATCH --ntasks=1' >> $batch
echo '#SBATCH --cpus-per-task=1' >> $batch
echo '#SBATCH --mem=60G' >> $batch
echo '#SBATCH --time=100:00:00' >> $batch
echo '#SBATCH --partition=volta-gpu' >> $batch
echo '#SBATCH --output=run-cnn-%j.log' >> $batch
echo '#SBATCH --gres=gpu:1' >> $batch
echo '#SBATCH --qos=gpu_access' >> $batch
echo 'source activate tf-gpu' >> $batch
echo '~/scripts/keras_CNN_MPHATE_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl TRAIN/train_simulation.Y.txt --tel TEST/test_simulation.Y.txt' >> $batch
sbatch $batch
