#!/bin/bash
batch=keras_lstmcnn.sh

if [ -e $batch ]
then
    continue
fi
echo '#!/bin/bash' > $batch
echo '#SBATCH --job-name=lstmcnn' >> $batch
echo '#SBATCH --ntasks=1' >> $batch
echo '#SBATCH --cpus-per-task=1' >> $batch
echo '#SBATCH --mem=60G' >> $batch
echo '#SBATCH --time=100:00:00' >> $batch
echo '#SBATCH --partition=volta-gpu' >> $batch
echo '#SBATCH --output=run-lstmcnn-%j.log' >> $batch
echo '#SBATCH --gres=gpu:1' >> $batch
echo '#SBATCH --qos=gpu_access' >> $batch
echo 'source activate tf-gpu' >> $batch
echo 'python3.9 ~/scripts/keras_LSTM_CNN_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl TRAIN/train_simulation.Y.txt --tel TEST/test_simulation.Y.txt' >> $batch
sbatch $batch
