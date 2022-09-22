#!/bin/bash
batch=keras_mlp.sh

echo '#!/bin/bash' > $batch
echo '#SBATCH --job-name=mlp' >> $batch
echo '#SBATCH --ntasks=1' >> $batch
echo '#SBATCH --cpus-per-task=1' >> $batch
echo '#SBATCH --mem=50G' >> $batch
echo '#SBATCH --time=100:00:00' >> $batch
echo '#SBATCH --partition=dschridelab' >> $batch
echo '#SBATCH --output=run-mlp-%j.log' >> $batch
echo '#SBATCH --gres=gpu:1' >> $batch
echo '#SBATCH --constraint=rhel8' >> $batch
echo 'source activate tf-gpu2' >> $batch
echo 'python ~/scripts/keras_MLP_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl TRAIN/train_simulation.Y.txt --tel TEST/test_simulation.Y.txt --trans sqrt' >> $batch
sbatch $batch
