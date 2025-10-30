```bash
conda create -n rocketenv python=3.14
conda activate rocketenv
mamba install hatchling build pytest numpy
pip install -e .
