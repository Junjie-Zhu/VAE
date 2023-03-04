# Enhancing Conformational Sampling for Intrinsically Disordered and Ordered Proteins by Variational Autoencoder

**Overview of VAE:** Variational Autoencoder (VAE) generates similar conformations based on a given MD trajectory, thus can be used to enhance the sampling of the diverse conformation for intrinsically disordered proteins (IDPs) and structural proteins.

### Dependencies

- Compatible with Python 3.x.
- Dependencies can be installed using the `requirements.txt` file.

### Training Model

1. Install all the requirements by executing `pip install -r requirements.txt.`
2. Install required protein `.pdb` processing library `biobox` from [this](https://github.com/Degiacomi-Lab/biobox) github repository.
3. Next execute `python preprocess.py pdb split` which aligns the conformations from input `pdb` and creates VAE-required dataset at certain `split`.
4. To start a training run:

```bash
python packed_vae.py pdb
```

where `pdb` denotes the filename of input MD trajectory.

5. To calculate the RMSD between generated conformations and original ones:

```shell
python rmsd_min.py pdb
```

### K-cluster

To conduct cluster analysis on generated conformations, please use the scripts in `k-cluster` folder.

1. Cluster analysis requires `MMTSB` tool set in this research, install the tool set from [this](https://mmtsb.org/) website.
2. execute `list.pl` to get file list for further analysis.
3. execute `run.sh` to conduct cluster analysis.

### Chemical Shift Calculation

To calculate secondary chemical shift of generated conformations, please use the scriptes in `chemical-shift` folder.

1. Calculation requires `SPARTA+` program from [this](https://spin.niddk.nih.gov/bax/software/SPARTA+/) website.
2. execute `split_pdb.py` to calculate chemical shift.
3. execute `avg_cs.py` to get average chemical shift of target conformation ensemble.

