# 3D MERMAID
3D Multishot Enhanced Recovery Motion Artifact Insensitive Diffusion (3D MERMAID) is a diffusion sequence proposed to improve 2D and 3D multi-slice/slab acquisitions, by enhacing SNR, solving slice/slab boundary artifacts, and improving spin history. The intial implementation includes additional inversion pulse immediately before the excitation to improve signal recovery, and TURBINE readout technique to correct phase errors between shots.
# Sequence simulation
Codes required for Bloch simulation can be found in ***3D_MERMAID_sim*** directory.
***MERMAID_sim.m*** includes examples that reproduces results of the Figures 1 and 3 of the paper.
# Image reconstruction
Image reconstruction pipeline for reconstructing scans acquired using 3D MERMAID sequence can be found in ***Imge_reconstruction*** directory. Currently, the code is adapted for using TWIX data of Siemens scanners.
Two different pipelines are presented as described in the paper:
* ***recon_pipeline.m*** is the pipeline of Figure 2.
* ***recon_pipeline_denoising.m*** is the pipeline in Figure S2 that includes denoising.
## Requirements:
* Siemens mapVBVD reader is required for importing TWIX files.
* BART toolbox is required for estimating coil sensitivity maps, and reconstructing scans using Compressed SENSE:
https://mrirecon.github.io/bart/
* NORDIC was integrated in the pipeline to include denoising. The provided code is adapted from the original implementation:
https://github.com/SteenMoeller/NORDIC_Raw
# Pulse sequence
The 3D MERMAID pulse sequence is available for these versions of Siemens scanners:
* VE11C
* XA30

# References
paper.

Please cite appropriate references for BART and NORDIC if they are used for the image reconstruction.

* For any question about the code and pulse sequence please contact sajjad.feizollah@mail.mcgill.ca
