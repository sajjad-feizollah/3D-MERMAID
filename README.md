# 3D MERMAID
3D Multishot Enhanced Recovery Motion Artifact Insensitive Diffusion (3D MERMAID) is a diffusion sequence designed to improve 2D and 3D multi-slice/slab acquisitions by enhancing SNR per unit time, eliminating slice/slab boundary artifacts, improving B1 uniformity, and addressing spin history effects. The initial implementation includes an additional inversion pulse immediately before excitation to improve signal recovery and the TURBINE readout technique to correct phase errors between shots.
# Sequence simulation
The code required for Bloch simulations can be found in the ***3D_MERMAID_sim*** directory.
***MERMAID_sim.m*** includes examples that reproduce the results of Figures 1 and 3 from the paper.
# Image reconstruction
The image reconstruction pipeline for processing scans acquired using the 3D MERMAID sequence can be found in the ***Imge_reconstruction*** directory. Currently, the code is adapted for use with TWIX data from Siemens scanners.
Two different pipelines are provided, as described in the paper:
* ***recon_pipeline.m:*** The image reconstruction pipeline for Figure 2.
* ***recon_pipeline_denoising.m:*** The image reconstruction pipeline for Figure S2, which includes denoising.
## Requirements:
* The Siemens mapVBVD reader is required for importing TWIX files.
* The BART toolbox is required for estimating coil sensitivity maps and reconstructing scans using Compressed SENSE:
https://mrirecon.github.io/bart/
* NORDIC is integrated into the pipeline for denoising. The provided code is adapted from the original implementation:
https://github.com/SteenMoeller/NORDIC_Raw
# Pulse sequence
The 3D MERMAID pulse sequence is available for the following Siemens scanner versions:
* VE11C
* XA30

# References
**Paper:** [Feizollah, S., Tardif, C. L. (2025). *3D MERMAID: 3D Multi‐shot enhanced recovery motion artifact insensitive diffusion for submillimeter, multi‐shell, and SNR‐efficient diffusion imaging*](https://doi.org/10.1002/mrm.30436)

Link to the dataset: [https://doi.org/10.5683/SP3/ULKLNY](https://doi.org/10.5683/SP3/ULKLNY)

Please cite the appropriate references for BART and NORDIC if they are used for image reconstruction.

For any questions regarding the code or pulse sequence, please contact sajjad.feizollah@mail.mcgill.ca
