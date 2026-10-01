## CHG_CD
Unsupervised SAR Image Change Detection using Coupled Heterogeneous Graph with Patch-level and Superpixel-level Attributes

## Introduction
Synthetic aperture radar (SAR) is an indispensable data source for change detection tasks due to its ability to operate under all-weather and all-illumination conditions. However, accurate extraction of changes remains challenging due to the inherent speckle noise in SAR images, coupled with the complex scattering characteristics of ground objects. To these ends, we proposes a coupled heterogeneous graph (CHG) for unsupervised SAR image change detection. In particular, each SAR image undergoes simultaneous segmentation into both a series of superpixels and a set of overlapping patches. Subsequently, a patch-level similarity graph (PSG) is constructed to capture detailed information, while a superpixel-level affinity graph (SAG) is developed to represent regional information. According to the affiliation between patches and superpixels, the PSG is coupled with the SAG, thereby constructing the coupled heterogeneous graph structure. With the support of CHG, the difference image (DI) generation relies on the matching of the bitemporal graphs, as well as the patch and superpixel attributes. Finally, the graph cuts algorithm is applied to segment the DI into changed and unchanged regions. Experimental results on three real SAR datasets demonstrate that the proposed method consistently outperforms competing approaches and presents a good candidate for change detection.

## Citation
If you use this code for your research, please cite our paper. Thank you!

@ARTICLE{**, author={Qiongjun Fu, Jun Wang, Sanku Niu, Xiangyu Yang, Chunyang Li},
journal={Remote Sensing Letters},
title={Unsupervised SAR Image Change Detection using Coupled Heterogeneous Graph with Patch-level and Superpixel-level Attributes},
year={2026}.}

## Running
Run the CHG_CD demo files (tested in Matlab 2024b)!

If you have any queries, please contact me (36110@qzc.edu.cn).
