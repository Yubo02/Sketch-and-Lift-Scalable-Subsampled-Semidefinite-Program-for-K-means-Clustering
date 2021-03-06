# Sketch-and-Lift-Scalable-Subsampled-Semidefinite-Program-for-K-means-Clustering
This repository contains Matlab codes and corresponding dataset generated by ourselves or obtained from some clustering benchmark data. 

The codes are sourced from the following paper:

- Yubo Zhuang, Xiaohui Chen, and Yun Yang. *Sketch-and-Lift: Scalable Subsampled Semidefinite Program for K-means Clustering.* 2022. arXiv: 2201.08226 [stat.ML].\
  https://arxiv.org/abs/2201.08226
  
  
The benchmark clustering data sets are sourced from the following resources:

- Levine et al. (2015). *Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis.* Cell, 162, pp. 184-197.\
  https://github.com/lmweber/benchmark-data-Levine-32-dim
- M. Rezaei and P. Fränti. *Set-matching methods for external cluster validity.* IEEE Trans. on Knowledge and Data Engineering, 28(8):2173–2186, 2016\
  http://cs.joensuu.fi/sipu/datasets/

# Background
Yubo Zhuang, Xiaohui Chen, and Yun Yang. *Sketch-and-Lift: Scalable Subsampled Semidefinite Program for K-means Clustering.* 
## The abstract of the paper
Semidefinite programming (SDP) is a powerful tool for tackling a wide range of computationally hard problems such as clustering. Despite the high accuracy, semidefinite programs are often too slow in practice with poor scalability on large (or even moderate) datasets. In this paper, we introduce a linear time complexity algorithm for approximating an SDP relaxed K-means clustering. The proposed *sketch-and-lift* (SL) approach solves an SDP on a subsampled dataset and then propagates the solution to all data points by a nearest-centroid rounding procedure. It is shown that the SL approach enjoys a similar exact recovery threshold as the K-means SDP on the full dataset, which is known to be information-theoretically tight under the Gaussian mixture model. The SL method can be made adaptive with enhanced theoretic properties when the cluster sizes are unbalanced. Our simulation experiments demonstrate that the statistical accuracy of the proposed method outperforms state-of-the-art fast clustering algorithms without sacrificing too much computational efficiency, and is comparable to the original K-means SDP with substantially reduced runtime.

# This repository
## Purpose
The codes are considered for one case (p changes under equal size condition). The codes for other conditions are similar. We give the codes for five methods here. (M1: *SL*, M2: *BCSL*, M3: *WSL*, M4: *ME-SL*, M5:  *MR-WSL*.) And we consider comparisons between Kmeans++ and MR-WSL method for two benchmark datasets.
## Steps
 - First you will need to install the optimization toolbox SDPNAL+ from the following website:\
    https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/ \
   Please download and install SDPNAL+ correctly following the instructions in the website. Then put all our files in the path of SDPNAL+ folder. (Or you may just download the file whole.zip, which contains the toolbox.)
 - Next run SDPNALplus_Demo.m in toolbox SDPNAL+ to properly install it.
 - Then we could run m1pc.m - m5pc.m to get the error rate and time cost for GMM models.
 - Finally we could run plot_e.m and plot_t.m to get the corresponding plots.
 - The codes le2015.m and re2016.m are for comparisons to kmeans++ for benchmark data.
## Contents
The files in this repository are:

- kmeans_SDP_1.m: The function we use to get the clusters by SDP method.
- optim.m: The function we will use in our methods for optimal matching.
- m1pc.m - m5pc.m:  The main code for methods M1-M5.
- plot_e.m and  plot_t.m: The codes of plotting error rate and time cost for methods M1-M5.
- kmeansplus.m: The function for kmeans++ method.
- le2015.m and re2016.m: The codes for comparisons of M5 to kmeans++ for benchmark data.
- Levine_32dim.fcs: The dataset from Levine et al. (2015).
- s1.txt and s1cb.txt: The dataset from M. Rezaei and P. Fränti. (2016)
- fca_readfcs.m: The function to read Levine_32dim.fcs data.
- whole.zip: The file contains all the files above and the toolbox SDPNAL+.
## Remarks
One could choose download whole.zip if one would not prefer to download and install the optimization toolbox SDPNAL+ through above steps. In practice, we will run the algorithm for blocks in m4pc.m parallelly with respect to M4. However, here for convenience we run the blocks one by one.
