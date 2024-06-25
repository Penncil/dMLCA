# dMLCA
Distributed Multi-site Latent Class Analysis

# Outline
1. Description
2. dMLCA illustration

# Description
The latent class analysis (LCA) model is a widely recognized and effective method for identifying disease subphenotypes. The emergence of clinical research networks enables multi-site data integration, providing more generalizable and reliable analysis. However, implementing LCA in a multi-site setting is challenging due to the prohibition of patient-level data sharing and the heterogeneity among different sites. The traditional divide-and-conquer strategy (e.g., meta-analysis) cannot be applied because the number and characteristics of latent classes learned on each site may be different and thus hard to combine. We develop a distributed multi-site latent class analysis (dMLCA) approach. dMLCA allows multiple sites to learn a common set of latent classes collaboratively without sharing patient-level data and can effectively account for the heterogeneity among different sites. We showed the validity and strength of the dMLCA through simulation. We applied it to identify subphenotypes of multisystem inflammatory syndrome in children (MIS-C) using data from 864 patients from nine children's hospitals.

# dMLCA illustration
<img width="835" alt="Screenshot 2024-06-24 at 23 52 11" src="https://github.com/Penncil/dMLCA/assets/70713739/81fa2e59-43b1-4c24-95c8-069db90e0a35">

