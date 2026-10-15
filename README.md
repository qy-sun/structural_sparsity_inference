# **Code for "Data-Adaptive Inference with FDR Control for Structural Sparsity Problems"**

Under review

---

This repository contains code for reproducing the results in the paper *"Data-Adaptive Inference with FDR Control for Structural Sparsity Problems"*. 

It implements a unified framework for data-adaptive inference with false discovery rate (FDR) control under structural sparsity, including:

* **TLasso: Transformational Lasso**
* **StatKnock: Stability Transformational Knockoff**
* **SATLasso: Sequential Adjusted Transformational Lasso**
* **Two-Stage "SATLasso + StatKnock" Procedure**

---

## **Description**

The repository is organized into the following main directories:

| Folder | Description |
|:--|:--|
| **function/** | Core implementation of the proposed methods, including TLasso, SATLasso, StatKnock, and SATLasso+StatKnock. |
| **function/other methods/** | Comparative baseline methods such as SplitKnockoff, GKnockoff, DS, and BY for structural sparsity problems. |
| **simulation/** | Simulation scripts for reproducing all experiments in the paper and its supplement. |
| **simulation/helper/** | Graph difference operators construction and corresponding coefficient generation utilities. |
| **real/** | Scripts for ADNI sMRI and Fama-French five-factor real-data analyses. |
| **real/ADNI/lookup/** | Lookup tables for hierarchical brain regions and the corresponding edge endpoints in the brain region graph. |

Each main file begins with a commented header describing its purpose.  

---

## **Simulation Studies**

- `simulation/TLasso_exp.R` — Reproduces **Table 1** and **Figure 4** in the main text.  
- `simulation/SATLasso+StatKnock_exp.R` — Reproduces **Figure 5** in the main text.
- `simulation/StatKnock_exp.R` — Reproduces **Figures S.8–S.11** in the Supplement.
- `simulation/SATLasso_exp.R` — Reproduces **Tables S.3 and S.4** in the Supplement.
- `simulation/StatKnock_ablation.R` — Reproduces **Figure S.1** and **Table S.1** in the Supplement.
- `simulation/TLasso_ablation.R` — Reproduces **Table S.2** in the Supplement.
- `simulation/SATLasso_ablation.R` — Reproduces **Figures S.2 and S.3** in the Supplement.
---

## **Real Data Analysis**

### **ADNI**

The ADNI analysis is based on the **Alzheimer's Disease Neuroimaging Initiative (ADNI)** structural MRI (sMRI) dataset, focusing on the **Level 5 segmentation results** derived from **MRICloud** processing.  
Due to data-sharing and privacy restrictions, the ADNI dataset is **not included** in this repository.  
Below we describe the **data acquisition, preprocessing, and integration workflow** used to obtain the analysis-ready dataset.  
Reproduction of this section requires **individual ADNI data access approval**.

- **ADNI data portal:** [https://adni.loni.usc.edu/data-samples/adni-data/](https://adni.loni.usc.edu/data-samples/adni-data/) (apply through the IDA: [https://ida.loni.usc.edu/collaboration/access/appApply.jsp?project=ADNI](https://ida.loni.usc.edu/collaboration/access/appApply.jsp?project=ADNI))  
- **MRICloud platform:** [https://mricloud.org/](https://mricloud.org/)  
- **MRIStudio platform:** [https://www.mristudio.org/installation.html](https://www.mristudio.org/installation.html)  
- **Reference tutorial:** [MRICloud T1 Tutorial](https://braingps.mricloud.org/docs/tutorials/mricloud.html)

The data preprocessing follows the official **MRICloud–T1 segmentation pipeline** with manual verification in **ROIEditor**. The following outlines the detailed workflow used to obtain the Level 5 regional volume data used in our analysis:

#### **Step 1: Download from ADNI Portal**  
   - Obtain sMRI data in **DICOM (dcm)** format from the ADNI portal. Our analysis uses **T1-weighted MP-RAGE** acquisitions from the **ADNI** study.  
   - Select relevant metadata (e.g., subject ID, age, sex, diagnosis group, and cognitive assessments).  
   - Download and retain the associated `.csv` files that record acquisition information and assessment scores.  
   - Use 3D structural MRI scans for subsequent processing.

#### **Step 2: Convert DICOM to Analyze Format**  
   - Use the conversion tool **`Dcm2Analyze_v3.exe`** (available via MRIStudio) to transform DICOM files into the **Analyze** `.img` / `.hdr` format pair required by MRICloud.  
   - Confirm that file naming conventions remain consistent (RID or subject ID).

#### **Step 3: Confirm Slice Type in ROIEditor**  
   - Open the converted images in **ROIEditor** (MRIStudio).  
   - Verify that the slice type is correct: *Sagittal*, *Axial*, or *Sagittal data converted to Axial*. In this study, the source files were **Sagittal** scans.

#### **Step 4: Perform Segmentation on MRICloud**  
   - Upload the preprocessed files to **MRICloud** for automatic segmentation. 
   - For batch uploads, combine up to **five subjects per .zip file** (no subfolders, and compressed using the system's default "Compress to ZIP" function while avoiding any third-party compression software).  
   - Choose the appropriate **Slice Type** and **Atlas** version. Recommended: `Adult_286labels_10atlases_V5L` (latest version of M2_252).  
   - Submit and wait for processing to complete.  
   - After segmentation, download the result package and extract:  
     - `*_286Labels_corrected_stats.txt` – contains the regional volumes in native subject space (used as feature matrix `X`).  
     - `multilevel_lookup_table.txt` – hierarchical label mapping file (used to construct transformation matrices).
   - The segmentation results can be visualized directly on the **MRICloud** web interface or inspected locally using **ROIEditor**.

#### **Step 5: Match Cognitive Outcomes**  
   - The cognitive outcome variable (response `y`) is the **ADNI_MEM** composite memory score from **"UW - Neuropsych Summary Scores [ADNI1, GO, 2, 3]"** (Crane et al., 2012). This variable serves as a standardized cognitive performance indicator widely used in ADNI literature.  
   - Match each MRI sample to `ADNI_MEM` by aligning:
     - The last four digits of the file name with the **RID** field.  
     - The **VISCODE2** visit code of the cognitive assessment (`m06`).

#### **Step 6: Integrate and Construct Transformation Matrices**  
   - Combine `*_286Labels_corrected_stats.txt` (regional measures) with the matched `ADNI_MEM` scores into a single dataset. Of the 286 atlas labels, the **279** brain structures with a full hierarchical assignment are retained as covariates.  
   - Each regional volume is centered across subjects and rescaled to Euclidean norm `sqrt(n)`, where `n` is the number of subjects, and the response is centered. The volumes are used in absolute form.  
   - Construct transformation matrices (`D1`, `D2`, `D3`) using the **lookup table**, which define the structural sparsity constraints used by the proposed SATLasso–StatKnock framework. Two Level-5 regions are treated as connected when their Level-3 labels agree.

### **FF5**

The FF5 analysis uses the Banks industry portfolio and the Fama-French five factors from the public Kenneth French Data Library. From `real/FF5/`, run `Real Data Preprocess_FF5.R` to prepare `financial_banks_ff5_model_data.csv`, then run `Real Data FF5.R`. The default analysis uses quarterly blocks and 200 iterations; outputs include `frequency_FF5.csv`, `TLasso_Fits_FF5.Rdata`, and `lookup/FF5_row_labels.csv`.

---

## **Additional Notes**

Some experiments were executed on the Duke Compute Cluster (DCC) using parallel jobs to accelerate replication. Minor numerical discrepancies may occur when running locally due to randomization, differences in numerical solvers, and parallel execution order.
