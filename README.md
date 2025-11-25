# Causal Inference: Partisan Sorting and Affective Polarization  
### ANES 2016–2020 Panel • Target Trial Emulation • AIPW • DML

<!-- Badges -->
<p align="left">
  <img src="https://img.shields.io/badge/Code-R%20%7C%20Python-blue.svg" alt="Code: R | Python">
  <img src="https://img.shields.io/badge/Field-Political%20Science%20%7C%20Causal%20Inference-brightgreen.svg" alt="Field: Political Science | Causal Inference">
  <img src="https://img.shields.io/badge/Methods-TTE%20%7C%20AIPW%20%7C%20DML%20%7C%20G--Formula-orange.svg" alt="Methods: TTE | AIPW | DML | G-Formula">
  <img src="https://img.shields.io/badge/Data-ANES%202016–2020%20Panel-critical.svg" alt="ANES Panel">
  <img src="https://img.shields.io/badge/License-Restricted%20(Data%20Not%20Included)-lightgrey.svg" alt="License: Restricted">
  <img src="https://img.shields.io/badge/Status-Work%20in%20Progress-yellow.svg" alt="WIP">
</p>

---

This repository contains replication materials for two connected causal inference studies examining  
**how partisan sorting shapes affective polarization** in American politics.  
Both projects apply **Target Trial Emulation (TTE)** to the **ANES 2016–2020 Panel**, using  
**Augmented Inverse Probability Weighting (AIPW)** and **Debiased Machine Learning (DML)**.  
A separate working paper additionally incorporates a **sequential g-formula mediation analysis**.

The repository includes:

- **Data Preprocessing.ipynb** — Variable construction & preprocessing  
- **IPSA Causal AIPW.R** — AIPW code for IPSA conference version  
- **IPSA Causal DML.ipynb** — DML code for IPSA conference version  
- **IPSA Conference PDF** — Full IPSA paper  
- **Working Paper AIPW.R** — Updated AIPW workflow  
- **Working Paper DML.ipynb** — Updated DML workflow  
- **README.md** — This file  

---

## 🔹 1. Data Availability and Copyright Notice  

All analyses use the **ANES 2016–2020 Panel Study**.  
Because ANES data are **copyrighted**, the original dataset **cannot be redistributed** here.

➡️ To replicate results, download ANES directly from:  
https://electionstudies.org/data-center/panel-studies/

All scripts assume the original ANES variable names and structure.

---

## 🔹 2. AI Assistance Disclosure  

Some preprocessing and modeling were generated with assistance from  
**AI-based code generation tools **.  
AI assistance was used **only for coding support**, not for producing or altering the underlying dataset.
All methodological decisions, TTE design choices, modeling assumptions, and interpretations are entirely my own.

---

## 🔹 3. ⚠️ Critical Methodological Warning (IPSA Version)

The **IPSA Conference** code and paper include a **serious causal identification error**:

> The IPSA version incorrectly adjusts for **Affective Polarization (2016)** when estimating  
> the effect of **Partisan Sorting (2016)** on **Affective Polarization (2020)**.

This introduces an **after-treatment variable** into the confounder set.

### Why is this a critical error?

- Violates **temporal ordering**  
- Violates **exchangeability**  
- Opens **collider bias**  
- Breaks **back-door criterion**  
- Invalidates the TTE structure  
- Produces **biased causal estimates**

The IPSA version is included **only for transparency**.  
It should **NOT** be treated as a valid causal estimate.

---

## 🔹 4. About the Working Paper Version  

The Working Paper builds on the IPSA project but corrects all major issues:

✔ Correct pre-treatment confounder set  
✔ Strict **time-zero alignment** (ANES 2016 pre-election wave)  
✔ Proper **ITT estimand** under TTE  
✔ Doubly robust AIPW + DML estimators  
✔ **Sequential g-formula** mediation for direct/indirect effects  
✔ No after-treatment conditioning  

Because the manuscript is not yet finalized, only the replication code is included.  
The PDF is intentionally not public.

---

## 🔹 5. Requesting the Full Working Paper  

If you wish to read the updated paper:

### **“From Identity to Emotion: Causal Evidence on How Partisan Sorting Generates and Sustains Affective Polarization”**

please contact:

📧 **jihun9965@gmail.com**

---

## 🔹 6. Repository Structure  

```

.
├── Data Preprocessing.ipynb
├── IPSA Causal AIPW.R
├── IPSA Causal DML.ipynb
├── IPSA_Leveraging Target Trial Emulation...pdf
├── Working Paper AIPW.R
├── Working Paper DML.ipynb
└── README.md

```

---

## 🔹 7. Citation  

If you reference this repository, please cite:

Kang, Ji Hun. *Leveraging Target Trial Emulation in Political Science: Assessing the Causal Effect of Partisan Sorting on Affective Polarization.*  
Presented at the 2024 International Political Science Association World Congress.

---

## 🔹 8. Contact  

For questions or working paper access:  
📧 **jihun9965@gmail.com**

---


