#!/usr/bin/env python3
"""
08_merge_annotate.py  —  v3

Merges all annotation layers into a single master table and exports
per-functional-group gene lists for Seurat/CellChat.


# Module-level default — overridden in main() per species
'gene_id' = "gene_id"

Key improvements over v2:
  - Fixed eggNOG column parsing (GO terms, COG cat, KEGG pathway, PFAMs)
  - Removed empty columns (tigrfam, ips_signalp, ips_tmhmm, eggnog_mouse_gene)
  - Fixed SignalP/IPS column collision (IPS no longer creates signalp/tmhmm cols)
  - Added DeepTMHMM topology (--tmhmm)
  - Added --proteome_fasta to seed full gene universe
  - DEAD_+Helicase_C → Stem_Cell_GMP/Vasa (not Antiviral)
  - Longest-key-first TIER3_NAME_TO_GROUP lookup
  - New protein group: AMP_Venom (antimicrobial peptides / venom)
  - New columns: all_go_terms, kegg_pathway_names, cog_description,
                 best_human_gene_name, protein_length_aa
  - Expanded TIER3_NAME_TO_GROUP covering 179 HMM profiles
  - Expanded EGGNOG_GENE_KEYWORDS covering all major gene families
  - assign_cellchat_role uses DeepTMHMM TM count + RECEPTOR_GROUPS

Usage:
    python3 08_merge_annotate.py \\
        --species         Gfas \\
        --run_dir         /scratch/.../runs/Gfas_20260605 \\
        --proteome_fasta  /nethome/.../gfas_1.0.proteins.fasta \\
        --eggnog          .../01_eggnog/Gfas.emapper.annotations \\
        --orthofinder     .../02_orthofinder/Results_OrthoFinder/Results_Jun05 \\
        --interproscan    .../03_interproscan/Gfas_interproscan.tsv \\
        --rbh             .../04_rbh_swissprot/Gfas_rbh_hits.tsv \\
        --signalp         .../05_signalp/prediction_results.txt \\
        --tmhmm           .../05_tmhmm/full/results/TMRs.gff3 \\
        --tier3           .../06_tier3_hmmsearch/Gfas_tier3_domains.domtblout \\
        --tier4           .../07_tier4_blast/Gfas_tier4_blast.tsv \\
        --out_table       .../08_final/Gfas_master_annotation.tsv \\
        --out_genelists   .../08_final/gene_lists/
"""

import argparse
import re
import urllib.request
import time
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
from tqdm import tqdm

# ============================================================
# COG CATEGORY DESCRIPTIONS
# ============================================================
COG_DESCRIPTIONS = {
    'A': 'RNA processing and modification',
    'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation, ribosomal structure and biogenesis',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'O': 'Post-translational modification, protein turnover, chaperones',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolite biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'T': 'Signal transduction mechanisms',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'X': 'Mobilome: prophages, transposons',
    'Z': 'Cytoskeleton',
}

# ============================================================
# TIER 3 HMM PROFILE → PROTEIN GROUP (Pfam accession-based)
# ============================================================
TIER3_TO_GROUP = {
    # ── Innate Immunity: PRR ──────────────────────────────────
    "PF01582": ("Innate_Immunity_PRR",        "TLR_TIR"),
    "PF00560": ("Innate_Immunity_PRR",        "TLR_LRR"),
    "PF13676": ("Innate_Immunity_PRR",        "TLR_LRR_extra"),
    "PF05729": ("Innate_Immunity_PRR",        "NLR_NACHT"),
    "PF13855": ("Innate_Immunity_PRR",        "NLR_LRR"),
    "PF29255": ("Innate_Immunity_PRR",        "NLR_new"),
    "PF00530": ("Innate_Immunity_PRR",        "SRCR"),
    "PF00059": ("Lectin_Receptor",            "CTLD"),
    "PF12662": ("Lectin_Receptor",            "CTLR_receptor"),
    "PF00337": ("Lectin_Receptor",            "Galectin"),
    "PF01476": ("Lectin_Receptor",            "Ficolin"),
    "PF04505": ("Lectin_Receptor",            "Siglec"),
    "PF15009": ("Innate_Immunity_Signalling", "STING"),
    "PF00354": ("Innate_Immunity_PRR",        "Pentraxin_CRP"),
    "PF03061": ("Innate_Immunity_PRR",        "Complement_C3_thioester"),
    "PF07703": ("Innate_Immunity_PRR",        "Alpha2_macroglobulin"),
    "PF01841": ("Innate_Immunity_PRR",        "Complement_anaphylatoxin"),
    "PF01130": ("Innate_Immunity_PRR",        "CD36_scavenger"),
    # ── Innate Immunity: Signalling ───────────────────────────
    "PF02022": ("Innate_Immunity_Signalling", "TRIM_antiviral"),
    "PF02176": ("Innate_Immunity_Signalling", "TRAF_domain"),
    "PF00619": ("Innate_Immunity_Signalling", "CARD_domain"),
    "PF00605": ("Innate_Immunity_Signalling", "IRF_TF"),
    "PF13920": ("Innate_Immunity_Signalling", "RING_E3_ubiquitin"),
    "PF00498": ("Innate_Immunity_Signalling", "F_box_E3"),
    "PF00240": ("Innate_Immunity_Signalling", "Ubiquitin"),
    "PF00627": ("Innate_Immunity_Signalling", "UBA_ubiquitin"),
    "PF00443": ("Innate_Immunity_Signalling", "UCH_deubiquitinase"),
    "PF10584": ("Innate_Immunity_Signalling", "Proteasome_alpha"),
    "PF00227": ("Innate_Immunity_Signalling", "Proteasome_beta"),
    "PF04182": ("Innate_Immunity_Signalling", "Beclin_autophagy"),
    "PF02991": ("Innate_Immunity_Signalling", "ATG8_LC3"),
    "PF10351": ("Innate_Immunity_Signalling", "p62_autophagy"),
    "PF03166": ("Innate_Immunity_Signalling", "SMAD_MH1"),
    "PF03105": ("Innate_Immunity_Signalling", "SMAD_MH2"),
    "PF00019": ("Innate_Immunity_Signalling", "TGFb_ligand"),
    "PF00732": ("Innate_Immunity_Signalling", "TNF_ligand"),
    "PF00020": ("Innate_Immunity_PRR",        "TNFR"),
    "PF01335": ("Innate_Immunity_Effector",   "Death_domain_FADD"),
    # ── Innate Immunity: Effectors ────────────────────────────
    "PF00656": ("Innate_Immunity_Effector",   "Caspase"),
    "PF01823": ("Innate_Immunity_Effector",   "MACPF_Perforin"),
    "PF00275": ("Innate_Immunity_Effector",   "Endonuclease"),
    "PF00264": ("Innate_Immunity_Effector",   "Tyrosinase_melanin"),
    "PF00089": ("Innate_Immunity_Effector",   "Serine_protease"),
    "PF00080": ("Innate_Immunity_Effector",   "SOD_CuZn"),
    "PF00081": ("Innate_Immunity_Effector",   "SOD_MnFe"),
    "PF00199": ("Innate_Immunity_Effector",   "Catalase"),
    "PF00255": ("Innate_Immunity_Effector",   "Glutathione_peroxidase"),
    "PF02798": ("Innate_Immunity_Effector",   "GST_N"),
    "PF00043": ("Innate_Immunity_Effector",   "GST_C"),
    "PF00578": ("Innate_Immunity_Effector",   "Peroxiredoxin"),
    "PF00085": ("Innate_Immunity_Effector",   "Thioredoxin"),
    "PF02898": ("Innate_Immunity_Effector",   "NOS_oxygenase"),
    "PF00899": ("Innate_Immunity_Effector",   "NOS_reductase"),
    "PF00112": ("Innate_Immunity_Effector",   "Cathepsin"),
    # ── Innate Immunity: Antiviral ────────────────────────────
    "PF00270": ("Innate_Immunity_Antiviral",  "DEAD_helicase_RIG"),
    "PF04851": ("Innate_Immunity_Antiviral",  "RIG_CARD"),
    "PF18119": ("Innate_Immunity_Antiviral",  "RIG_helicase_C"),
    "PF17406": ("Innate_Immunity_Antiviral",  "OAS_antiviral"),
    "PF10609": ("Innate_Immunity_Antiviral",  "OAS_2_5A"),
    # ── AMP / Venom ───────────────────────────────────────────
    "PF01823": ("AMP_Venom",                  "MACPF_perforin"),   # also Effector
    "PF06758": ("AMP_Venom",                  "Actinoporin"),
    "PF05730": ("AMP_Venom",                  "Cytolysin_RTX"),
    "PF01355": ("AMP_Venom",                  "Hemolysin"),
    "PF06271": ("AMP_Venom",                  "Aerolysin"),
    "PF00301": ("AMP_Venom",                  "Defensin_beta"),
    "PF01097": ("AMP_Venom",                  "Defensin_invertebrate"),
    "PF08030": ("AMP_Venom",                  "Defensin_propeptide"),
    "PF00278": ("AMP_Venom",                  "Saposin_SAPLIP"),
    "PF05764": ("AMP_Venom",                  "Saposin_B"),
    "PF01342": ("AMP_Venom",                  "Saposin_A"),
    "PF00188": ("AMP_Venom",                  "CAP_SCP_CRISP"),
    "PF01403": ("AMP_Venom",                  "Snaclec_venom"),
    "PF07654": ("AMP_Venom",                  "Toxin_3finger"),
    "PF00555": ("AMP_Venom",                  "Huwentoxin"),
    "PF00014": ("AMP_Venom",                  "Kunitz_BPTI"),
    "PF00068": ("AMP_Venom",                  "Phospholipase_A2"),
    "PF09056": ("AMP_Venom",                  "PLA2_inhibitor"),
    "PF00062": ("AMP_Venom",                  "Lysozyme_C"),
    "PF00959": ("AMP_Venom",                  "Lysozyme_G"),
    "PF31277": ("AMP_Venom",                  "Damicornin"),
    "PF00556": ("AMP_Venom",                  "ShK_toxin"),
    "PF06369": ("AMP_Venom",                  "Actinoporin_alt"),
    "PF07365": ("AMP_Venom",                  "Conotoxin"),
    "PF00666": ("AMP_Venom",                  "Cathelicidin"),
    # ── HSP / Chaperones / UPR ───────────────────────────────
    "PF00012": ("HSP_Chaperone",              "HSP70"),
    "PF02976": ("HSP_Chaperone",              "HSP70_substrate"),
    "PF00118": ("HSP_Chaperone",              "HSP60_GroEL"),
    "PF00226": ("HSP_Chaperone",              "HSP40_DnaJ"),
    "PF02518": ("HSP_Chaperone",              "HSP90"),
    "PF13589": ("HSP_Chaperone",              "HSP90_middle"),
    "PF00011": ("HSP_Chaperone",              "HSP20_small"),
    "PF00160": ("HSP_Chaperone",              "Cyclophilin_PPIase"),
    "PF00254": ("HSP_Chaperone",              "FKBP_PPIase"),
    "PF01920": ("HSP_Chaperone",              "Prefoldin"),
    "PF02996": ("HSP_Chaperone",              "Prefoldin_2"),
    "PF01323": ("HSP_Chaperone",              "PDI_thioredoxin"),
    "PF00262": ("HSP_Chaperone",              "Calreticulin_calnexin"),
    # ── Stem Cells: GMP ───────────────────────────────────────
    "PF02171": ("Stem_Cell_GMP",              "Piwi_Argonaute"),
    "PF00271": ("Stem_Cell_GMP",              "Helicase_C_Vasa"),
    "PF03417": ("Stem_Cell_GMP",              "Nanos"),
    "PF05741": ("Stem_Cell_GMP",              "Nanos_zinc"),
    "PF00806": ("Stem_Cell_GMP",              "Pumilio_PUF"),
    "PF00567": ("Stem_Cell_GMP",              "Tudor_domain"),
    "PF08699": ("Stem_Cell_GMP",              "Argonaute_PAZ"),
    # ── Stem Cells: Proliferation ─────────────────────────────
    "PF00705": ("Stem_Cell_Proliferation",    "PCNA"),
    "PF00307": ("Stem_Cell_Proliferation",    "Cyclin_N"),
    "PF02032": ("Stem_Cell_Proliferation",    "Cyclin_C"),
    "PF05195": ("Stem_Cell_Proliferation",    "Aurora_kinase"),
    "PF00514": ("Stem_Cell_Proliferation",    "Armadillo_beta_catenin"),
    # ── Calcification ─────────────────────────────────────────
    "PF00194": ("Calcification_Universal",    "Carbonic_anhydrase_alpha"),
    "PF07836": ("Calcification_Universal",    "Carbonic_anhydrase_beta"),
    "PF24748": ("Calcification_Universal",    "Galaxin"),
    "PF01391": ("Calcification_Universal",    "Collagen_fibrillar"),
    "PF00093": ("Calcification_Universal",    "VWD_matrix"),
    "PF00187": ("Calcification_Universal",    "Chitin_binding"),
    "PF01607": ("Calcification_Universal",    "Chitin_bind_2"),
    "PF00704": ("Calcification_Universal",    "Chitinase"),
    "PF01644": ("Calcification_Universal",    "Chitin_synthase"),
    "PF00090": ("Calcification_Universal",    "Thrombospondin"),
    "PF00822": ("Calcification_IonTransport", "Claudin"),
    "PF07565": ("Calcification_IonTransport", "SLC4_bicarbonate"),
    "PF01699": ("Calcification_IonTransport", "SLC26_bicarbonate"),
    "PF00916": ("Calcification_IonTransport", "SLC26_sulfate"),
    "PF00122": ("Calcification_IonTransport", "Ca_ATPase"),
    "PF00137": ("Calcification_IonTransport", "V_ATPase"),
    "PF00909": ("Calcification_IonTransport", "Ammonium_transport"),
    "PF01052": ("Calcification_IonTransport", "NCX_Na_Ca_exchanger"),
    # ── Symbiosis ─────────────────────────────────────────────
    "PF00071": ("Symbiosis",                  "Rab_GTPase"),
    "PF00566": ("Symbiosis",                  "Rab_GTPase_2"),
    "PF00957": ("Symbiosis",                  "SNARE"),
    "PF05739": ("Symbiosis",                  "SNARE_assoc"),
    "PF01852": ("Symbiosis",                  "START_lipid"),
    "PF02770": ("Symbiosis",                  "AcylCoA_synthetase"),
    "PF12349": ("Symbiosis",                  "NPC_sterol"),
    "PF00019": ("Symbiosis",                  "TGFb_symbiosis"),  # dual role
    # ── ECM / Adhesion ────────────────────────────────────────
    "PF00008": ("ECM_Adhesion",               "EGF_like"),
    "PF00041": ("ECM_Adhesion",               "Fibronectin_III"),
    "PF00053": ("ECM_Adhesion",               "Laminin_N"),
    "PF03160": ("ECM_Adhesion",               "Integrin_alpha"),
    "PF00028": ("ECM_Adhesion",               "Cadherin"),
    "PF05986": ("ECM_Adhesion",               "ADAM_disintegrin"),
    "PF05679": ("ECM_Adhesion",               "Hyaluronan_binding"),
    # ── Transcription Factors ─────────────────────────────────
    "PF00096": ("Transcription_Factor",       "C2H2_zinc_finger"),
    "PF00250": ("Transcription_Factor",       "Forkhead_FOX"),
    "PF00178": ("Transcription_Factor",       "ETS_domain"),
    "PF00105": ("Transcription_Factor",       "Nuclear_receptor_C4"),
    "PF03736": ("Transcription_Factor",       "Gli_zinc_finger"),
    "PF07645": ("Transcription_Factor",       "EGF_Notch"),
    "PF00853": ("Transcription_Factor",       "Runt_RUNX"),
    "PF00870": ("Transcription_Factor",       "p53_family"),
    "PF00319": ("Transcription_Factor",       "MADS_box"),
    "PF00249": ("Transcription_Factor",       "Myb_DNA_bind"),
    "PF00170": ("Transcription_Factor",       "bZIP_AP1"),
    "PF00907": ("Transcription_Factor",       "T_box"),
    "PF00554": ("Transcription_Factor",       "RHD_NF_kB"),
    # ── Kinases ───────────────────────────────────────────────
    "PF00069": ("Kinase",                     "Ser_Thr_kinase"),
    "PF07714": ("Kinase",                     "Tyr_RTK"),
    # ── GPCR ──────────────────────────────────────────────────
    "PF00001": ("GPCR",                       "7tm_class_A"),
    "PF10320": ("GPCR",                       "7tm_class_C"),
    "PF00003": ("GPCR",                       "7tm_odorant"),
    # ── Ion Channels ──────────────────────────────────────────
    "PF00520": ("Ion_Channel",                "Voltage_gated"),
    "PF00060": ("Ion_Channel",                "Ligand_gated"),
    "PF01007": ("Ion_Channel",                "Inward_rectifier_K"),
    "PF04547": ("Ion_Channel",                "CLC_chloride"),
    "PF00625": ("Ion_Channel",                "VG_calcium"),
    # ── Lectin Receptor ───────────────────────────────────────
    "PF02171": ("Lectin_Receptor",            "Piwi"),   # also GMP
    "PF00076": ("Neuron_Neuropeptide",        "ELAV_RRM"),
    # ── Neuron ────────────────────────────────────────────────
    "PF00355": ("Neuron_Neuropeptide",        "ELAV_HuD"),
}

# ============================================================
# TIER 3 HMM PROFILE NAME → PROTEIN GROUP
# (Matches on the NAME field stored in tier3_domain_hits)
# Sorted longest-key-first at lookup time to handle multi-domain keys
# ============================================================
TIER3_NAME_TO_GROUP = {
    # Multi-domain keys first (longest)
    "DEAD_,Helicase_C":        ("Stem_Cell_GMP",              "Vasa_DDX4"),
    "Helicase_C,DEAD_":        ("Stem_Cell_GMP",              "Vasa_DDX4"),
    "Arm_,DEAD_,Helicase_C":   ("Stem_Cell_GMP",              "Vasa_DDX4"),
    "GST_N_,GST_C_":           ("Innate_Immunity_Effector",   "GST"),
    "GST_C_,GST_N_":           ("Innate_Immunity_Effector",   "GST"),
    "Proteasome_,Proteasome_A_N": ("Innate_Immunity_Signalling", "Proteasome"),
    "AhpC-TSA_,Thioredoxin_":  ("Innate_Immunity_Effector",   "Peroxiredoxin_Thioredoxin"),
    # Single-domain keys
    "TIR_":           ("Innate_Immunity_PRR",        "TLR_TIR"),
    "LRR_1":          ("Innate_Immunity_PRR",        "TLR_LRR"),
    "LRR_8":          ("Innate_Immunity_PRR",        "NLR_LRR"),
    "NACHT_":         ("Innate_Immunity_PRR",        "NLR_NACHT"),
    "SRCR_":          ("Innate_Immunity_PRR",        "SRCR"),
    "Lectin_C":       ("Lectin_Receptor",            "CTLD"),
    "Gal-bind_lectin":("Lectin_Receptor",            "Galectin"),
    "LysM_":          ("Lectin_Receptor",            "Ficolin"),
    "F5_F8_type_C":   ("Lectin_Receptor",            "Lectin_carbohydrate"),
    "ig_":            ("Innate_Immunity_PRR",        "Immunoglobulin"),
    "STING_LBD":      ("Innate_Immunity_Signalling", "STING"),
    "TNFR_c6":        ("Innate_Immunity_PRR",        "TNFR"),
    "DED_":           ("Innate_Immunity_Effector",   "Death_domain"),
    "Peptidase_C14":  ("Innate_Immunity_Effector",   "Caspase"),
    "MACPF_":         ("AMP_Venom",                  "MACPF_Perforin"),
    "Tyrosinase_":    ("Innate_Immunity_Effector",   "Tyrosinase_melanin"),
    "GSHPx_":         ("Innate_Immunity_Effector",   "Glutathione_peroxidase"),
    "GST_N_":         ("Innate_Immunity_Effector",   "GST"),
    "GST_C_":         ("Innate_Immunity_Effector",   "GST"),
    "AhpC-TSA":       ("Innate_Immunity_Effector",   "Peroxiredoxin"),
    "Thioredoxin_":   ("Innate_Immunity_Effector",   "Thioredoxin"),
    "Sod_Cu_":        ("Innate_Immunity_Effector",   "SOD_CuZn"),
    "Sod_Fe_N":       ("Innate_Immunity_Effector",   "SOD_MnFe"),
    "Catalase_":      ("Innate_Immunity_Effector",   "Catalase"),
    "Peptidase_C1":   ("Innate_Immunity_Effector",   "Cathepsin"),
    "Peptidase_C6":   ("Innate_Immunity_Effector",   "Cathepsin"),
    "NO_synthase":    ("Innate_Immunity_Signalling", "NOS"),
    "RHD_DNA_bind":   ("Transcription_Factor",       "NF-kB_RHD"),
    "ubiquitin_":     ("Innate_Immunity_Signalling", "Ubiquitin"),
    "UCH_":           ("Innate_Immunity_Signalling", "Deubiquitinase"),
    "UBA_":           ("Innate_Immunity_Signalling", "UBA_ubiquitin"),
    "Proteasome_":    ("Innate_Immunity_Signalling", "Proteasome"),
    "ATG8_":          ("Innate_Immunity_Signalling", "Autophagy_LC3"),
    "DEAD_":          ("Innate_Immunity_Antiviral",  "DEAD_helicase"),  # alone = RIG-I
    "RNA_helicase":   ("Innate_Immunity_Antiviral",  "RNA_helicase"),
    "RRM_1":          ("Stem_Cell_GMP",              "RNA_binding"),
    "Piwi_":          ("Stem_Cell_GMP",              "Piwi_Argonaute"),
    "PCNA_N":         ("Stem_Cell_Proliferation",    "PCNA"),
    "HLH_":           ("Transcription_Factor",       "bHLH"),
    "Homeodomain_":   ("Transcription_Factor",       "Homeodomain"),
    "HMG_box":        ("Transcription_Factor",       "HMG_SOX"),
    "Myb_Cef":        ("Transcription_Factor",       "MYB"),
    "GATA_":          ("Transcription_Factor",       "GATA"),
    "PHD_":           ("Transcription_Factor",       "PHD_chromatin"),
    "Arm_":           ("Transcription_Factor",       "Armadillo_Wnt"),
    "PH_":            ("Signalling",                 "PH_domain"),
    "SH2_":           ("Signalling",                 "SH2_domain"),
    "C2_":            ("Signalling",                 "C2_calcium"),
    "HEAT_":          ("Signalling",                 "HEAT_repeat"),
    "HATPase_c":      ("HSP_Chaperone",              "ATPase_chaperone"),
    "Cpn60_TCP1":     ("HSP_Chaperone",              "HSP60_chaperonin"),
    "DnaJ_":          ("HSP_Chaperone",              "HSP40_DnaJ"),
    "HSP70_":         ("HSP_Chaperone",              "HSP70"),
    "7TM_GPCR_Srsx":  ("GPCR",                      "7tm_class_A"),
    "Ras_":           ("Symbiosis",                  "Ras_GTPase"),
    "CRAL_TRIO_2":    ("Symbiosis",                  "CRAL_TRIO_lipid"),
    "Laminin_G_1":    ("ECM_Adhesion",               "Laminin_G"),
    "PK_Tyr_Ser-Thr": ("Kinase",                    "Protein_kinase"),
    # AMP/Venom domains
    "Trypsin_":       ("Innate_Immunity_Effector",   "Serine_protease"),
    "Kazal_1":        ("AMP_Venom",                  "Kazal_inhibitor"),
    "Pentraxin_":     ("Innate_Immunity_PRR",        "Pentraxin_CRP"),
    "Toxin_3":        ("AMP_Venom",                  "Toxin_3finger"),
    "Kunitz_":        ("AMP_Venom",                  "Kunitz_toxin"),
    "PF01549":        ("AMP_Venom",                  "ShK_toxin"),
    "Phospholipase_A2": ("AMP_Venom",               "Phospholipase_A2"),
    "Lysozyme_C":     ("AMP_Venom",                  "Lysozyme_C"),
    "Lysozyme_G":     ("AMP_Venom",                  "Lysozyme_G"),
    "Saposin_":       ("AMP_Venom",                  "Saposin_SAPLIP"),
    "CAP_":           ("AMP_Venom",                  "CAP_SCP_CRISP"),
    "CAP":            ("AMP_Venom",                  "CAP_SCP_CRISP"),
    "Defensin_":      ("AMP_Venom",                  "Defensin"),
    "Actinoporin":    ("AMP_Venom",                  "Actinoporin"),
    "Hemolysin_":     ("AMP_Venom",                  "Hemolysin"),
    "Galaxin_":       ("Calcification_Universal",    "Galaxin"),
    "Galaxin_repeat": ("Calcification_Universal",    "Galaxin"),
    "Pentaxin":       ("Innate_Immunity_PRR",        "Pentraxin_CRP"),
    "ResIII":         ("Innate_Immunity_Effector",   "Restriction_endonuclease"),
    "Anemone_cytotox":("AMP_Venom",                  "Actinoporin_cytotoxin"),
    "ShK":            ("AMP_Venom",                  "ShK_toxin"),
    "Phospholip_A2_1":("AMP_Venom",                  "Phospholipase_A2"),
    "Phospholip_A2_3":("AMP_Venom",                  "Phospholipase_A2"),
    "Damicornin":     ("AMP_Venom",                  "Damicornin"),
    "Trypsin":        ("Innate_Immunity_Effector",   "Serine_protease"),
    "DEAD":           ("Innate_Immunity_Antiviral",  "DEAD_helicase"),
    "Ion_trans":      ("Ion_Channel",                "Ion_transporter"),
    "Lig_chan":        ("Ion_Channel",                "Ligand_gated_channel"),
    "IRK":            ("Ion_Channel",                "Inward_rectifier_K"),
    "Anoctamin":      ("Ion_Channel",                "Anoctamin_TMEM16"),
    "Ank_2":          ("Signalling",                 "Ankyrin_repeat"),
    "EGF_CA":         ("ECM_Adhesion",               "EGF_calcium_binding"),
    "cEGF":           ("ECM_Adhesion",               "EGF_calcium_binding"),
    "fn3":            ("ECM_Adhesion",               "Fibronectin_III"),
    "TSP_1":          ("ECM_Adhesion",               "Thrombospondin_1"),
    "Collagen":       ("Calcification_Universal",    "Collagen"),
    "Cadherin":       ("ECM_Adhesion",               "Cadherin"),
    "Laminin_G_1":    ("ECM_Adhesion",               "Laminin_G"),
    "zf-C3HC4_3":     ("Innate_Immunity_Signalling", "RING_E3_ubiquitin"),
    "zf-TRAF":        ("Innate_Immunity_Signalling", "TRAF_zinc_finger"),
    "NHL":            ("Innate_Immunity_Antiviral",  "TRIM_NHL_antiviral"),
    "zf-C2H2":        ("Transcription_Factor",       "C2H2_zinc_finger"),
    "zf-C4":          ("Transcription_Factor",       "Nuclear_receptor_C4"),
    "Ets":            ("Transcription_Factor",       "ETS_domain"),
    "bZIP_1":         ("Transcription_Factor",       "bZIP_AP1"),
    "SRF-TF":         ("Transcription_Factor",       "MADS_SRF"),
    "IRF":            ("Innate_Immunity_Signalling", "IRF_TF"),
    "CARD":           ("Innate_Immunity_Signalling", "CARD_domain"),
    "DED":            ("Innate_Immunity_Effector",   "Death_domain"),
    "T-box":          ("Transcription_Factor",       "T_box"),
    "Forkhead":       ("Transcription_Factor",       "Forkhead_FOX"),
    "Runt":           ("Transcription_Factor",       "Runt_RUNX"),
    "P53":            ("Transcription_Factor",       "p53_family"),
    "Myb_DNA-binding":("Transcription_Factor",       "MYB"),
    "HMG_box":        ("Transcription_Factor",       "HMG_SOX"),
    "GATA":           ("Transcription_Factor",       "GATA"),
    "PHD":            ("Transcription_Factor",       "PHD_chromatin"),
    "bZIP_1":         ("Transcription_Factor",       "bZIP_AP1"),
    "Notch":          ("Transcription_Factor",       "Notch_TF"),
    "RIG-I_C":        ("Innate_Immunity_Antiviral",  "RIG_I_helicase_C"),
    "STING_LBD":      ("Innate_Immunity_Signalling", "STING"),
    "TIR_2":          ("Innate_Immunity_PRR",        "TLR_TIR"),
    "TNFR_c6":        ("Innate_Immunity_PRR",        "TNFR"),
    "LRR_NLRC3":      ("Innate_Immunity_PRR",        "NLR_LRR"),
    "MH2":            ("Innate_Immunity_Signalling", "SMAD_MH2"),
    "TGF_beta":       ("Calcification_Signalling",   "TGFb_BMP_ligand"),
    "SNARE":          ("Symbiosis",                  "SNARE_vesicle"),
    "Synaptobrevin":  ("Symbiosis",                  "SNARE_vesicle"),
    "START":          ("Symbiosis",                  "START_lipid_transfer"),
    "CRAL_TRIO_2":    ("Symbiosis",                  "CRAL_TRIO_lipid"),
    "RabGAP-TBC":     ("Symbiosis",                  "Rab_GAP"),
    "Ras":            ("Symbiosis",                  "Ras_GTPase"),
    "Sterol-sensing": ("Symbiosis",                  "Sterol_sensing"),
    "CHGN":           ("Symbiosis",                  "Chitin_glucan"),
    "Ammonium_transp":("Calcification_IonTransport", "Ammonium_transport"),
    "Sulfate_transp": ("Calcification_IonTransport", "Sulfate_transporter"),
    "Carb_anhydrase": ("Calcification_Universal",    "Carbonic_anhydrase"),
    "CBM_14":         ("Calcification_Universal",    "Chitin_binding"),
    "PMP22_Claudin":  ("Calcification_IonTransport", "Claudin"),
    "CD36":           ("Innate_Immunity_PRR",        "CD36_scavenger"),
    "F5_F8_type_C":   ("Lectin_Receptor",            "Lectin_carbohydrate"),
    "SRCR":           ("Innate_Immunity_PRR",        "SRCR"),
    "Lectin_C":       ("Lectin_Receptor",            "CTLD"),
    "LRR_8":          ("Innate_Immunity_PRR",        "NLR_LRR"),
    "NACHT":          ("Innate_Immunity_PRR",        "NLR_NACHT"),
    "TUDOR":          ("Stem_Cell_GMP",              "Tudor_domain"),
    "Piwi":           ("Stem_Cell_GMP",              "Piwi_Argonaute"),
    "ArgoL1":         ("Stem_Cell_GMP",              "Argonaute"),
    "zf-nanos":       ("Stem_Cell_GMP",              "Nanos_zinc"),
    "PUF":            ("Stem_Cell_GMP",              "Pumilio_PUF"),
    "PCNA_N":         ("Stem_Cell_Proliferation",    "PCNA"),
    "HSP20":          ("HSP_Chaperone",              "HSP20_small"),
    "HSP70":          ("HSP_Chaperone",              "HSP70"),
    "Cpn60_TCP1":     ("HSP_Chaperone",              "HSP60_chaperonin"),
    "DnaJ":           ("HSP_Chaperone",              "HSP40_DnaJ"),
    "Pro_isomerase":  ("HSP_Chaperone",              "Cyclophilin_PPIase"),
    "FKBP_C":         ("HSP_Chaperone",              "FKBP_PPIase"),
    "Calreticulin":   ("HSP_Chaperone",              "Calreticulin_calnexin"),
    "Prefoldin":      ("HSP_Chaperone",              "Prefoldin"),
    "Prefoldin_2":    ("HSP_Chaperone",              "Prefoldin_2"),
    "Proteasome":     ("Innate_Immunity_Signalling", "Proteasome"),
    "Proteasome_A_N": ("Innate_Immunity_Signalling", "Proteasome"),
    "ATG8":           ("Innate_Immunity_Signalling", "Autophagy_LC3"),
    "AhpC-TSA":       ("Innate_Immunity_Effector",   "Peroxiredoxin"),
    "GSHPx":          ("Innate_Immunity_Effector",   "Glutathione_peroxidase"),
    "GST_C":          ("Innate_Immunity_Effector",   "GST"),
    "GST_N":          ("Innate_Immunity_Effector",   "GST"),
    "Catalase":       ("Innate_Immunity_Effector",   "Catalase"),
    "Tyrosinase":     ("Innate_Immunity_Effector",   "Tyrosinase_melanin"),
    "MACPF":          ("AMP_Venom",                  "MACPF_Perforin"),
    "Peptidase_C1":   ("Innate_Immunity_Effector",   "Cathepsin"),
    "Peptidase_C14":  ("Innate_Immunity_Effector",   "Caspase"),
    "MATH":           ("Innate_Immunity_Signalling", "TRAF_MATH"),
    "ubiquitin":      ("Innate_Immunity_Signalling", "Ubiquitin"),
    "UBA":            ("Innate_Immunity_Signalling", "UBA_ubiquitin"),
    "UCH":            ("Innate_Immunity_Signalling", "Deubiquitinase"),
    "ThiF":           ("Innate_Immunity_Signalling", "E1_ubiquitin_activating"),
    "Frizzled":       ("Calcification_Signalling",   "Frizzled_Wnt_receptor"),
    "NO_synthase":    ("Innate_Immunity_Signalling", "NOS"),
    "SH2":            ("Signalling",                 "SH2_domain"),
    "PH":             ("Signalling",                 "PH_domain"),
    "C2":             ("Signalling",                 "C2_calcium"),
    "HEAT":           ("Signalling",                 "HEAT_repeat"),
    "Helicase_C":     ("Signalling",                 "Helicase_C_general"),
    "RRM_1":          ("Neuron_Neuropeptide",         "RNA_binding_neuronal"),
    "Lys":            ("AMP_Venom",                  "Lysozyme"),
    "Sema":           ("ECM_Adhesion",               "Semaphorin"),
}

# ============================================================
# TIER 4 SOURCE → PROTEIN GROUP
# ============================================================
TIER4_TO_GROUP = {
    "SOMP_Calc":   ("Calcification_Universal",  "SOMP_Galaxin"),
    "SCRiP_Venom": ("AMP_Venom",                "SCRiP"),
    "Neuropep_NP": ("Neuron_Neuropeptide",       "Neuropeptide_prepro"),
    "Symbiosis":   ("Symbiosis",                "LePin_lectin"),
    "Cnidocyte":   ("AMP_Venom",                "Minicollagen"),
}

# ============================================================
# EGGNOG GENE NAME KEYWORDS → PROTEIN GROUP
# ============================================================
EGGNOG_GENE_KEYWORDS = {
    # ── HSP / UPR ─────────────────────────────────────────────
    r'^HSP\d+|^HSPA|^HSPB|^HSPC|^HSPD|^HSPE|^HSPH|^DNAJ|^DNAK': ("HSP_Chaperone", "HSP_family"),
    r'^GRP\d+|^HSPA5$|^HSPA1|^CANX$|^CALR$|^ERP\d+|^PDIA': ("HSP_Chaperone", "ER_chaperone"),
    r'^IRE1$|^ERN1$|^EIF2AK3$|^PERK$|^ATF6$|^XBP1$': ("HSP_Chaperone", "UPR_sensor"),
    # ── Transcription Factors ─────────────────────────────────
    r'^SOX\d+|^PAX\d+|^GATA\d+|^RUNX\d+|^ATF\d+|^NFAT|^IRF\d+|^NFKB|^RELA|^RELB|^CREL': ("Transcription_Factor", "TF_immune"),
    r'^TWIST|^MSX\d+|^SP7|^RUNX2|^DLX|^SP3$|^SP1$|^KLF': ("Transcription_Factor", "TF_developmental"),
    r'^POU|^MYC$|^MYCN|^NANOG': ("Transcription_Factor", "TF_pluripotency"),
    r'^FOX[A-Z]|^FOXA|^FOXO|^FOXP|^FOXN': ("Transcription_Factor", "Forkhead_FOX"),
    r'^ETS\d*|^ETV|^ELF|^FLI1|^ERG$': ("Transcription_Factor", "ETS_family"),
    r'^NR\d+|^RORA|^RORB|^PPARA|^PPARB|^PPARG|^LXR|^FXR': ("Transcription_Factor", "Nuclear_receptor"),
    r'^GLI\d+|^SUFU$': ("Transcription_Factor", "Gli_Hedgehog"),
    r'^NOTCH\d+|^HES\d+|^HEY\d+': ("Transcription_Factor", "Notch_TF"),
    r'^JUN$|^JUNB|^JUND|^FOS$|^FOSB|^CREB\d*': ("Transcription_Factor", "bZIP_AP1"),
    r'^TBX\d+|^BRACHYURY|^TBXT$': ("Transcription_Factor", "T_box"),
    r'^SMAD\d+': ("Transcription_Factor", "SMAD"),
    r'^TP53$|^TP63$|^TP73$': ("Transcription_Factor", "p53_family"),
    r'^MYB$|^MYBL\d+': ("Transcription_Factor", "MYB"),
    r'^MEF2|^SRF$': ("Transcription_Factor", "MADS"),
    # ── Immunity ──────────────────────────────────────────────
    r'^TLR\d+|^MYD88$|^IRAK\d+|^TRAF\d+|^STING\d*|^CGAS$|^TMEM173$': ("Innate_Immunity_Signalling", "TLR_signalling"),
    r'^CASP\d+|^BCL2$|^BCLXL|^BAX$|^MCL1$|^BID$|^BAD$': ("Innate_Immunity_Effector", "Apoptosis"),
    r'^DDX58$|^IFIH1$|^DHX58$|^RIG': ("Innate_Immunity_Antiviral", "RIG_I_like"),
    r'^TRIM\d+': ("Innate_Immunity_Antiviral", "TRIM_antiviral"),
    r'^OAS\d+': ("Innate_Immunity_Antiviral", "OAS_antiviral"),
    r'^NLRP\d+|^NOD\d+|^NLRC\d+': ("Innate_Immunity_PRR", "NLR"),
    r'^COMP\d*|^C3$|^C4[AB]|^C5$|^CFB$|^MASP\d+': ("Innate_Immunity_PRR", "Complement"),
    # ── AMP / Venom ───────────────────────────────────────────
    r'^DEFB\d+|^DEFA\d+|^DEFM\d+': ("AMP_Venom", "Defensin"),
    r'^CAMP$|^CATHL\d+': ("AMP_Venom", "Cathelicidin"),
    r'^PLA2G\d+|^SPLA2': ("AMP_Venom", "Phospholipase_A2"),
    r'^LYZ$|^LYSC$|^LYSB': ("AMP_Venom", "Lysozyme"),
    # ── Stem cells ────────────────────────────────────────────
    r'^PIWIL|^PIWI': ("Stem_Cell_GMP", "Piwi"),
    r'^DDX4$|^VASA$': ("Stem_Cell_GMP", "Vasa"),
    r'^NANOS\d': ("Stem_Cell_GMP", "Nanos"),
    r'^TDRD\d+': ("Stem_Cell_GMP", "Tudor"),
    r'^CDC42$|^CBX\d+|^MEIS\d+|^RUNX1$': ("Stem_Cell_HSC", "HSC_self_renewal"),
    r'^MKI67$|^TOP2A$|^PCNA$|^AURKB$|^AURKA$|^CDK1$|^CCNB': ("Stem_Cell_Proliferation", "Proliferation_markers"),
    r'^CTNNB1$|^FZD\d+': ("Calcification_Signalling", "Wnt_beta_catenin"),
    r'^WNT\d+': ("Calcification_Signalling", "Wnt_ligand"),
    # ── Calcification ─────────────────────────────────────────
    r'^CA\d+$|^CAR\d+$|^STPCA': ("Calcification_Universal", "Carbonic_anhydrase"),
    r'^SLC4A|^SLC26A': ("Calcification_IonTransport", "Bicarbonate_transporter"),
    r'^ATP2A|^ATP2B': ("Calcification_IonTransport", "Ca_ATPase"),
    r'^BMP\d+|^BMPR': ("Calcification_Signalling", "BMP_pathway"),
    r'^MMP\d+|^ADAMTS': ("Calcification_Universal", "Matrix_metalloprotease"),
    # ── Symbiosis ─────────────────────────────────────────────
    r'^RAB\d+': ("Symbiosis", "Rab_GTPase"),
    r'^NPC\d+|^STARD\d+|^STAR$': ("Symbiosis", "Sterol_lipid_transport"),
    r'^LGALS\d+': ("Lectin_Receptor", "Galectin"),
    r'^ELOVL\d+|^ACSL\d+|^FASN$': ("Symbiosis", "Lipid_metabolism"),
    r'^RHBG$|^RHCG$|^SLC42A': ("Symbiosis", "Ammonium_transport"),
    # ── Neuron ────────────────────────────────────────────────
    r'^SYT\d+|^SNAP\d+|^SYN\d+|^STX\d+|^VAMP\d+': ("Neuron_Neuropeptide", "Synaptic_vesicle"),
    r'^ELAVL\d+|^HUD$|^HUC$|^HUB$': ("Neuron_Neuropeptide", "ELAV_neuronal"),
    r'^SCN\d+[A-Z]|^KCNA|^KCNB|^CACNA': ("Ion_Channel", "Voltage_gated_neural"),
    r'^GRIA\d+|^GRIN\d+|^GRIP\d+': ("Ion_Channel", "Glutamate_receptor"),
    r'^GABRA|^GABRB|^GABRG': ("Ion_Channel", "GABA_receptor"),
    r'^CHRNA|^CHRNB': ("Ion_Channel", "AChR"),
    r'^HTR\d+|^DRD\d+|^ADRA|^ADRB': ("GPCR", "Monoamine_GPCR"),
    # ── GPCR ──────────────────────────────────────────────────
    r'^ADRA|^ADRB|^HTR\d+|^DRD\d+|^OPRM|^OPRD|^OPRK': ("GPCR", "Monoamine_GPCR"),
    r'^PTHR\d*|^PTH1R|^PTH2R|^CALCR': ("GPCR", "Peptide_GPCR"),
    # ── Kinases ───────────────────────────────────────────────
    r'^EGFR$|^FGFR\d+|^VEGFR|^KIT$|^FLT\d+|^PDGFR': ("Kinase", "RTK"),
    r'^MAPK\d+|^MAP2K|^MAP3K|^ERK\d+|^JNK\d+|^P38|^MAPK8': ("Kinase", "MAPK_cascade"),
    r'^PIK3|^AKT\d+|^MTOR$|^PTEN$': ("Kinase", "PI3K_AKT"),
    # ── ECM ───────────────────────────────────────────────────
    r'^COL\d+A|^LAMA|^LAMB|^LAMC|^FN1$|^FBN\d+': ("ECM_Adhesion", "ECM_structural"),
    r'^ITGA|^ITGB|^CDH\d+|^VCAM|^ICAM': ("ECM_Adhesion", "Cell_adhesion"),
}


def match_gene_keyword(gene_name):
    """Return (protein_group, subgroup) if gene name matches keyword patterns."""
    if not gene_name or pd.isna(gene_name):
        return None, None
    # Clean up: strip _HUMAN/_MOUSE suffix, take first if multiple
    name = str(gene_name).split(',')[0].strip()
    name = re.sub(r'_HUMAN$|_MOUSE$|_RAT$|_BOVIN$|_XENLA$|_DANRE$', '', name)
    for pattern, (group, subgroup) in EGGNOG_GENE_KEYWORDS.items():
        if re.search(pattern, name, re.IGNORECASE):
            return group, subgroup
    return None, None


# ============================================================
# KEGG PATHWAY NAME LOOKUP
# ============================================================
def fetch_kegg_pathway_names(ko_series, max_unique=2000):
    """
    Map KEGG KO numbers to pathway names via KEGG REST API.
    Returns a dict: ko_number → comma-separated pathway names.
    Only fetches unique KO numbers. Skips if too many unique KOs.
    """
    print("  Fetching KEGG pathway names...")
    # Collect unique KO numbers
    all_kos = set()
    for val in ko_series.dropna():
        for ko in str(val).split(','):
            ko = ko.strip()
            if ko.startswith('ko:') or ko.startswith('K'):
                ko_clean = ko.replace('ko:', '').strip()
                if ko_clean:
                    all_kos.add(ko_clean)

    if len(all_kos) > max_unique:
        print(f"    {len(all_kos)} unique KOs — using KEGG batch API (10 per call)...")
        ko_list = sorted(all_kos)
        ko_to_pathway = {}
        failed = 0
        for i in range(0, len(ko_list), 10):
            batch = ko_list[i:i+10]
            url = "https://rest.kegg.jp/get/" + "+".join(f"ko:{k}" for k in batch)
            try:
                with urllib.request.urlopen(url, timeout=15) as resp:
                    text = resp.read().decode('utf-8')
                current_ko = None
                pathways = []
                in_pathway = False
                for line in text.splitlines():
                    if line.startswith('ENTRY'):
                        if current_ko and pathways:
                            ko_to_pathway[current_ko] = '; '.join(pathways[:5])
                        current_ko = line.split()[1].strip()
                        pathways = []
                        in_pathway = False
                    elif line.startswith('PATHWAY'):
                        in_pathway = True
                        parts = line.split(None, 2)
                        if len(parts) >= 3:
                            pathways.append(parts[2].strip())
                    elif in_pathway and line.startswith(' '):
                        parts = line.split(None, 2)
                        if len(parts) >= 2:
                            pathways.append(parts[-1].strip())
                    elif in_pathway and not line.startswith(' '):
                        in_pathway = False
                if current_ko and pathways:
                    ko_to_pathway[current_ko] = '; '.join(pathways[:5])
            except Exception:
                failed += 1
            time.sleep(0.2)
            if i % 500 == 0 and i > 0:
                print(f"    Progress: {i}/{len(ko_list)} KOs processed...")
        print(f"    Mapped {len(ko_to_pathway)} KOs ({failed} batches failed)")
        return ko_to_pathway

    print(f"    Fetching names for {len(all_kos)} unique KO numbers...")
    ko_to_pathway = {}
    failed = 0

    for ko in sorted(all_kos):
        url = f"https://rest.kegg.jp/get/ko:{ko}"
        try:
            with urllib.request.urlopen(url, timeout=10) as resp:
                text = resp.read().decode('utf-8')
            # Parse pathway lines from KEGG entry
            pathways = []
            in_pathway = False
            for line in text.splitlines():
                if line.startswith('PATHWAY'):
                    in_pathway = True
                    parts = line.split(None, 2)
                    if len(parts) >= 3:
                        pathways.append(parts[2].strip())
                elif in_pathway and line.startswith(' '):
                    parts = line.split(None, 2)
                    if len(parts) >= 2:
                        pathways.append(parts[-1].strip())
                elif in_pathway and not line.startswith(' '):
                    in_pathway = False
            if pathways:
                ko_to_pathway[ko] = '; '.join(pathways[:5])  # cap at 5 pathways
        except Exception:
            failed += 1
        time.sleep(0.1)  # KEGG rate limit

    print(f"    Mapped {len(ko_to_pathway)} KOs to pathway names ({failed} failed)")
    return ko_to_pathway


def map_kegg_pathways(ko_col, ko_to_pathway):
    """Map a column of KO numbers to pathway names."""
    def _map(val):
        if not val or pd.isna(val):
            return ''
        pathways = []
        seen = set()
        for ko in str(val).split(','):
            ko = ko.strip().replace('ko:', '').strip()
            if ko.startswith('map') or not ko.startswith('K'):
                continue
            if ko in ko_to_pathway:
                # Split on semicolon and deduplicate individual pathway names
                for pw in ko_to_pathway[ko].split('; '):
                    pw = pw.strip()
                    # Remove generic trailing words
                    pw = pw.rstrip('; pathways; metabolism; signaling')
                    pw = pw.strip().rstrip(';').strip()
                    if pw and pw not in seen and len(pw) > 3:
                        seen.add(pw)
                        pathways.append(pw)
        return '; '.join(pathways[:10])  # cap at 10 pathways per gene
    return ko_col.apply(_map)


# ============================================================
# PARSING FUNCTIONS
# ============================================================

def parse_gff3(gff3_file, gene_to_protein_csv=''):
    """Parse GFF3 to extract genomic coordinates and exon counts per transcript.

    For standard GFF3 (Gfas, Pdam): mRNA ID matches protein FASTA ID directly.
    For NCBI GFF3 (Acer, Ofav, Amur): mRNA ID uses gnl|WGS|... format; protein
    FASTA IDs are KAK accessions. A gene_to_protein CSV is used to bridge them
    via the shared locus_tag field.

    The CSV must have columns: gene_id (protein FASTA accession), protein_id (locus tag).
    """
    print("  Parsing GFF3 for genomic coordinates...")
    if not gff3_file or not Path(gff3_file).exists():
        print(f"    WARNING: GFF3 not found: {gff3_file}")
        return pd.DataFrame()

    # Build locus_tag -> protein FASTA ID mapping if CSV provided
    locus_to_fasta = {}
    if gene_to_protein_csv and Path(gene_to_protein_csv).exists():
        import csv as _csv
        with open(gene_to_protein_csv) as _f:
            for row in _csv.DictReader(_f):
                # gene_id = KAK254xxxx (protein FASTA), protein_id = P5673_xxxxxx (locus tag)
                locus_to_fasta[row['protein_id']] = row['gene_id']
        print(f"    Gene-to-protein map: {len(locus_to_fasta):,} locus tag -> FASTA ID entries")

    rows = []
    exon_counts = {}
    # Track mRNA internal ID -> protein FASTA ID for exon parent matching
    mrna_id_to_fasta = {}

    with open(gff3_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            attr_dict = {}
            for attr in attrs.split(';'):
                attr = attr.strip()
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attr_dict[k] = v

            if ftype in ('mRNA', 'transcript'):
                mrna_id = attr_dict.get('ID', '')
                locus_tag = attr_dict.get('locus_tag', '')

                # Determine protein FASTA ID
                if locus_to_fasta and locus_tag:
                    fasta_id = locus_to_fasta.get(locus_tag, '')
                    gff_protein_id = locus_tag  # store locus tag as bridge ID
                else:
                    fasta_id = mrna_id
                    gff_protein_id = ''

                if not fasta_id:
                    continue

                mrna_id_to_fasta[mrna_id] = fasta_id
                rows.append({
                    'gene_id':          fasta_id,
                    'acer_locus_tag':   locus_tag if locus_to_fasta else '',
                    'scaffold':         seqid,
                    'gene_start':       int(start),
                    'gene_end':         int(end),
                    'strand':           strand,
                    'gene_length_bp':   int(end) - int(start) + 1,
                })
                exon_counts[fasta_id] = 0

            elif ftype == 'exon':
                parent_mrna = attr_dict.get('Parent', '')
                fasta_id = mrna_id_to_fasta.get(parent_mrna, '')
                if fasta_id and fasta_id in exon_counts:
                    exon_counts[fasta_id] += 1

    df = pd.DataFrame(rows)
    if not df.empty:
        df['n_exons'] = df['gene_id'].map(exon_counts).fillna(0).astype(int)
        mapped = (df['acer_locus_tag'] != '').sum() if 'acer_locus_tag' in df.columns else 0
        if mapped:
            print(f"    GFF3: {len(df):,} transcripts, {df['scaffold'].nunique():,} scaffolds ({mapped:,} mapped via locus tag)")
        else:
            print(f"    GFF3: {len(df):,} transcripts, {df['scaffold'].nunique():,} scaffolds")
    return df



def parse_fasta_lengths(fasta_file):
    """Return dicts of gene_id → protein length gene_id → sequence."""
    if not fasta_file or not Path(fasta_file).exists():
        return {}, {}
    lengths = {}
    sequences = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id:
                    seq = ''.join(current_seq)
                    lengths[current_id]   = len(seq)
                    sequences[current_id] = seq
                current_id  = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seq = ''.join(current_seq)
        lengths[current_id]   = len(seq)
        sequences[current_id] = seq
    return lengths, sequences


def parse_eggnog(eggnog_file):
    """Parse eggNOG-mapper v2 annotations.

    Column order (0-indexed):
    0  #query
    1  seed_ortholog
    2  evalue
    3  score
    4  eggNOG_OGs
    5  max_annot_lvl
    6  COG_category       ← single-letter codes e.g. K, T, S
    7  Description
    8  Preferred_name     ← gene name
    9  GOs                ← GO terms pipe-separated
    10 EC
    11 KEGG_ko
    12 KEGG_Pathway
    13 KEGG_Module
    14 KEGG_Reaction
    15 KEGG_rclass
    16 BRITE
    17 KEGG_TC
    18 CAZy
    19 BiGG_Reaction
    20 PFAMs
    """
    print("  Parsing eggNOG annotations...")
    if not Path(eggnog_file).exists():
        print(f"    WARNING: {eggnog_file} not found. Returning empty.")
        return pd.DataFrame()

    rows = []
    with open(eggnog_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue
            rows.append({
                'gene_id':       parts[0],
                'eggnog_ortholog':    parts[1]  if parts[1]  != '-' else '',
                'eggnog_human_gene':  parts[8]  if len(parts) > 8  and parts[8]  != '-' else '',
                'eggnog_description': parts[7]  if len(parts) > 7  and parts[7]  != '-' else '',
                'eggnog_go_terms':    parts[9]  if len(parts) > 9  and parts[9]  != '-' else '',
                'eggnog_kegg_ko':     parts[11] if len(parts) > 11 and parts[11] != '-' else '',
                'eggnog_kegg_pathway':(','.join(sorted(set(
                    k.replace('ko:','').replace('map','').strip()
                    for k in parts[12].split(',')
                    if k.strip() and k.strip() != '-'
                    and not k.strip().startswith('map')
                ))) if len(parts) > 12 and parts[12] != '-' else ''),
                'eggnog_cog_cat':     parts[6]  if len(parts) > 6  and parts[6]  != '-' else '',
                'eggnog_pfams':       parts[20] if len(parts) > 20 and parts[20] != '-' else '',
            })

    df = pd.DataFrame(rows)
    print(f"    eggNOG: {len(df):,} annotations")
    return df


def parse_orthofinder(of_results_dir, species):
    """Parse OrthoFinder ortholog tables."""
    print("  Parsing OrthoFinder results...")
    of_dir = Path(of_results_dir)
    # Auto-detect dated subdirectory (e.g. Results_Jun05, Results_Jun16)
    if not (of_dir / "Orthologues").exists():
        dated = sorted([d for d in of_dir.iterdir() if d.is_dir() and d.name.startswith("Results_")])
        if dated:
            of_dir = dated[-1]
            print(f"    OrthoFinder: auto-detected results dir: {of_dir.name}")
    orthologue_dir = of_dir / "Orthologues" / f"Orthologues_{species}"

    if not orthologue_dir.exists():
        print(f"    WARNING: {orthologue_dir} not found. Returning empty.")
        return pd.DataFrame()

    result = {}

    for tsv_file in orthologue_dir.glob("*.tsv"):
        other_species = tsv_file.stem.replace(f"{species}__v__", "").replace(f"__v__{species}", "")
        try:
            df = pd.read_csv(tsv_file, sep='\t', header=0)
        except Exception:
            continue

        cols = df.columns.tolist()
        query_col = next((c for c in cols if species in c), None)
        other_col = next((c for c in cols if other_species in c and c != query_col), None)
        og_col    = next((c for c in cols if 'orthogroup' in c.lower() or c == cols[0]), cols[0])

        if not query_col or not other_col:
            continue

        for _, row in df.iterrows():
            orthogroup  = str(row[og_col])
            query_genes = str(row[query_col]).split(',') if pd.notna(row[query_col]) else []
            other_genes = str(row[other_col]).split(',') if pd.notna(row[other_col]) else []
            query_genes = [g.strip() for g in query_genes if g.strip()]
            other_genes = [g.strip() for g in other_genes if g.strip()]

            relationship = "one-to-one"   if len(query_genes) == 1 and len(other_genes) == 1 else \
                           "one-to-many"  if len(query_genes) == 1 else "many-to-many"

            for qg in query_genes:
                if qg not in result:
                    result[qg] = {'gene_id': qg, 'of_orthogroup': orthogroup,
                                  'of_human_ortholog': '', 'of_mouse_ortholog': '',
                                  'of_nvec_ortholog': '', 'of_relationship': relationship}
                if 'Hsap' in other_species:
                    result[qg]['of_human_ortholog'] = ','.join(other_genes[:3])
                elif 'Mmus' in other_species:
                    result[qg]['of_mouse_ortholog'] = ','.join(other_genes[:3])
                elif 'Nvec' in other_species:
                    result[qg]['of_nvec_ortholog']  = ','.join(other_genes[:3])

    df_out = pd.DataFrame(list(result.values())) if result else pd.DataFrame()
    print(f"    OrthoFinder: {len(df_out):,} genes with orthologs")
    return df_out


def parse_interproscan(ips_file):
    """Parse InterProScan TSV output.
    NOTE: IPS 5.78 does not include SignalP or TMHMM — those are run separately.
    """
    print("  Parsing InterProScan results...")
    if not Path(ips_file).exists():
        print(f"    WARNING: {ips_file} not found. Returning empty.")
        return pd.DataFrame()

    pfam_hits    = defaultdict(set)
    panther_hits = defaultdict(set)
    go_hits      = defaultdict(set)

    with open(ips_file) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue

            gene_id = parts[0]
            db      = parts[3]
            db_accn = parts[4] if len(parts) > 4 else ''
            go_str  = parts[13] if len(parts) > 13 else ''

            if db == 'Pfam':
                pfam_hits[gene_id].add(db_accn)
            elif db == 'PANTHER':
                panther_hits[gene_id].add(db_accn)

            if go_str and go_str != '-':
                for go in go_str.split('|'):
                    # Strip source annotations e.g. GO:0001234(PANTHER) → GO:0001234
                    go_clean = re.sub(r'\(.*?\)', '', go).strip()
                    if go_clean.startswith('GO:'):
                        go_hits[gene_id].add(go_clean)

    all_genes = set(pfam_hits) | set(panther_hits)

    rows = []
    for gene_id in all_genes:
        rows.append({
            'gene_id':      gene_id,
            'interpro_pfam':     ','.join(sorted(pfam_hits[gene_id])),
            'interpro_panther':  ','.join(sorted(panther_hits[gene_id])),
            'interpro_go_terms': ','.join(sorted(go_hits[gene_id])),
        })

    df = pd.DataFrame(rows)
    print(f"    InterProScan: {len(df):,} annotated genes")
    return df


def parse_rbh(rbh_file):
    """Parse RBH SwissProt hits."""
    print("  Parsing RBH hits...")
    if not Path(rbh_file).exists():
        print(f"    WARNING: {rbh_file} not found. Returning empty.")
        return pd.DataFrame()
    df = pd.read_csv(rbh_file, sep='\t')
    print(f"    RBH: {len(df):,} reciprocal best hits")
    return df


def parse_signalp(signalp_file):
    """Parse SignalP 6.0 prediction_results.txt.
    Format: ID  Prediction  OTHER_prob  SP_prob  [CS_Position]
    Column 3 (index 3) is the SP probability.
    """
    print("  Parsing SignalP predictions...")
    if not Path(signalp_file).exists():
        print(f"    WARNING: {signalp_file} not found. Returning empty.")
        return pd.DataFrame()

    rows = []
    with open(signalp_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split("	")
            if len(parts) >= 2:
                rows.append({
                    'gene_id':       parts[0].split()[0],
                    'signalp_prediction': parts[1],
                    'signalp_prob':       float(parts[3]) if len(parts) > 3 else None,
                })

    df = pd.DataFrame(rows)
    print(f"    SignalP: {len(df):,} predictions")
    return df


def parse_tmhmm(tmhmm_gff3):
    """Parse DeepTMHMM TMRs.gff3 output."""
    print("  Parsing DeepTMHMM predictions...")
    if not tmhmm_gff3 or not Path(tmhmm_gff3).exists():
        print(f"    WARNING: {tmhmm_gff3} not found. Returning empty.")
        return pd.DataFrame()

    rows = []
    current_gene = None
    n_tmrs = 0
    topology_types = set()

    with open(tmhmm_gff3) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("##") or not line:
                continue
            if line.startswith("# ") and "Number of predicted TMRs" in line:
                parts = line.split()
                current_gene = parts[1]
                n_tmrs = int(parts[-1])
            elif line.startswith("//"):
                if current_gene is not None:
                    if n_tmrs == 0:
                        topo = "cytoplasmic" if "inside" in topology_types else "extracellular"
                    elif "signal" in topology_types:
                        topo = f"membrane_signal_{n_tmrs}TM"
                    else:
                        topo = f"membrane_{n_tmrs}TM"
                    rows.append({
                        'gene_id':   current_gene,
                        'tmhmm_topology': topo,
                        'tmhmm_n_tmrs':   int(n_tmrs),
                    })
                current_gene = None
                n_tmrs = 0
                topology_types = set()
            elif not line.startswith("#"):
                parts = line.split('\t')
                if len(parts) >= 2:
                    topology_types.add(parts[1].lower())

    df = pd.DataFrame(rows)
    if not df.empty:
        membrane = (df['tmhmm_n_tmrs'] > 0).sum()
        print(f"    DeepTMHMM: {len(df):,} genes, {membrane:,} with >=1 TM helix")
    return df


def parse_tier3(tier3_file):
    """Parse HMMER domain table output (--domtblout).
    Stores HMM profile NAME (not Pfam accession) in tier3_domain_hits.
    """
    print("  Parsing Tier 3 HMM hits...")
    if not Path(tier3_file).exists():
        print(f"    WARNING: {tier3_file} not found. Returning empty.")
        return pd.DataFrame()

    hits = defaultdict(lambda: {'domains': set(), 'best_evalue': 999.0})

    with open(tier3_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 20:
                continue
            target   = parts[0]    # gene/protein
            hmm_name = parts[3]    # HMM profile NAME field
            dom_eval = float(parts[12])

            if dom_eval < 1e-3:
                hits[target]['domains'].add(hmm_name)
                hits[target]['best_evalue'] = min(hits[target]['best_evalue'], dom_eval)

    rows = []
    for gene_id, data in hits.items():
        rows.append({
            'gene_id':      gene_id,
            'tier3_domain_hits': ','.join(sorted(data['domains'])),
            'tier3_best_evalue': data['best_evalue'],
        })

    df = pd.DataFrame(rows)
    print(f"    Tier 3 HMM: {len(df):,} genes with domain hits")
    return df


def parse_tier4(tier4_file):
    """Parse Tier 4 BLAST hits."""
    print("  Parsing Tier 4 BLAST hits...")
    if not Path(tier4_file).exists():
        print(f"    WARNING: {tier4_file} not found. Returning empty.")
        return pd.DataFrame()
    df = pd.read_csv(tier4_file, sep='\t')
    df = df.sort_values('tier4_evalue').drop_duplicates(subset=list(df.columns[:1]), keep='first')
    print(f"    Tier 4 BLAST: {len(df):,} genes with hits")
    return df


# ============================================================
# CLASSIFICATION LOGIC
# ============================================================

def nonempty(val):
    """True only if val is a non-empty, non-NaN string."""
    if val is None:
        return False
    if isinstance(val, float):
        return False
    s = str(val).strip()
    return bool(s) and s not in ('', 'nan', 'NaN', 'None', '-')


def assign_classification(row):
    """
    Assign confidence tier and protein group for one gene.

    Priority:
      1. Tier 4 coral-specific BLAST (highest specificity)
      2. OrthoFinder ortholog + RBH agreement → Gold
      3. OrthoFinder ortholog alone → Silver
      4. RBH alone → Silver
      5. Tier 3 HMM domain (accession-based, then name-based) → Bronze
      6. eggNOG gene name keyword → Bronze
      7. InterProScan Pfam alone → Bronze
      8. Unclassified
    """
    confidence    = "Unclassified"
    protein_group = "Unclassified"
    subgroup      = ""

    has_of_human = nonempty(row.get('of_human_ortholog', ''))
    has_of_mouse = nonempty(row.get('of_mouse_ortholog', ''))
    has_rbh      = nonempty(row.get('rbh_gene_name', ''))
    has_eggnog   = nonempty(row.get('eggnog_human_gene', ''))
    has_tier3    = nonempty(row.get('tier3_domain_hits', ''))
    has_tier4    = nonempty(row.get('tier4_source', ''))
    has_pfam     = nonempty(row.get('interpro_pfam', ''))

    # ── Tier 4: coral-specific families ──────────────────────
    if has_tier4:
        t4_source = str(row.get('tier4_source', ''))
        if t4_source in TIER4_TO_GROUP:
            protein_group, subgroup = TIER4_TO_GROUP[t4_source]
        else:
            protein_group = f"Tier4_{t4_source}"
            subgroup = t4_source
        confidence = "Gold" if (has_rbh or has_of_human) else "Informative"
        return confidence, protein_group, subgroup

    # ── Gold: OrthoFinder + RBH ───────────────────────────────
    if (has_of_human or has_of_mouse) and has_rbh:
        confidence = "Gold"
        gene_name = row.get('of_human_ortholog') or row.get('rbh_gene_name', '')
        pg, sg = match_gene_keyword(gene_name)
        if pg:
            protein_group, subgroup = pg, sg

    # ── Silver: OrthoFinder only ──────────────────────────────
    elif has_of_human or has_of_mouse:
        confidence = "Silver"
        gene_name = row.get('of_human_ortholog') or row.get('of_mouse_ortholog', '')
        pg, sg = match_gene_keyword(gene_name)
        if pg:
            protein_group, subgroup = pg, sg

    # ── Silver: RBH only ──────────────────────────────────────
    elif has_rbh:
        confidence = "Silver"
        pg, sg = match_gene_keyword(row.get('rbh_gene_name', ''))
        if pg:
            protein_group, subgroup = pg, sg

    # ── Tier 3: HMM domain annotation ────────────────────────
    if has_tier3 and protein_group == "Unclassified":
        domain_hits = str(row.get('tier3_domain_hits', ''))

        # 1. Pfam accession-based lookup
        matched = False
        pfam_ids = str(row.get('interpro_pfam', '')).split(',')
        for pfam_id in pfam_ids:
            pfam_id = pfam_id.strip()
            if pfam_id in TIER3_TO_GROUP:
                protein_group, subgroup = TIER3_TO_GROUP[pfam_id]
                if confidence not in ("Gold", "Silver"):
                    confidence = "Bronze"
                matched = True
                break

        # 2. HMM name-based lookup (longest key first)
        if not matched:
            for hmm_name, (pg, sg) in sorted(
                    TIER3_NAME_TO_GROUP.items(), key=lambda x: len(x[0]), reverse=True):
                if hmm_name in domain_hits:
                    protein_group, subgroup = pg, sg
                    if confidence not in ("Gold", "Silver"):
                        confidence = "Bronze"
                    break

    # ── Bronze: eggNOG gene name keyword ─────────────────────
    if has_eggnog and protein_group == "Unclassified":
        pg, sg = match_gene_keyword(row.get('eggnog_human_gene', ''))
        if pg:
            protein_group, subgroup = pg, sg
            if confidence == "Unclassified":
                confidence = "Bronze"

    # ── Bronze: Pfam accession from InterProScan alone ───────
    if has_pfam and protein_group == "Unclassified":
        for pfam_id in str(row.get('interpro_pfam', '')).split(','):
            pfam_id = pfam_id.strip()
            if pfam_id in TIER3_TO_GROUP:
                protein_group, subgroup = TIER3_TO_GROUP[pfam_id]
                confidence = "Bronze"
                break

    if protein_group != "Unclassified" and confidence == "Unclassified":
        confidence = "Bronze"

    return confidence, protein_group, subgroup


def assign_annotation_status(row):
    """
    Light/Partial/Dark classification per Stephens et al. 2026 GBE (evag072).
    Light:   named ortholog OR informative functional description
    Partial: named Pfam domain or Tier 3/4 hit but no ortholog/description
    Dark:    no ortholog, no description, no named domain
    """
    AMBIGUOUS = re.compile(
        r"uncharacterized|hypothetical|unknown function|DUF\d|domain of unknown|"
        r"predicted protein|unnamed protein|function unknown|putative uncharacterized|"
        r"no function|conserved protein|novel protein", re.IGNORECASE)
    DOMAIN_ONLY = re.compile(
        r"^(helix.loop.helix|helicase|kinase|zinc finger|WD repeat|ankyrin|"
        r"TPR repeat|armadillo|HEAT repeat|leucine.rich|EGF.like|"
        r"immunoglobulin|fibronectin|cadherin|collagen|repeat|domain)"
        r"[\s\-]*(domain|repeat|family|motif|like)?$", re.IGNORECASE)

    best_name  = str(row.get("best_human_gene_name", "")).strip()
    desc       = str(row.get("eggnog_description", "")).strip()
    pfam_names = str(row.get("interpro_pfam_names", "")).strip()
    pfam_acc   = str(row.get("interpro_pfam", "")).strip()
    t3         = str(row.get("tier3_domain_hits", "")).strip()
    t4         = str(row.get("tier4_source", "")).strip()

    has_ortholog = (best_name not in ("", "nan", "None") and
                    not AMBIGUOUS.search(best_name))
    has_informative_desc = (desc not in ("", "nan", "None", "-") and
                            not AMBIGUOUS.search(desc) and
                            not DOMAIN_ONLY.match(desc) and len(desc) > 10)
    has_named_pfam = (pfam_names not in ("", "nan", "None") and
                      not re.match(r"^Domain of unknown function|^DUF|^Uncharacterised",
                                   pfam_names, re.IGNORECASE))
    has_any_domain = any(x not in ("", "nan", "None") for x in [pfam_acc, t3, t4])

    if has_ortholog or has_informative_desc:
        return "Light"
    elif has_named_pfam or has_any_domain:
        return "Partial"
    else:
        return "Dark"



# Secondary group rules: (primary_group) -> {keyword_in_annotation: secondary_group}
# Only genuinely dual-role genes — keywords must indicate structural/functional duality
# Excludes: TRAF/RING zinc fingers (not TFs), guanylate kinases (not signal kinases)
SECONDARY_GROUP_RULES = {
    "Innate_Immunity_PRR": {
        "c-type lectin":        "Lectin_Receptor",
        "lectin c-type":        "Lectin_Receptor",
        "galectin":             "Lectin_Receptor",
        "fibronectin":          "ECM_Adhesion",
        "collagen":             "ECM_Adhesion",
        "protein kinase":       "Kinase",
        "traf":                 "Innate_Immunity_Signalling",
    },
    "Lectin_Receptor": {
        "nacht":                "Innate_Immunity_PRR",
        "toll":                 "Innate_Immunity_PRR",
        "complement":         "Innate_Immunity_PRR",
        "fibronectin":          "ECM_Adhesion",
        "collagen":             "ECM_Adhesion",
        "protein kinase":       "Kinase",
    },
    "Kinase": {
        "tir domain":           "Innate_Immunity_PRR",
        "toll-interleukin":     "Innate_Immunity_PRR",
        "nacht":                "Innate_Immunity_PRR",
        "homeobox":             "Transcription_Factor",
        "fork head":            "Transcription_Factor",
        "dna-binding":          "Transcription_Factor",
        "traf":                 "Innate_Immunity_Signalling",
        "caspase":              "Innate_Immunity_Effector",
    },
    "ECM_Adhesion": {
        "c-type lectin":        "Lectin_Receptor",
        "lectin c-type":        "Lectin_Receptor",
        "galectin":             "Lectin_Receptor",
        "nacht":                "Innate_Immunity_PRR",
        "protein kinase":       "Kinase",
        "ion channel":          "Ion_Channel",
    },
    "Innate_Immunity_Signalling": {
        "protein kinase":       "Kinase",
        "homeobox":             "Transcription_Factor",
        "fork head":            "Transcription_Factor",
        "helix-turn-helix":     "Transcription_Factor",
        "dna-binding":          "Transcription_Factor",
    },
    "GPCR": {
        "neuropeptide":         "Neuron_Neuropeptide",
        "ion channel":          "Ion_Channel",
    },
    "Ion_Channel": {
        "neuropeptide":         "Neuron_Neuropeptide",
        "protein kinase":       "Kinase",
    },
    "Innate_Immunity_Effector": {
        "c-type lectin":        "Lectin_Receptor",
        "lectin c-type":        "Lectin_Receptor",
        "complement":           "Innate_Immunity_PRR",
        "protein kinase":       "Kinase",
    },
    "Calcification_Universal": {
        "protein kinase":       "Kinase",
        "ion channel":          "Ion_Channel",
        "fibronectin":          "ECM_Adhesion",
    },
    "Transcription_Factor": {
        "protein kinase":       "Kinase",
        "traf":                 "Innate_Immunity_Signalling",
    },
    "Neuron_Neuropeptide": {
        "7 transmembrane":      "GPCR",
        "rhodopsin":            "GPCR",
        "ion channel":          "Ion_Channel",
    },
}

def assign_secondary_group(row):
    """Assign a secondary protein group where a gene genuinely has dual roles."""
    primary = str(row.get("protein_group", "")).strip()
    if primary in ("Unclassified", "", "nan"):
        return ""
    rules = SECONDARY_GROUP_RULES.get(primary, {})
    if not rules:
        return ""
    text = " ".join([
        str(row.get("interpro_pfam_names", "") or ""),
        str(row.get("eggnog_description", "") or ""),
        str(row.get("tier3_domain_hits", "") or ""),
        str(row.get("best_human_gene_name", "") or ""),
    ]).lower()
    for keyword, secondary in rules.items():
        if keyword in text and secondary != primary:
            return secondary
    return ""



def assign_cellchat_role(row):
    """Assign CellChat role based on SignalP 6 + DeepTMHMM + protein group."""
    sp    = str(row.get('signalp_prediction', '')).upper()
    tm    = str(row.get('tmhmm_topology', '')).lower()
    group = str(row.get('protein_group', ''))

    n_tm = row.get('tmhmm_n_tmrs', 0)
    try:
        n_tm = int(n_tm) if n_tm == n_tm else 0
    except (ValueError, TypeError):
        n_tm = 0

    RECEPTOR_GROUPS = {
        'GPCR', 'Ion_Channel', 'Innate_Immunity_PRR',
        'Lectin_Receptor', 'Kinase'
    }

    is_membrane = n_tm >= 1 or 'membrane' in tm or 'tm' in tm

    if sp in ('SP', 'LIPO', 'TAT', 'TATLIPO'):
        if group in RECEPTOR_GROUPS or is_membrane:
            return "Receptor"
        return "Ligand"

    if group in RECEPTOR_GROUPS:
        return "Receptor"

    if is_membrane:
        return "Receptor_candidate"

    return "Neither"


def make_best_human_gene_name(row):
    """Single clean human gene symbol for heatmaps/volcano plots."""
    # Priority: OrthoFinder human > eggNOG name > RBH gene name
    for col in ['of_human_ortholog', 'eggnog_human_gene', 'rbh_gene_name']:
        val = str(row.get(col, '')).strip()
        if val and val not in ('', 'nan', 'None', '-'):
            # Take first entry if multiple, strip species suffix
            name = val.split(',')[0].strip()
            name = re.sub(r'_HUMAN$|_MOUSE$|_RAT$|_BOVIN$|_XENLA$|_DANRE$', '', name)
            if name:
                return name
    return ''


def make_all_go_terms(row):
    """Union of eggNOG GO terms + InterProScan GO terms, deduplicated."""
    go_set = set()
    for col in ['eggnog_go_terms', 'interpro_go_terms']:
        val = str(row.get(col, '')).strip()
        if val and val not in ('', 'nan', 'None', '-'):
            for term in re.split(r'[,|]', val):
                term = re.sub(r'\(.*?\)', '', term).strip()
                if term.startswith('GO:') and len(term) >= 10:
                    go_set.add(term)
    return '|'.join(sorted(go_set)) if go_set else ''


def make_cog_description(cog_cat):
    """Expand single-letter COG category codes to descriptions."""
    if not cog_cat or pd.isna(cog_cat) or str(cog_cat) in ('', 'nan', '-'):
        return ''
    codes = str(cog_cat).strip()
    descs = []
    for c in codes:
        if c in COG_DESCRIPTIONS:
            descs.append(COG_DESCRIPTIONS[c])
    return '; '.join(descs) if descs else cog_cat


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--species',        required=True,  help='Species code e.g. Gfas')
    parser.add_argument('--run_dir',        required=True,  help='Run directory')
    parser.add_argument('--proteome_fasta', default='',     help='Proteome FASTA (seeds gene universe + lengths)')
    parser.add_argument('--gff3', default='', help='GFF3 annotation file for genomic coordinates')
    parser.add_argument('--eggnog',         required=True)
    parser.add_argument('--orthofinder',    required=True)
    parser.add_argument('--interproscan',   required=True)
    parser.add_argument('--rbh',            required=True)
    parser.add_argument('--signalp',        required=True)
    parser.add_argument('--tmhmm',          default='',     help='DeepTMHMM TMRs.gff3')
    parser.add_argument('--tier3',          required=True)
    parser.add_argument('--tier4',          required=True)
    parser.add_argument('--out_table',      required=True)
    parser.add_argument('--out_genelists',  required=True)
    parser.add_argument('--skip_kegg_api',  action='store_true',
                        help='Skip KEGG REST API calls (faster, no pathway names)')
    parser.add_argument('--genome_version', default='',
                        help='Genome version string for annotation summary (e.g. Gfas_v1.0)')
    parser.add_argument('--proteome_source', default='',
                        help='Proteome source URL for annotation summary')
    parser.add_argument('--staging_dir', default='',
                        help='Root staging directory for GitHub upload '
                             '(e.g. /nethome/kxw755/github_upload). '
                             'If provided, outputs are copied here after completion.')
    parser.add_argument('--gene_to_protein', default='',
                        help='CSV mapping protein FASTA IDs to locus tags (for NCBI genomes). '
                             'Required columns: gene_id (protein FASTA accession), '
                             'protein_id (locus tag e.g. P5673_xxxxxx). '
                             'Enables GFF3 coordinate parsing for NCBI-formatted genomes.')
    args = parser.parse_args()

    species = args.species
    global gene_id_col
    gene_id_col = f"{species.lower()}_gene_id"

    # ── Set up log tee: write to terminal and to a log file ────
    import sys as _sys
    out_final_dir = Path(args.out_table).parent
    out_final_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_final_dir / f"{species}_merge_run.out"

    class _Tee:
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush()
        def flush(self):
            for f in self.files:
                f.flush()
        def fileno(self):
            return self.files[0].fileno()

    _log_fh = open(log_path, 'w')
    _sys.stdout = _Tee(_sys.__stdout__, _log_fh)

    print(f"\n{'='*60}")
    print(f"  08_merge_annotate.py v3 — Species: {species}")

    # Load Pfam accession → name lookup
    PFAM_NAMES = {}
    pfam_names_file = Path("/scratch/dark_genes/annotation_pipeline/databases/pfam_names.tsv")
    if pfam_names_file.exists():
        with open(pfam_names_file) as _f:
            for _line in _f:
                _parts = _line.strip().split('\t')
                if len(_parts) == 2:
                    PFAM_NAMES[_parts[0]] = _parts[1]
        print(f"  Loaded {len(PFAM_NAMES):,} Pfam name mappings")
    else:
        print(f"  WARNING: pfam_names.tsv not found — interpro_pfam_names will be empty")
    print(f"{'='*60}\n")

    # ── Step 1: Parse annotation layers ──────────────────────
    print("--- Parsing annotation layers ---")
    df_eggnog = parse_eggnog(args.eggnog)
    df_of     = parse_orthofinder(args.orthofinder, species)
    df_ips    = parse_interproscan(args.interproscan)
    df_rbh    = parse_rbh(args.rbh)
    df_sp     = parse_signalp(args.signalp)
    df_gff    = parse_gff3(args.gff3, args.gene_to_protein) if args.gff3 else pd.DataFrame()
    df_tm     = parse_tmhmm(args.tmhmm) if args.tmhmm else pd.DataFrame()
    df_t3     = parse_tier3(args.tier3)
    df_t4     = parse_tier4(args.tier4)

    # ── Step 2: Build gene universe ───────────────────────────
    print("\n--- Building gene universe ---")
    all_genes = set()
    if args.proteome_fasta and Path(args.proteome_fasta).exists():
        protein_lengths, protein_sequences = parse_fasta_lengths(args.proteome_fasta)
        all_genes = set(protein_lengths.keys())
        print(f"  Seeded from proteome FASTA: {len(all_genes):,} genes")
    else:
        protein_lengths = {}
        for df in [df_eggnog, df_of, df_ips, df_rbh, df_sp, df_t3, df_t4]:
            if not df.empty and gene_id_col in df.columns:
                all_genes.update(df[gene_id_col].dropna().astype(str))
        print(f"  Seeded from annotation outputs: {len(all_genes):,} genes")

    master = pd.DataFrame({gene_id_col: sorted(all_genes)})
    print(f"  Total genes: {len(master):,}")

    # ── Step 3: Merge layers ──────────────────────────────────
    # Rename any legacy gene_id column to species-specific gene_id_col
    for _df in [df_eggnog, df_of, df_ips, df_rbh, df_sp, df_t3, df_t4, df_gff, df_tm]:
        if _df.empty:
            continue
        for _old in list(_df.columns):
            if _old.endswith("_gene_id") or _old == "gene_id":
                if _old != gene_id_col:
                    _df.rename(columns={_old: gene_id_col}, inplace=True)
                break
    print("\n--- Merging annotation layers ---")

    def safe_merge(base, right, on=None, how="left"):
        if on is None: on = gene_id_col
        if right is None or (hasattr(right, 'empty') and right.empty):
            return base
        overlap = [c for c in right.columns if c in base.columns and c != on]
        if overlap:
            right = right.drop(columns=overlap)
        return base.merge(right, on=on, how=how)

    master = safe_merge(master, df_eggnog)
    master = safe_merge(master, df_of)
    master = safe_merge(master, df_ips)
    master = safe_merge(master, df_rbh)
    master = safe_merge(master, df_sp)    # real SignalP — no IPS collision now
    master = safe_merge(master, df_t3)
    master = safe_merge(master, df_t4)
    if not df_gff.empty:
        # Drop acer_locus_tag column if empty (non-NCBI species where no mapping was used)
        if 'acer_locus_tag' in df_gff.columns and (df_gff['acer_locus_tag'].fillna('') == '').all():
            df_gff = df_gff.drop(columns=['acer_locus_tag'])
        master = safe_merge(master, df_gff)
    if not df_tm.empty:
        master = safe_merge(master, df_tm)
        if 'tmhmm_n_tmrs' in master.columns:
            master['tmhmm_n_tmrs'] = master['tmhmm_n_tmrs'].fillna(0).astype(int)

    # Add protein length
    if protein_lengths:
        master['protein_length_aa'] = master[gene_id_col].map(protein_lengths).fillna(0).astype(int)
    if protein_sequences:
        master['protein_sequence']  = master[gene_id_col].map(protein_sequences).fillna('')

    # Fill NaN in string columns
    str_cols = [c for c in master.columns if master[c].dtype == object]
    master[str_cols] = master[str_cols].fillna('')
    for col in ['tier4_source', 'tier4_blast_hit', 'rbh_gene_name',
                'of_human_ortholog', 'of_mouse_ortholog', 'of_nvec_ortholog',
                'eggnog_human_gene', 'tier3_domain_hits',
                'signalp_prediction', 'tmhmm_topology']:
        if col in master.columns:
            master[col] = master[col].fillna('').astype(str).replace('nan', '')

    print(f"  Merged table: {master.shape[0]:,} genes × {master.shape[1]} columns")

    # ── Step 4: Classify genes ────────────────────────────────
    print("\n--- Classifying genes ---")
    classifications = []
    for _, row in tqdm(master.iterrows(), total=len(master), desc="  Classifying"):
        conf, pg, sg = assign_classification(row)
        classifications.append((conf, pg, sg))

    master['confidence_tier']  = [c[0] for c in classifications]
    master['protein_group']    = [c[1] for c in classifications]
    master['protein_subgroup'] = [c[2] for c in classifications]
    master['protein_group_secondary'] = master.apply(assign_secondary_group, axis=1).replace('', pd.NA)
    n_sec = master['protein_group_secondary'].notna().sum()
    print(f"  protein_group_secondary: {n_sec:,} genes with dual roles")
    master['cellchat_role']    = master.apply(assign_cellchat_role, axis=1)
    # ── Step 5: Post-processing columns ──────────────────────
    print("\n--- Adding derived columns ---")

    # best_human_gene_name — single clean gene symbol
    master['best_human_gene_name'] = master.apply(make_best_human_gene_name, axis=1)
    print(f"  best_human_gene_name: {(master['best_human_gene_name'] != '').sum():,} genes")

    # all_go_terms — union of eggNOG + IPS GO terms, clean format
    master['all_go_terms'] = master.apply(make_all_go_terms, axis=1)
    print(f"  all_go_terms: {(master['all_go_terms'] != '').sum():,} genes")

    # cog_description — expand single-letter COG codes
    if 'eggnog_cog_cat' in master.columns:
        master['cog_description'] = master['eggnog_cog_cat'].apply(make_cog_description)
        print(f"  cog_description: {(master['cog_description'] != '').sum():,} genes")


    # interpro_pfam_names — human-readable names for Pfam accessions
    def make_pfam_names(val, lookup):
        if not val or pd.isna(val) or str(val).strip() in ('', 'nan'):
            return ''
        names = []
        for acc in str(val).split(','):
            acc = acc.strip()
            if acc in lookup:
                names.append(lookup[acc])
        return '; '.join(names)

    master['interpro_pfam_names'] = master['interpro_pfam'].apply(
        lambda x: make_pfam_names(x, PFAM_NAMES))

    master['annotation_status'] = master.apply(assign_annotation_status, axis=1)
    n_l = (master['annotation_status']=='Light').sum()
    n_p = (master['annotation_status']=='Partial').sum()
    n_d = (master['annotation_status']=='Dark').sum()
    print(f'  annotation_status: Light={n_l:,} ({100*n_l/len(master):.1f}%), Partial={n_p:,} ({100*n_p/len(master):.1f}%), Dark={n_d:,} ({100*n_d/len(master):.1f}%)')

    n_pfam = (master['interpro_pfam_names'] != '').sum()
    print(f'  interpro_pfam_names: {n_pfam:,} genes')

    # isoelectric_point and molecular_weight_kda from protein sequence
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        def compute_pi_mw(seq):
            if not seq or pd.isna(seq) or str(seq).strip() in ('', 'nan'):
                return pd.Series({'isoelectric_point': None, 'molecular_weight_kda': None})
            try:
                s = str(seq).strip().replace('*','').replace('-','')
                s = s.replace('B','D').replace('Z','E').replace('X','A').replace('U','C')
                a = ProteinAnalysis(s)
                return pd.Series({'isoelectric_point': round(a.isoelectric_point(), 2),
                                  'molecular_weight_kda': round(a.molecular_weight() / 1000, 2)})
            except Exception:
                return pd.Series({'isoelectric_point': None, 'molecular_weight_kda': None})
        pi_mw = master['protein_sequence'].apply(compute_pi_mw)
        master['isoelectric_point']    = pi_mw['isoelectric_point']
        master['molecular_weight_kda'] = pi_mw['molecular_weight_kda']
        n_pi = master['isoelectric_point'].notna().sum()
        print(f'  isoelectric_point / MW: {n_pi:,} genes')
    except ImportError:
        print("  WARNING: BioPython not available — skipping pI/MW")
        master['isoelectric_point']    = None
        master['molecular_weight_kda'] = None

    # kegg_pathway_names — map KO numbers to pathway names via KEGG
    if not args.skip_kegg_api and 'eggnog_kegg_ko' in master.columns:
        ko_to_pathway = fetch_kegg_pathway_names(master['eggnog_kegg_ko'])
        if ko_to_pathway:
            master['kegg_pathway_names'] = map_kegg_pathways(master['eggnog_kegg_ko'], ko_to_pathway)
            n = (master['kegg_pathway_names'] != '').sum()
            print(f"  kegg_pathway_names: {n:,} genes")
        else:
            master['kegg_pathway_names'] = ''
    else:
        master['kegg_pathway_names'] = ''
        if args.skip_kegg_api:
            print("  kegg_pathway_names: skipped (--skip_kegg_api)")

    # ── Step 6: Summary ───────────────────────────────────────
    print("\n--- Classification summary ---")
    tier_order = ['Gold', 'Silver', 'Bronze', 'Informative', 'Unclassified']
    tier_counts = master['confidence_tier'].value_counts()
    print("\nConfidence tiers:")
    for tier in tier_order:
        n = tier_counts.get(tier, 0)
        pct = 100 * n / len(master)
        print(f"  {tier:<15}: {n:>6,}  ({pct:.1f}%)")

    print("\nTop 25 protein groups:")
    group_counts = master['protein_group'].value_counts()
    for group, n in group_counts.head(25).items():
        print(f"  {group:<40}: {n:>5,}")

    print("\nCellChat roles:")
    for role, n in master['cellchat_role'].value_counts().items():
        print(f"  {role:<30}: {n:>6,}")

    # ── Step 7: Write master table ────────────────────────────
    print(f"\n--- Writing master annotation table ---")
    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)
    # ── Reorder columns: move {species}_locus_tag to col 2 if present ──
    locus_col = f"{species.lower()}_locus_tag"
    if locus_col in master.columns:
        cols = list(master.columns)
        cols.remove(locus_col)
        cols.insert(1, locus_col)
        master = master[cols]

    master.to_csv(args.out_table, sep='\t', index=False)
    print(f"  Written: {args.out_table}")
    print(f"  Shape: {master.shape[0]:,} genes × {master.shape[1]} columns")
    print(f"  Columns: {list(master.columns)}")

    # ── Step 8: Export gene lists ─────────────────────────────
    print("\n--- Exporting gene lists ---")
    out_dir = Path(args.out_genelists)
    out_dir.mkdir(parents=True, exist_ok=True)

    groups = master[master['protein_group'] != 'Unclassified']['protein_group'].unique()
    n_written = 0
    for group in sorted(groups):
        genes = master[master['protein_group'] == group][gene_id_col].tolist()
        with open(out_dir / f"{group}.txt", 'w') as f:
            f.write('\n'.join(genes) + '\n')
        n_written += 1

    # Unclassified
    unclassified = master[master['protein_group'] == 'Unclassified'][gene_id_col].tolist()
    with open(out_dir / 'Unclassified.txt', 'w') as f:
        f.write('\n'.join(unclassified) + '\n')

    # CellChat lists

    # Annotation status gene lists (Light/Partial/Dark)
    for status in ['Light', 'Partial', 'Dark']:
        status_genes = master[master['annotation_status']==status][gene_id_col].tolist()
        with open(out_dir / f'{status}_genes.txt', 'w') as f:
            f.write('\n'.join(status_genes) + '\n')
        print(f'  {status}_genes: {len(status_genes):,} genes')

    for role in ['Ligand', 'Receptor', 'Receptor_candidate']:
        genes = master[master['cellchat_role'] == role][gene_id_col].tolist()
        with open(out_dir / f'CellChat_{role}.txt', 'w') as f:
            f.write('\n'.join(genes) + '\n')
        print(f"  CellChat_{role}: {len(genes):,} genes")

    print(f"\n  Gene lists written: {n_written} functional groups")
    print(f"  Gene lists directory: {out_dir}/")

    # ── Step 9: R helper script ───────────────────────────────
    r_helper = Path(args.out_table).parent / f"{species}_load_gene_lists.R"
    with open(r_helper, 'w') as f:
        f.write(f'''# ============================================================
# {species} Functional Annotation — Seurat/CellChat helper
# Generated by 08_merge_annotate.py v3
# ============================================================

library(Seurat)
library(dplyr)

GENE_LISTS_DIR <- "{out_dir}"
MASTER_TABLE   <- "{args.out_table}"

# ── Load master annotation ────────────────────────────────
annot <- read.table(MASTER_TABLE, sep="\\t", header=TRUE,
                    stringsAsFactors=FALSE, quote="")
cat("Master table:", nrow(annot), "genes x", ncol(annot), "columns\\n")

# ── Useful lookup vectors ─────────────────────────────────
# Best human gene name for each {species} gene
{species.lower()}_to_human <- setNames(annot$best_human_gene_name, annot${species.lower()}_gene_id)
{species.lower()}_to_human <- {species.lower()}_to_human[{species.lower()}_to_human != ""]

# Protein group per gene
{species.lower()}_to_group <- setNames(annot$protein_group, annot${species.lower()}_gene_id)

# CellChat role per gene
{species.lower()}_to_cellchat <- setNames(annot$cellchat_role, annot${species.lower()}_gene_id)

# ── Load gene lists ───────────────────────────────────────
list_files <- list.files(GENE_LISTS_DIR, pattern="\\\\.txt$", full.names=TRUE)
gene_lists <- lapply(list_files, readLines)
names(gene_lists) <- gsub("\\\\.txt$", "", basename(list_files))
cat("Gene lists loaded:", length(gene_lists), "\\n")
for (nm in names(gene_lists)) {{
  cat("  ", nm, ":", length(gene_lists[[nm]]), "genes\\n")
}}

# ── Score cells with AddModuleScore ──────────────────────
# Replace `seurat_obj` with your Seurat object name
# Genes present in the atlas
# present <- rownames(seurat_obj)
# gene_lists_filt <- lapply(gene_lists, function(g) g[g %in% present])
#
# for (group in names(gene_lists_filt)) {{
#   n <- length(gene_lists_filt[[group]])
#   if (n >= 5) {{
#     seurat_obj <- AddModuleScore(seurat_obj,
#                                   features = list(gene_lists_filt[[group]]),
#                                   name = paste0(group, "_score"))
#     cat("Scored:", group, "(", n, "genes)\\n")
#   }}
# }}

# ── CellChat custom DB setup ──────────────────────────────
# Ligands and receptors for custom CellChatDB construction
# ligands   <- gene_lists[["CellChat_Ligand"]]
# receptors <- gene_lists[["CellChat_Receptor"]]
# candidates <- gene_lists[["CellChat_Receptor_candidate"]]

cat("\\nDone.\\n")
''')

    print(f"\n  R helper: {r_helper}")

    # ── Step 10: Write annotation summary markdown ────────────
    summary_path = Path(args.out_table).parent / f"{species}_annotation_summary.md"

    # Compute summary stats from master table
    n_genes       = len(master)
    n_scaffolds   = master['scaffold'].nunique() if 'scaffold' in master.columns else 0
    n_eggnog      = (master['eggnog_ortholog'].fillna('') != '').sum()
    n_of          = (master['of_human_ortholog'].fillna('') != '').sum()
    n_ips         = (master['interpro_pfam'].fillna('') != '').sum()
    n_rbh         = (master['rbh_gene_name'].fillna('') != '').sum()
    n_t3          = (master['tier3_domain_hits'].fillna('') != '').sum()
    n_t4          = (master['tier4_blast_hit'].fillna('') != '').sum()
    n_tm          = int((master['tmhmm_n_tmrs'] >= 1).sum()) if 'tmhmm_n_tmrs' in master.columns else 0
    n_sp          = int((master['signalp_prediction'] == 'SP').sum())
    n_human       = (master['best_human_gene_name'].fillna('') != '').sum()
    n_go          = (master['all_go_terms'].fillna('') != '').sum()
    n_kegg        = (master['kegg_pathway_names'].fillna('') != '').sum()
    med_len       = int(master['protein_length_aa'].median()) if 'protein_length_aa' in master.columns else 0

    tier_counts   = master['confidence_tier'].value_counts()
    status_counts = master['annotation_status'].value_counts()
    group_counts  = master['protein_group'].value_counts()
    role_counts   = master['cellchat_role'].value_counts()

    n_gold        = tier_counts.get('Gold', 0)
    n_silver      = tier_counts.get('Silver', 0)
    n_bronze      = tier_counts.get('Bronze', 0)
    n_info        = tier_counts.get('Informative', 0)
    n_unclass     = tier_counts.get('Unclassified', 0)

    n_light       = status_counts.get('Light', 0)
    n_partial     = status_counts.get('Partial', 0)
    n_dark        = status_counts.get('Dark', 0)

    n_ligand      = role_counts.get('Ligand', 0)
    n_receptor    = role_counts.get('Receptor', 0)
    n_rec_cand    = role_counts.get('Receptor_candidate', 0)
    n_neither     = role_counts.get('Neither', 0)
    n_signalling  = n_ligand + n_receptor + n_rec_cand

    n_sec         = int(master['protein_group_secondary'].notna().sum()) if 'protein_group_secondary' in master.columns else 0

    proteome_file = Path(args.proteome_fasta).name if args.proteome_fasta else ''
    gff3_file     = Path(args.gff3).name if args.gff3 else ''
    run_date      = Path(args.out_table).parent.parent.name.split('_')[-1] if '_' in Path(args.out_table).parent.parent.name else ''
    genome_ver    = args.genome_version if args.genome_version else f'{species}_v1.0'
    prot_source   = args.proteome_source if args.proteome_source else ''

    # Build protein groups table (all groups, sorted by count descending)
    group_rows = ''
    for grp, cnt in group_counts.items():
        group_rows += f'| {grp} | {cnt:,} |\n'

    import datetime
    today = datetime.date.today().isoformat()

    md = f"""# *{species}* Annotation Summary

## Contents

- [Run Details](#run-details)
- [Gene Universe](#gene-universe)
- [Confidence Tiers](#confidence-tiers)
- [Annotation Status](#annotation-status)
- [Protein Groups](#protein-groups)
- [CellChat Roles](#cellchat-roles)
- [Tool Versions Used](#tool-versions-used)
- [Gene Lists](#gene-lists)

---

## Run Details

| Field | Value |
|-------|-------|
| Genome version | {genome_ver} |
| Proteome source | {prot_source} |
| Proteome file | `{proteome_file}` |
| GFF3 file | `{gff3_file}` |
| Pipeline run date | {today} |
| Merge script | `08_merge_annotate.py` |
| Master table | `{species}_master_annotation.tsv` ({n_genes:,} genes x {master.shape[1]} columns) |

See [`docs/column_descriptions.md`](../../docs/column_descriptions.md) for full column schema and [`docs/classification_system.md`](../../docs/classification_system.md) for classification logic and protein group definitions.

---

## Gene Universe

| Metric | Value |
|--------|-------|
| Total genes | {n_genes:,} |
| Scaffolds | {n_scaffolds:,} |
| Genes with eggNOG annotation | {n_eggnog:,} ({100*n_eggnog/n_genes:.1f}%) |
| Genes with OrthoFinder ortholog | {n_of:,} ({100*n_of/n_genes:.1f}%) |
| Genes with InterProScan annotation | {n_ips:,} ({100*n_ips/n_genes:.1f}%) |
| Genes with RBH SwissProt hit | {n_rbh:,} ({100*n_rbh/n_genes:.1f}%) |
| Genes with Tier 3 HMM domain hit | {n_t3:,} ({100*n_t3/n_genes:.1f}%) |
| Genes with Tier 4 BLAST hit | {n_t4:,} ({100*n_t4/n_genes:.1f}%) |
| Genes with >=1 TM helix (DeepTMHMM) | {n_tm:,} ({100*n_tm/n_genes:.1f}%) |
| Genes with signal peptide SP (SignalP) | {n_sp:,} ({100*n_sp/n_genes:.1f}%) |
| Genes with best human gene name | {n_human:,} ({100*n_human/n_genes:.1f}%) |
| Genes with GO terms (any source) | {n_go:,} ({100*n_go/n_genes:.1f}%) |
| Genes with KEGG pathway names | {n_kegg:,} ({100*n_kegg/n_genes:.1f}%) |
| Median protein length | {med_len} aa |

---

## Confidence Tiers

Evidence strength from orthology and homology sources. See [`docs/classification_system.md`](../../docs/classification_system.md) for tier definitions.

| Tier | n genes | % |
|------|---------|---|
| Gold | {n_gold:,} | {100*n_gold/n_genes:.1f}% |
| Silver | {n_silver:,} | {100*n_silver/n_genes:.1f}% |
| Bronze | {n_bronze:,} | {100*n_bronze/n_genes:.1f}% |
| Informative | {n_info:,} | {100*n_info/n_genes:.1f}% |
| Unclassified | {n_unclass:,} | {100*n_unclass/n_genes:.1f}% |

---

## Annotation Status

Functional characterisation in public databases, following Stephens et al. 2026 (*Genome Biol Evol* evag072).

| Status | n genes | % | Definition |
|--------|---------|---|------------|
| Light | {n_light:,} | {100*n_light/n_genes:.1f}% | Named human ortholog or informative eggNOG description |
| Partial | {n_partial:,} | {100*n_partial/n_genes:.1f}% | Named Pfam domain; no ortholog or functional description |
| Dark | {n_dark:,} | {100*n_dark/n_genes:.1f}% | No ortholog, no description, no named domain |

---

## Protein Groups

24 functional groups assigned based on orthology, domain, and HMM evidence. For group definitions and assignment criteria see [`docs/classification_system.md`](../../docs/classification_system.md).

| Protein group | n genes |
|---------------|---------|
{group_rows}
Dual-role genes with a secondary group assignment (`protein_group_secondary`): **{n_sec:,} genes**

---

## CellChat Roles

Predicted intercellular communication roles based on signal peptide and transmembrane topology.

| Role | n genes | % |
|------|---------|---|
| Neither | {n_neither:,} | {100*n_neither/n_genes:.1f}% |
| Receptor_candidate | {n_rec_cand:,} | {100*n_rec_cand/n_genes:.1f}% |
| Receptor | {n_receptor:,} | {100*n_receptor/n_genes:.1f}% |
| Ligand | {n_ligand:,} | {100*n_ligand/n_genes:.1f}% |

Total putative signalling genes (Ligand + Receptor + Receptor_candidate): **{n_signalling:,} ({100*n_signalling/n_genes:.1f}%)**

---

## Tool Versions Used

| Tool | Version | Notes |
|------|---------|-------|
| eggNOG-mapper | v2 | eggNOG 5.0 database |
| OrthoFinder | v2 | MSA mode, 4 species: {species}, human, mouse, *N. vectensis* |
| InterProScan | 5.78-109.0 | 17 member databases |
| DIAMOND | v2.2.1 | RBH vs SwissProt (Jan 2026) |
| SignalP | 6.0 | slow-sequential mode |
| DeepTMHMM | 1.0 | 2,000-sequence chunks |
| HMMER | 3.4 | 179 Tier 3 profiles |
| BioPython | v1.87 | pI and MW computation |
| 08_merge_annotate.py | | Integrates all 7 annotation layers; assigns confidence tiers, protein groups, CellChat roles, pI/MW, KEGG pathway names |

---

## Gene Lists

Pre-computed gene lists are in the `gene_lists/` directory:

- `Light_genes.txt` ({n_light:,} genes)
- `Partial_genes.txt` ({n_partial:,} genes)
- `Dark_genes.txt` ({n_dark:,} genes)
- `CellChat_Ligand.txt` ({n_ligand:,} genes)
- `CellChat_Receptor.txt` ({n_receptor:,} genes)
- `CellChat_Receptor_candidate.txt` ({n_rec_cand:,} genes)
- `by_protein_group/` (one file per functional group; 24 files total)

An R helper script (`{species}_load_gene_lists.R`) is also provided for loading gene lists directly into Seurat or other R workflows.
"""

    with open(summary_path, 'w') as f:
        f.write(md)
    print(f"\n  Annotation summary: {summary_path}")

    print(f"\n{'='*60}")
    print(f"  Pipeline complete for {species}")
    print(f"  Master table : {args.out_table}")
    print(f"  Gene lists   : {out_dir}/")
    print(f"  R helper     : {r_helper}")
    print(f"  Summary MD   : {summary_path}")
    print(f"  Merge log    : {log_path}")
    print(f"{'='*60}\n")

    # ── Step 11: Stage outputs to GitHub upload directory ─────
    if args.staging_dir:
        import shutil
        print(f"\n--- Staging outputs to {args.staging_dir} ---")
        stage = Path(args.staging_dir) / 'species' / species
        stage_int = stage / 'intermediate'
        stage_of  = stage_int / 'orthofinder' / 'Orthologues'
        stage_gl  = stage / 'gene_lists' / 'by_protein_group'

        for d in [stage_int, stage_of, stage_gl]:
            d.mkdir(parents=True, exist_ok=True)

        # Top-level files
        shutil.copy2(args.out_table,  stage / f"{species}_master_annotation.tsv")
        shutil.copy2(str(r_helper),   stage / f"{species}_load_gene_lists.R")
        shutil.copy2(str(log_path),   stage / f"{species}_merge_run.out")
        shutil.copy2(str(summary_path), stage / f"{species}_annotation_summary.md")
        print(f"  Copied master table, R helper, merge log, summary MD")

        # Gene lists
        if out_dir.exists():
            for f in out_dir.iterdir():
                if f.is_file():
                    shutil.copy2(str(f), stage / 'gene_lists' / f.name)
                elif f.is_dir() and f.name == 'by_protein_group':
                    for gf in f.iterdir():
                        shutil.copy2(str(gf), stage_gl / gf.name)
        print(f"  Copied gene lists")

        # Intermediate: eggNOG
        if args.eggnog and Path(args.eggnog).exists():
            shutil.copy2(args.eggnog, stage_int / 'eggnog_results.tsv')
            print(f"  Copied eggnog_results.tsv")

        # Intermediate: IPS filtered to Pfam + PANTHER, cols 1,4,5,6,14
        if args.interproscan and Path(args.interproscan).exists():
            try:
                print(f"  Filtering InterProScan to Pfam + PANTHER...")
                ips_chunks = []
                for chunk in pd.read_csv(args.interproscan, sep='\t', header=None,
                                         comment='#', chunksize=100000,
                                         dtype=str, low_memory=False):
                    chunk = chunk[chunk[3].isin(['Pfam', 'PANTHER'])]
                    ips_chunks.append(chunk[[0, 3, 4, 5, 13]])
                if ips_chunks:
                    ips_filtered = pd.concat(ips_chunks, ignore_index=True)
                    ips_filtered.to_csv(stage_int / 'interproscan_results.tsv',
                                        sep='\t', index=False, header=False)
                    print(f"  Copied interproscan_results.tsv ({len(ips_filtered):,} rows, Pfam+PANTHER only)")
                else:
                    print(f"  WARNING: No Pfam/PANTHER rows found in InterProScan output")
            except Exception as e:
                print(f"  WARNING: IPS filtering failed ({e}) — skipping")

        # Intermediate: SignalP, Tier3, Tier4, DeepTMHMM
        for src, dst in [
            (args.signalp,                        stage_int / 'signalp_summary.txt'),
            (args.tier3,                          stage_int / 'tier3_hmm_hits.domtblout'),
            (args.tier4,                          stage_int / 'tier4_blast_hits.tsv'),
            (args.tmhmm if args.tmhmm else '',    stage_int / 'deeptmhmm_merged.gff3'),
        ]:
            if src and Path(src).exists():
                shutil.copy2(src, dst)
                print(f"  Copied {dst.name}")
            elif src:
                print(f"  WARNING: {src} not found — skipping {dst.name}")

        # Intermediate: OrthoFinder
        of_results = Path(args.orthofinder)
        if of_results.exists():
            orthogroups = of_results / 'Orthogroups' / 'Orthogroups.tsv'
            stats       = of_results / 'Comparative_Genomics_Statistics' / 'Statistics_Overall.tsv'
            orthologues = of_results / 'Orthologues'
            if orthogroups.exists():
                shutil.copy2(str(orthogroups), stage_int / 'orthofinder' / 'Orthogroups.tsv')
                print(f"  Copied orthofinder/Orthogroups.tsv")
            if stats.exists():
                shutil.copy2(str(stats), stage_int / 'orthofinder' / 'Statistics_Overall.tsv')
                print(f"  Copied orthofinder/Statistics_Overall.tsv")
            if orthologues.exists():
                dst_orth = stage_int / 'orthofinder' / 'Orthologues'
                if dst_orth.exists():
                    shutil.rmtree(dst_orth)
                shutil.copytree(str(orthologues), str(dst_orth))
                print(f"  Copied orthofinder/Orthologues/")

        print(f"\n  Staging complete: {stage}")
        print(f"  To download locally, run:")
        print(f"  scp -r <user>@pegasus.ccs.miami.edu:{stage} ~/Desktop/coral-func-annotation-pipeline/species/")
        print()

    # Close log file
    _log_fh.close()
    import sys as _sys
    _sys.stdout = _sys.__stdout__


if __name__ == '__main__':
    main()
