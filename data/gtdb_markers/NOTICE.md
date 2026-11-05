# GTDB marker HMMs

These files hold the hidden Markov models of the 120 bacterial (bac120) and 53
archaeal (ar53) GTDB marker genes: 152 TIGRFAM and 16 Pfam models. They are
copied unmodified from the GTDB-Tk R232 reference data package
(`markers/tigrfam/tigrfam.hmm`, `markers/pfam/Pfam-A.hmm` and
`markers/pfam/Pfam-A.hmm.dat`), with the two libraries gzip-compressed.
`aeon classify-db` uses them when no GTDB-Tk package is given.

They are third-party data, distributed under their own licences, not aeon's
MIT licence.

## TIGRFAM (`tigrfam.hmm.gz`)

TIGRFAMs release 15.0. Licensed under the Creative Commons
Attribution-ShareAlike 4.0 International licence
(https://creativecommons.org/licenses/by-sa/4.0/). Rights were transferred
from the J. Craig Venter Institute to the National Center for Biotechnology
Information (NCBI) in April 2018; see
https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/license_and_availability.txt.

Use of TIGRFAMs must cite:

- Haft DH, Loftus BJ, Richardson DL, Yang F, Eisen JA, Paulsen IT, White O.
  TIGRFAMs: a protein family resource for the functional identification of
  proteins. Nucleic Acids Res. 2001;29(1):41-3.
- Haft DH, Selengut JD, White O. The TIGRFAMs database of protein families.
  Nucleic Acids Res. 2003;31(1):371-3.
- Selengut JD, Haft DH, Davidsen T, Ganapathy A, Gwinn-Giglio M, Nelson WC,
  Richter AR, White O. TIGRFAMs and Genome Properties: tools for the
  assignment of molecular function and biological process in prokaryotic
  genomes. Nucleic Acids Res. 2007;35(Database issue):D260-4.
- Haft DH, Selengut JD, Richter RA, Harkins D, Basu MK, Beck E. TIGRFAMs and
  Genome Properties in 2013. Nucleic Acids Res. 2013;41(Database
  issue):D387-95.

## Pfam (`Pfam-A.hmm.gz`, `Pfam-A.hmm.dat`)

Pfam-A models PF00380.20, PF00410.20, PF00466.21, PF00687.22, PF00827.18,
PF00900.21, PF01000.27, PF01015.19, PF01025.20, PF01090.20, PF01200.19,
PF01280.21, PF02576.18, PF03726.15, PF04919.13 and PF07541.13. Pfam is
released under the Creative Commons Zero (CC0) public domain dedication
(https://creativecommons.org/publicdomain/zero/1.0/).

- Mistry J, Chuguransky S, Williams L, Qureshi M, Salazar GA, Sonnhammer ELL,
  Tosatto SCE, Paladin L, Raj S, Richardson LJ, Finn RD, Bateman A. Pfam: The
  protein families database in 2021. Nucleic Acids Res. 2021;49(D1):D412-9.

## Marker selection

The bac120 and ar53 marker sets are those of the Genome Taxonomy Database:

- Parks DH, Chuvochina M, Rinke C, Mussig AJ, Chaumeil PA, Hugenholtz P. GTDB:
  an ongoing census of bacterial and archaeal diversity through a
  phylogenetically consistent, rank normalized and complete genome-based
  taxonomy. Nucleic Acids Res. 2022;50(D1):D785-94.
