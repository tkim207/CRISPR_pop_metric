# PDI_PI_IDI
Scripts that make PDI IDI calculations

## Getting Started

These instructions will get let you caculate PI, PDI, and IDI based on a set of genomes/

### Prerequisites

BLASTn

### Installing

You will need to download two scripts create_index.py and run_PDI_total.py. They will have to be in the same folder to work. Also, you have to run these programs within the same directory.

## Running the create index

Strains with CRISPRs contain spacers that are shared between one another. However, in order to find which spacers are shared, you need to group them based on similarity and assign them a number or identifier that is shared between all strains. 
To do this, you will need a file that is tab delimited and has genome name you want to give a strain with CRISPRs and filename of spacers in fasta format. The fasta files can be located anywhere; you just have to type the directory path as well. The spacers do not have to be in order.

For example genome2fileYellowstone file contains:
```
NL01B_C01_01	NL01B_C01_01.fa
NL01B_C01_03	NL01B_C01_03.fa
NL01B_C01_04	NL01B_C01_04.fa
NL01B_C01_05	NL01B_C01_05.fa
NL01B_C01_06	NL01B_C01_06.fa
NL01B_C01_07	NL01B_C01_07.fa
NL01B_C01_08	NL01B_C01_08.fa
NL01B_C01_09	NL01B_C01_09.fa
NL01B_C01_10	NL01B_C01_10.fa
NL01B_C01_11	NL01B_C01_11.fa
NL01B_C01_12	NL01B_C01_12.fa
NL01B_C01_13	NL01B_C01_13.fa
NL01B_C01_14	NL01B_C01_14.fa
NL01B_C01_15	NL01B_C01_15.fa
NL01B_C01_16	NL01B_C01_16.fa
NL01B_C01_17	NL01B_C01_17.fa
NL01B_C01_18	NL01B_C01_18.fa
NL01B_C01_19	NL01B_C01_19.fa
NL01B_C01_20	NL01B_C01_20.fa
NL01B_C01_21	NL01B_C01_21.fa
NL01B_C01_22	NL01B_C01_22.fa
NL01B_C01_23	NL01B_C01_23.fa
NL01B_C01_24	NL01B_C01_24.fa
NL01B_C01_25	NL01B_C01_25.fa
NL03_C02_01	NL03_C02_01.fa
NL03_C02_02	NL03_C02_02.fa
NL03_C02_03	NL03_C02_03.fa
NL03_C02_04	NL03_C02_04.fa
NL03_C02_05	NL03_C02_05.fa
NL03_C02_06	NL03_C02_06.fa
NL03_C02_07	NL03_C02_07.fa
NL03_C02_08	NL03_C02_08.fa
NL03_C02_09	NL03_C02_09.fa
NL03_C02_10	NL03_C02_10.fa
NL03_C02_11	NL03_C02_11.fa
NL03_C02_12	NL03_C02_12.fa
NL13_C01_01	NL13_C01_01.fa
NL13_C01_02	NL13_C01_02.fa
NL13_C01_03	NL13_C01_03.fa
NL13_C01_04	NL13_C01_04.fa
```
This file is a list of strains of S. islandicus found in yellowstone with CRISPR spacers.

You can run create_index.py once you have this genome to fasta file (genome2fileYellowstone). Then, you can specificy the number of mismatches between spacers. In the example below, -m indicates number of mismatches between spacers to group them between strains. In this case, I wanted perfectly matched spacers.

```
python2 create_index.py -g genome2fileYellowstone -o yellowout -m 0
```

