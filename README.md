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

You can run create_index.py once you have this genome to fasta file (genome2fileYellowstone). Then, you can specificy the number of mismatches between spacers. 
```
usage: create_index.py [-h] -g GENOMEFILE -o OUTPUTINDEX -m MAX_MISMATCH
```
In the example below, -m indicates number of mismatches between spacers to group them between strains. In this case, I wanted perfectly matched spacers. 
```
python2 create_index.py -g genome2fileYellowstone -o yell.index -m 0
```
The example above also spits out a fasta of consolidated spacers (consolidatedspacers.fasta) which is a set of numbered spacers shared between strains. This repository has the example fasta that is created. This example also spits out a file called yell.index. The content of this file contains a strain name <tab> spacer <tab> spacer<tab>... It has all the strains with corresponding spacer found in the consolidated spacer file. An snippet is shown below.
  
```  
NL01B_C01_23	 13 24 28 47 68 78 103 131 141 157 166 193 212 220 223 237 252 255 270 284 293 309 327 332 338 350 369 377 407 438 445 463 467 472 482 543 546 567 648 673 691 698 702 703 716 753 755 781 802 812 830 844 876 883 892 895 902 906 948 949 953 959 991 994 1021 1024 1035 1058 1071 1075 1109 1114 1115 1136 1137 1154 1168 1180 1244 1255 1264 1267 1287 1302 1312 1328 1338 1344 1346 1355 1359 1364 1372 1387 1389 1422 1428 1434 1439 1477 1495 1580 1590 1695 1705 1708
NL01B_C01_14	 5 6 9 40 49 84 107 109 114 118 134 140 146 158 161 181 184 203 217 227 230 233 243 272 282 301 308 312 333 347 366 397 398 403 409 422 423 504 527 552 558 575 595 612 614 652 663 683 685 689 692 693 769 778 813 818 823 870 901 911 926 942 979 1005 1011 1041 1053 1065 1076 1080 1093 1112 1126 1128 1151 1157 1165 1171 1194 1209 1231 1242 1270 1275 1279 1283 1292 1298 1300 1305 1339 1340 1341 1350 1353 1369 1373 1380 1390 1393 1395 1396 1397 1405 1414 1426 1440 1450 1451 1460 1462 1463 1465 1482 1488 1504 1526 1527 1535 1567 1577 1585 1671 1692 1711
NL03_C02_11	 89 104 120 125 134 139 149 179 191 192 208 269 271 302 303 308 311 313 333 337 347 357 392 396 403 405 410 414 417 442 443 458 483 486 505 510 513 514 518 520 532 548 555 557 558 563 564 571 577 579 581 590 593 604 607 623 625 637 641 662 664 665 679 683 695 729 743 748 759 763 768 776 782 800 801 804 806 809 821 826 831 834 845 848 866 874 888 890 891 928 935 962 1009 1010 1029 1037 1065 1073 1091 1097 1100 1145 1148 1205 1232 1245 1246 1262 1272 1276 1288 1292 1297 1310 1322 1324 1340 1353 1362 1366 1404 1408 1479 1486 1494 1511 1529 1534 1553 1565 1574 1593 1612 1614 1625 1644 1669 1676 1679 1685
```
An example consolidated spacers fasta snippet is shown below. Each of these fasta headers match the spacers in the index file that each genome has.

```
>1
AAATTAGTTAAGGCACGGCTTAAAAGTGATCTTAAAAAT
>2
ACGGCGTAAAGTTGTAAATGCCAGGCTGGAAGTACATGTTG
>3
AAATCTAAGTATTTTTCTTCCTCCTTAAGCTCATATTTATCAGT
>4
TTCAGGTACAGGTTTAAGATAGTATGAGAAGAAATAAAG
>5
ACTTGTACTAAGATATGATCGGATGGCGATGTTATCATA
>6
TCCGTGGTGCCCAACATGTATTACGGTTACAGTGTGTTCAA
>7
ATGATAAGAAATAACCTATCTTCCCAAGTTTGAAGTTGA
>8
CAACATTGTATATTCCTTTATTATATATTATGATATTATAGCCA
>9
TTTTCAGTTTCTATATTAATTGTATTTGCTGAGGTATAAT
>10
GCAAATCAACTTCTAACCCAACTTCTTGCAAATTTTCAGC
>11
GTTTAAAACTTTTAGTTACTAAACATATTCGTAATTTATTTAAAT
>12
TATTTTCTAGTAGCAATTAGAAACACTATCACTATAATCA
>13
ACAGTTGTTTACAATTTATCACATTTTGTCTAATAGTGCTT
>14
TGAAGGGGAACCTACTTTATATATGTGTGCTATTTCT

```
