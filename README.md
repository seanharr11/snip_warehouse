# Warehousing DbSNP's JSON Data into PostgreSQL

### Intro
NCBI hosts a large, open-sourced dataset of human SNPs (Single-nucleotide Polymorphisms). Further, they store a good deal of auxillary data that is related to each SNP. The data is hosted on an FTP server here:

> ftp://ftp.ncbi.nlm.nih.gov/snp/.redesign/latest_release/JSON

and is split across 25 gzipped JSON files (Chromosomes 1-22, X, Y and Mitochondrial DNA), amassing a total compressed size of **~100GB** (**~2TB uncompressed!**).

### Further Reading
More details can be found in [this series of blog posts](http://localhost:4000/warehousing-DbSNP-Part-I-download-and-create-db), detailing a three-part walkthrough, breaking the development of this application down in three steps:
1. [Downloading JSON SNP Data & Initilizing the Database](http://thelaziestprogrammer.com/warehousing-DbSNP-Part-I-download-and-create-db)
2. [Extracting ClinVar Disease & Frequency Study Data](http://thelaziestprogrammer.com/warehousing-DbSNP-Part-II-parsing-RefSNP-JSON)
3. [Efficiently Writing Data to PostgreSQL Database](http://thelaziestprogrammer.com/warehousing-DbSNP-Part-III-bulk-inserting-SNP-data)
