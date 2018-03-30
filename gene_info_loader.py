import csv
from sqlalchemy import create_engine, insert
from snp_matcher.schema import Gene

e = create_engine("sqlite:///data/snps.sql")

with open("data/gene_info_smaller.tsv", "r") as fp:
    dr = csv.DictReader(fp, delimiter='\t')
    buff = []
    cnt = 0
    for d in dr:
        buff.append(d)
        cnt += 1
        if len(buff) == 100:
            print(f"Processed {cnt} genes")
            e.execute(insert(Gene).values(buff))
            buff = []
