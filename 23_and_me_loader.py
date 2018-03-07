import csv
from sqlalchemy import create_engine


with open("data/23_and_me_snps.txt", "r") as fp:
    cnt = 1
    variants = 0
    engine = create_engine("sqlite:///data/snps.sql")
    dr = csv.DictReader(filter(lambda r: r[0] != '#', fp), delimiter='\t')
    # dr.fieldnames = ('rsid', 'chromosome', 'position' 'genotype',)
    fp = open("my_diseases.csv", "w")
    for snp in dr:
        if cnt % 1000 == 0:
            print(f"Count: {cnt}")
        if snp['rsid'][0] == "i":
            # 23andMe's internal Snp
            continue
        rsid = snp['rsid'].replace("rs", "")
        q = f"""
            SELECT (snv_freq.allele_count / snv_freq.total_count) freq,
                   scdn.disease_name_csv,
                   scdn.clinical_significance_csv,
                   snvs.rsnp_id,
                   UPPER(alt_seq) alt_seq,
                   scdn.citation_csv
            FROM snvs
            INNER JOIN snv_clinical_disease_names scdn
             ON scdn.snv_id = snvs.id
            LEFT JOIN snv_frequencies snv_freq
             ON snv_freq.snv_id = snvs.id
            WHERE rsnp_id = '{rsid}'
             AND '{snp['genotype']}' LIKE '%' || UPPER(alt_seq) || '%'
             AND disease_name_csv != ''
             AND disease_name_csv != 'not specified'
        """
        # print(q)
        rows = engine.execute(q).fetchall()
        """INNER JOIN gene_snvs
             ON gene_snvs.snv_id = snvs.id
            INNER JOIN genes
             ON UPPER(gene_snvs.locus) = UPPER(genes.locus)
        """
        if rows:
            variants += 1
            for row in rows:
                fp.write("|".join([
                    '1' if row.alt_seq * 2 == snp['genotype'] else '0',
                    str(row),
                    f"<{snp['genotype']}>\n\r"]))
        cnt += 1
    print(f"Variant Count: '{variants}'")
