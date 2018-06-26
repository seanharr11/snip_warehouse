from collections import namedtuple


# DTO's
RefSnpCopyFromData = namedtuple("RefSnpCopyFromData", [
    "ref_snp_alleles",
    "ref_snp_allele_freq_studies",
    "ref_snp_allele_clin_diseases"])

RefSnpAllele = namedtuple("RefSnpAllele", [
    'ins_seq', 'del_seq', 'position', 'chromosome',
    'ref_snp_allele_idx', 'ref_snp_id', 'gene_locii'])


RefSnpAlleleFreqStudy = namedtuple("RefSnpAlleleFreqStudy", [
    'name', 'allele_count', 'total_count', 'ref_snp_allele_idx',
    'ref_snp_id'])

RefSnpAlleleClinDisease = namedtuple("RefSnpAlleleClinDisease", [
    'disease_name_csv', 'clinical_significance_csv', 'citation_list',
    'ref_snp_allele_idx', 'ref_snp_id'])
