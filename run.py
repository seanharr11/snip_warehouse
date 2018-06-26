import os
from snip_warehouse import SnipLoader
from snip_warehouse.schema import init_db

DB_NAME = os.environ["SNIP_DB_NAME"]


input("Are you sure you want to re-init DB?"
      "This takes 2-3 hrs per chromosome")
init_db(DB_NAME)
snip_loader = SnipLoader(DB_NAME)
chr_suffixes = [i for i in range(1, 23)]
chr_suffixes += ["X", "Y", "MT"]
for chr_suffix in chr_suffixes:
    snip_loader.download_dbsnp_file(f"refsnp-chr{chr_suffix}.json.gz",
                                    chr_suffix)
    snip_loader.load_ref_snps(
        f"refsnp-chr{chr_suffix}.json.gz", str(chr_suffix))
    os.system(f"rm refsnp-chr{chr_suffix}.json")
