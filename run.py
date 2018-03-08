import asyncio
from snp_matcher import SnvLoader

loop = asyncio.get_event_loop()


async def main():
    snp_loadr = SnvLoader("sqlite:///data/snps.sql")
    chr_suffixes = [i for i in range(1, 24)]
    chr_suffixes += ["X", "Y", "MT"]
    for chr_suffix in chr_suffixes:
        snp_loadr.download_dbsnp_file(f"refsnp-chr{chr_suffix}.json.gz")
        snp_loadr.load_ref_snps(f"refsnp-chr{chr_suffix}.json")
loop.run_until_complete(main())
