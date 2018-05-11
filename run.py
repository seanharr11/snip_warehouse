import asyncio
import os
from snip_warehouse import SnvLoader, SnvUploader
from snip_warehouse.schema import init_db

DB_NAME = os.environ["SNVS_DB_NAME"]


def main():
    input("Are you sure you want to re-init DB?"
          "This takes 2-3 hrs per chromosome")
    # init_db(DB_NAME)
    snp_loadr = SnvLoader(DB_NAME)
    chr_suffixes = [i for i in range(24, 24)]
    chr_suffixes += ["X", "Y", "MT"]
    for chr_suffix in chr_suffixes:
        snp_loadr.download_dbsnp_file(f"refsnp-chr{chr_suffix}.json.gz")
        snp_loadr.load_ref_snps(f"refsnp-chr{chr_suffix}.json")
        os.system(f"rm refsnp-chr{chr_suffix}.json")


async def async_upload_23_and_me(loop):
    snv_uploader = SnvUploader(DB_NAME)
    await snv_uploader.connect()
    with open("data/23_and_me_snps.txt", "r") as fp:
        await snv_uploader.upload(fp, user_id=1)


main()
loop = asyncio.get_event_loop()
loop.run_until_complete(async_upload_23_and_me(loop))
