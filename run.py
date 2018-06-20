import asyncio
import os
from snip_warehouse import SnipLoader, SnipUploader
from snip_warehouse.schema import init_db

DB_NAME = os.environ["SNIP_DB_NAME"]


def main():
    input("Are you sure you want to re-init DB?"
          "This takes 2-3 hrs per chromosome")
    init_db(DB_NAME)
    snip_loader = SnipLoader(DB_NAME)
    # chr_suffixes = [i for i in range(1, 2)]
    # chr_suffixes += ["X", "Y", "MT"]
    for chr_suffix in [1]:  # chr_suffixes:
        snip_loader.download_dbsnp_file(f"refsnp-chr{chr_suffix}.json.gz",
                                        chr_suffix)
        snip_loader.load_ref_snps(f"refsnp-chr{chr_suffix}.json")
        os.system(f"rm refsnp-chr{chr_suffix}.json")


async def async_upload_23_and_me(loop):
    snv_uploader = SnvUploader(DB_NAME)
    await snv_uploader.connect()
    with open("data/23_and_me_snps.txt", "r") as fp:
        await snv_uploader.upload(fp, user_id=1)


main()
# loop = asyncio.get_event_loop()
# loop.run_until_complete(async_upload_23_and_me(loop))
