import asyncio
import asyncpg
import csv
import os


class SnipUploader:
    def __init__(self, username, password, host):
        self.username = username
        self.password = password
        self.host = host

    async def connect(self):
        self.pool = await asyncpg.create_pool(
            user=self.username,
            password=self.password,
            host=self.host,
            database='snip')
        self.conn = await self.pool.acquire()

    async def upload(self, iterable_buff, user_id):
        """ `iterable_buff`: Iterable
             -  will never be loaded entirely into memory!
             `user_id`: int
             -  the id of the user uploading SNP data
        """
        # await self.conn.execute(
        #    "SET session_replication_role = replica;")
        await self._upload(iterable_buff, user_id)
        # await self.conn.execute(
        #    "SET session_replication_role = DEFAULT;")

    async def _upload(self, buff, user_id=1):
        cr = csv.reader(filter(lambda ln: not ln.startswith("#"), buff),
                        delimiter="\t")
        next(cr)
        record_gen = ((int(r[0].replace("rs", "").replace("i", "")),
                       # int(r[1]),
                       r[3],
                       user_id,) for r in cr)
        """
        await self.conn.executemany(
            "INSERT INTO user_ref_snps (ref_snp_id, genotype, user_id)"
            " VALUES($1, $2, $3)", record_gen)
        """
        await self.conn.copy_records_to_table(
            "user_ref_snps",
            columns=("ref_snp_id", "genotype", "user_id",),
            records=record_gen)


async def main():
    snv_uploader = SnipUploader(
        username=os.environ['SNIP_DB_INSERT_USERNAME'],
        password=os.environ['SNIP_DB_INSERT_PASSWORD'],
        host=os.environ['SNIP_DB_HOST']
    )
    await snv_uploader.connect()
    with open("data/23_and_me_snps.txt", "r") as fp:
        await snv_uploader.upload(fp, 1)

if __name__ == "__main__":
    asyncio.get_event_loop().run_until_complete(main())
