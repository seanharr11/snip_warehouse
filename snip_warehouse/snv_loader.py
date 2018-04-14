import aiofiles
import asyncio
import asyncpg
from collections import namedtuple
import ftplib
from multiprocessing import cpu_count, Value, Process
import os
import time
from sqlalchemy import create_engine
import threading
import ujson as json
import zlib

# from .schema import (
#     Snv, SnvFrequency, SnvClinicalDiseaseName, GeneSnv)

FrequencyStudy = namedtuple('FrequencyStudy', [
    'name',
    'allele_count',
    'total_count'])
# https://www.ncbi.nlm.nih.gov/books/NBK21088/table/ch5.ch5_t3/?report=objectonly


class SnvLoader:
    def __init__(self, database_name):  # , db_conn_string):
        self._variant_output_filename = "./variants.csv"
        self.rows_processed = Value('i', 0)
        self.new_ref_snp_allele_id = None
        self.database_name = database_name
        self.file_blocksize = 1024 * 1024
        # self.db_conn_string = db_conn_string

    def download_dbsnp_file(self, dbsnp_filename):
        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.set_debuglevel(2)
        ftp.login()
        ftp.cwd("snp/.redesign/latest_release/JSON")
        size = ftp.size(dbsnp_filename)
        print(f"Filesize: {round(size / 1024 / 1024 / 1024, 2)} GB")
        sock = None
        while not sock:
            print("Trying to establish FTP conn")
            sock = ftp.transfercmd(f"RETR {dbsnp_filename}")
            time.sleep(5)

        def download():
            decompressor = zlib.decompressobj(32 + zlib.MAX_WBITS)
            transferred = 0
            fp = open(dbsnp_filename.replace(".gz", ""), "wb")
            blocks = 0
            while True:
                byte_chunk = sock.recv(1024*1024*8)
                if not byte_chunk:
                    break
                blocks += 1
                transferred = transferred + len(byte_chunk)
                transferred_mb = round(transferred / 1024 / 1024, 2)
                if blocks % 1000 == 0:
                    print(
                        f"Transferred {transferred_mb}MB / "
                        f"{round(size / 1024 / 1024, 2)}MB")
                fp.write(decompressor.decompress(byte_chunk))
            sock.close()
        t = threading.Thread(target=download)
        t.start()
        while t.is_alive():
            t.join(60)
            ftp.voidcmd("NOOP")

    def _print_status(self):
        if self.rows_processed.value % 10000 == 0:
            print(f"Processed '{self.rows_processed.value}' Ref SNPs")
        self.rows_processed.value += 1

    def load_ref_snps(self, dbsnp_filename):
        num_processes = cpu_count()
        print(f"Found '{num_processes}' CPUs")
        conn = create_engine(os.environ['SNVS_DB_URL'])
        row = conn.execute(
            "SELECT MAX(id) FROM ref_snp_alleles").fetchone()
        self.new_ref_snp_allele_id = Value('i', (row[0] or 0) + 1)
        start_offset = 0
        end_offset = 0
        file_size = os.path.getsize(dbsnp_filename)
        processes = []
        for i in range(num_processes):
            with open(dbsnp_filename, "r") as fp:
                end_offset = self._find_newline_offset(
                    fp, end_offset + file_size // num_processes)
            p = Process(target=self._load_ref_snps,
                        args=(dbsnp_filename, start_offset, end_offset))
            start_offset = end_offset
            print(f"starting process '{i}'")
            p.start()
            processes.append(p)
        for p in processes:
            p.join()

    @staticmethod
    def _find_newline_offset(fp, line_offset):
        if line_offset == 0:
            return fp
        else:
            fp.seek(line_offset)
            fp.readline()
            return fp.tell()

    def _load_ref_snps(self, dbsnp_filename, start_offset, end_offset):
        print(f"Handling bytes '{start_offset}' => '{end_offset}'")
        loop = asyncio.get_event_loop()
        loop.run_until_complete(
            self._async_load_ref_snps(
                dbsnp_filename, start_offset, end_offset, loop))

    async def _async_load_ref_snps(self, dbsnp_filename,
                                   start_offset, end_offset, loop):
        dbsnp_fp = await aiofiles.open(dbsnp_filename, "r", loop=loop)
        conn = await asyncpg.connect(
            user='SeanH', database=self.database_name, loop=loop)
        await conn.execute(f"SET session_replication_role = replica")
        await dbsnp_fp.seek(start_offset)
        json_ls, offset = await self.read_json(dbsnp_fp, self.file_blocksize)
        insert_stmts = []
        while json_ls and offset < end_offset:
            futures = [self.generate_insert_stmts(json_ls),
                       self.read_json(dbsnp_fp, self.file_blocksize),
                       self.insert_records(insert_stmts, conn)]
            insert_stmts, (json_ls, offset), _ = await asyncio.gather(*futures)
        print(f"Bytes '{start_offset}' -> '{end_offset}' parsed!")
        await dbsnp_fp.close()
        await conn.close()

    @staticmethod
    async def read_json(dbsnp_fp, blocksize):
        # TODO: Couldn't blocksize be "num_lines" ?
        init_tell = tell = prev_tell = await dbsnp_fp.tell()
        raw_json = ""
        while True:
            line = await dbsnp_fp.readline()
            prev_tell = tell
            tell = await dbsnp_fp.tell()
            if not line:
                break  # EOF
            if (tell - init_tell) > blocksize:
                # Went past our blocksize...
                await dbsnp_fp.seek(prev_tell)
                tell = prev_tell
                break
            raw_json += (line + ",")
        final_json = raw_json[:-1] if raw_json else ""
        return json.loads("[" + final_json + "]"), tell

    async def generate_insert_stmts(self, json_ls):
        (ref_snp_alleles, gene_ref_snp_alleles,
         snv_freqs, snv_clinical_diseases) = [], [], [], []
        if not json_ls:
            return []
        for rsnp_json in json_ls:
            self._print_status()
            rsnp_placements = rsnp_json['primary_snapshot_data'][
                                    'placements_with_allele']
            ref_snp_id = int(rsnp_json['refsnp_id'])
            if not rsnp_placements:
                continue
            alleles = self.find_alleles_from_assembly(rsnp_placements)
            if not alleles:
                continue
            variant_alleles = self.get_variant_allele(alleles)
            if not variant_alleles:
                continue
            for variant_allele in variant_alleles:
                # ref_snp_alleles.id is auto-incrementing, follow it!
                with self.new_ref_snp_allele_id.get_lock():
                    self.new_ref_snp_allele_id.value += 1
                    new_ref_snp_allele_id = self.new_ref_snp_allele_id.value
                allele_annotation = self.get_allele_annotation(
                    rsnp_json, variant_allele['allele_idx'])
                ref_snp_alleles.append((new_ref_snp_allele_id,
                                        ref_snp_id,
                                        variant_allele['del_seq'],
                                        variant_allele['ins_seq'],
                                        variant_allele['position'],))
                gene_ref_snp_alleles += [(new_ref_snp_allele_id,
                                          gene['locus'],
                                          gene['id'],)
                                         for gene in allele_annotation[
                                             'genes']]
                snv_freqs += [(new_ref_snp_allele_id,
                               fs.name,
                               fs.allele_count,
                               fs.total_count)
                              for fs in allele_annotation[
                                  'frequency_studies']]
                snv_clinical_diseases += [
                    (new_ref_snp_allele_id,
                     clin['disease_names'],
                     clin['clinical_significances'],
                     clin['citation_list'],)
                    for clin in allele_annotation['clinical_entries']]
        insert_queries = [
            ('ref_snp_alleles',
                ('id', 'ref_snp_id', 'del_seq', 'ins_seq', 'position'),
                ref_snp_alleles),
            ('gene_ref_snp_alleles',
                ('ref_snp_allele_id', 'locus', 'gene_id'),
                gene_ref_snp_alleles),
            ('ref_snp_allele_freq_studies',
             ('ref_snp_allele_id', 'name', 'allele_count',
              'total_count'),
             snv_freqs),
            ('ref_snp_allele_clin_diseases',
                ('ref_snp_allele_id', 'disease_name_csv',
                 'clinical_significance_csv', 'citation_list'),
                snv_clinical_diseases)]
        return insert_queries

    async def insert_records(self, insert_stmts, conn):
        # futures = []
        for j, (table_name, columns, records) in enumerate(insert_stmts):
            """futures.append(connections[j].copy_records_to_table(
                           table_name,
                           columns=columns,
                           records=records))
            """
            await conn.copy_records_to_table(
                               table_name,
                               columns=columns,
                               records=records)
        # await asyncio.gather(*futures)

    def find_alleles_from_assembly(self,
                                   rsnp_placements,
                                   assembly_name="GRCh38"):
        for rsnp_placement in rsnp_placements:
            annot = rsnp_placement.get('placement_annot')
            if not annot or not annot.get('seq_id_traits_by_assembly'):
                return
            assembly_info_ls = annot['seq_id_traits_by_assembly']
            if len(assembly_info_ls) > 1:
                print(f"Assembly Info ls len g.t. 1: {assembly_info_ls}")
            assembly_info = assembly_info_ls[0]
            # TODO: Why is this a list
            this_assembly_name = assembly_info.get("assembly_name") or None
            if assembly_name in this_assembly_name:
                alleles = rsnp_placement['alleles']
                return alleles

    def get_variant_allele(self, alleles):
        # Find the allele that represents the variation
        variant_allele_tups = [(a, i) for (i, a) in enumerate(alleles)
                               if a['allele']['spdi']['inserted_sequence'] !=
                               a['allele']['spdi']['deleted_sequence']]
        if not variant_allele_tups:
            return {}
        var_alleles = []
        for allele_tup in variant_allele_tups:
            var_spdi = allele_tup[0]['allele']['spdi']
            allele_idx = allele_tup[1]
            var_alleles.append({
                'del_seq': var_spdi['deleted_sequence'],
                'ins_seq': var_spdi['inserted_sequence'],
                'position': var_spdi['position'],
                'allele_idx': allele_idx
            })
        return var_alleles

    def get_allele_annotation(self, rsnp_obj, allele_idx):
        var_allele_annotation = rsnp_obj['primary_snapshot_data'][
                                    'allele_annotations'][
                                    allele_idx]
        assembly_annot = var_allele_annotation['assembly_annotation']
        frequencies = var_allele_annotation['frequency']
        fs = [FrequencyStudy(
            name=freq['study_name'],
            allele_count=freq['allele_count'],
            total_count=freq['total_count']) for freq in frequencies if freq]
        clinical_entries = [{
            'citation_list': clin['citations'],
            'disease_names': ",".join(clin['disease_names']),
            'clinical_significances': ",".join(clin['clinical_significances'])}
                for clin in var_allele_annotation['clinical']]
        return {'frequency_studies': fs,
                'clinical_entries': clinical_entries,
                # TODO: Multiple assembly annotations?
                'genes': assembly_annot[0]['genes']
                if assembly_annot else None,
                # TODO: Just store gene 'names' and 'locus'
                'seq_id': assembly_annot[0]['seq_id']}
