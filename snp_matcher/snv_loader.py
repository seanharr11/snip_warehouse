# import asyncio
import asyncpg
from collections import namedtuple
import ftplib
import ujson as json
import sys
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
        self.rows_processed = 0
        self.database_name = database_name
        # self.db_conn_string = db_conn_string

    def download_dbsnp_file(self, dbsnp_filename):
        decompressor = zlib.decompressobj(32 + zlib.MAX_WBITS)
        transferred = 0
        fp = open(dbsnp_filename.replace(".gz", ""), "wb")
        blocksize = 8192
        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()
        ftp.cwd("snp/.redesign/latest_release/JSON")
        size = ftp.size(dbsnp_filename)

        def callback(byte_chunk):
            nonlocal transferred
            transferred = transferred + sys.getsizeof(byte_chunk)
            transferred_mb = transferred / 1024 / 1024
            if transferred_mb % 10 == 0:
                print(
                    f"Transferred {transferred_mb}MB /"
                    "{round(size / 1024 / 1024, 2)}MB")
            fp.write(decompressor.decompress(byte_chunk))

        print(f"Filesize: {size / 1024 / 1024 / 1024}GB")
        ftp.retrbinary(f"RETR {dbsnp_filename}", callback, blocksize=blocksize)

    def _print_status(self):
        if self.rows_processed % 10000 == 0:
            print(f"Processed '{self.rows_processed}' Ref SNPs")
        self.rows_processed += 1

    async def load_ref_snps(self, dbsnp_filename):
        self.ref_variant_generator = open(dbsnp_filename, "r")
        pool = await asyncpg.create_pool(
            user='SeanH', database=self.database_name)
        conn = await pool.acquire()
        conn.execute(f"SET session_replication_role = replica")
        row = await conn.fetchrow(
            "SELECT MAX(id) FROM ref_snp_alleles")
        new_ref_snp_allele_id = row[0] or 0
        while True:
            # TODO: Multi-processing w/ file offsets, and Lock on
            # 'ref_snp_allele_id'
            # TODO: Async read from file, and write to DB.
            lines = ",\n".join(self.ref_variant_generator.readlines(
                1024*1024*8))  # 8MB
            if not lines:
                return
            json_ls = json.loads(f"[{lines}]")
            (ref_snp_alleles, gene_ref_snp_alleles,
             snv_freqs, snv_clinical_diseases) = [], [], [], []
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
                    new_ref_snp_allele_id += 1
                    allele_annotation = self.get_allele_annotation(
                        rsnp_json, variant_allele['allele_idx'])
                    ref_snp_alleles.append((ref_snp_id,
                                            variant_allele['del_seq'],
                                            variant_allele['ins_seq'],
                                            variant_allele['position'],))
                    gene_ref_snp_alleles += [(new_ref_snp_allele_id,
                                              gene['locus'],
                                              gene['id'],)
                                             for gene in allele_annotation[
                                                 'genes']]
                    # NOTE: The 'observation' is stored in JSON redundantly...
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
                    ('ref_snp_id', 'del_seq', 'ins_seq', 'position'),
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
            # futures = []
            for j, (table_name, columns, records) in enumerate(insert_queries):
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
        conn.commit()
        conn.close()

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
