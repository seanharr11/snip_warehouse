# import asyncio
import abc
import asyncpg
from collections import namedtuple
import ftplib
import ujson as json
import zlib

# from .schema import (
#     Snv, SnvFrequency, SnvClinicalDiseaseName, GeneSnv)

FrequencyStudy = namedtuple('FrequencyStudy', [
    'project_name',
    'allele_count',
    'total_count'])
# https://www.ncbi.nlm.nih.gov/books/NBK21088/table/ch5.ch5_t3/?report=objectonly


class Variant(metaclass=abc.ABCMeta):
    def __init__(self, variant_object):
        self.variant_object = variant_object
        # TODO: parse this into self attribuets!

    @abc.abstractmethod
    def get_id():
        pass


class SNPVariant(Variant):
    pass


# This takes a variant source, which contains a collection of various
# variants of different types, and standardizes the data
class DbSNPVariant(Variant):
    pass


class SnvLoader:
    def __init__(self, db_conn_string):
        self._variant_output_filename = "./variants.csv"
        self.rows_processed = 0
        self.db_conn_string = db_conn_string

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
            fp.write(decompressor.decompress(byte_chunk))
            transferred = transferred + 8192
            transferred_mb = transferred / 1024 / 1024
            if transferred_mb % 10 == 0:
                print(
                    f"Transferred {transferred_mb}MB / {size / 1024 / 1024}MB")

        print(f"Filesize: {size / 1024 / 1024 / 1024}GB")
        ftp.retrbinary(f"RETR {dbsnp_filename}", callback, blocksize=blocksize)

    def _print_status(self):
        if self.rows_processed % 10000 == 0:
            print(f"Processed '{self.rows_processed}' Ref SNPs")
        self.rows_processed += 1

    async def load_ref_snps(self, dbsnp_filename):
        self.ref_variant_generator = open(dbsnp_filename, "rb")
        snv_id = 0
        pool = await asyncpg.create_pool(user='SeanH', database="snvs")
        connections = [await pool.acquire() for _ in range(4)]
        conn = connections[0]
        for table_name in ["snvs", "gene_snvs", "genes",
                           "snv_frequencies", "snv_clinical_disease_names"]:
            await connections[0].execute(
                f"ALTER TABLE {table_name} DISABLE TRIGGER ALL;")
        while True:
            if snv_id == 50000:
                exit(0)
            lines = ",\n".join(self.ref_variant_generator.readlines(
                1024*1024*8))  # 8MB
            if not lines:
                return
            json_ls = json.loads(f"[{lines}]")
            snvs, gene_snvs, snv_freqs, snv_disease_names = [], [], [], []
            for rsnp_json in json_ls:
                self._print_status()
                rsnp_placements = rsnp_json['primary_snapshot_data'][
                                        'placements_with_allele']
                refsnp_id = rsnp_json['refsnp_id']
                if not rsnp_placements:
                    continue
                alleles = self.find_alleles_from_assembly(rsnp_placements)
                if not alleles:
                    continue
                variant_allele = self.get_variant_allele(alleles)
                # Save to Database
                snv_id += 1
                allele_annotation = self.get_allele_annotation(
                    rsnp_json, variant_allele['allele_idx'])
                snvs.append((refsnp_id,
                             variant_allele['ref_seq'],
                             variant_allele['alt_seq'],
                             variant_allele['position'],))
                gene_snvs += [(snv_id,
                               gene['locus'],
                               gene['id'],)
                              for gene in allele_annotation['genes']]
                snv_freqs += [(af.project_name,
                               af.allele_count,
                               af.total_count,
                               snv_id,)
                              for af in allele_annotation['frequency_studies']]
                snv_disease_names += [
                    (snv_id,
                     clin['disease_names'],
                     clin['clinical_significances'],
                     clin['citation_csv'],)
                    for clin in allele_annotation['clinical_entries']]

            insert_queries = [
                ('snvs', ('rsnp_id', 'ref_seq', 'alt_seq', 'position'), snvs),
                ('gene_snvs', ('snv_id', 'locus', 'gene_id'), gene_snvs),
                ('snv_frequencies',
                 ('project_name', 'allele_count', 'total_count'), snv_freqs),
                ('snv_clinical_disease_names',
                    ('snv_id', 'disease_name_csv', 'clinical_significance_csv',
                     'citation_csv'), snv_disease_names)]
            for j, (table_name, columns, records) in enumerate(insert_queries):
                conn.copy_records_to_table(table_name, columns=columns,
                                           records=records)

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
        allele_tup = [(a, i) for (i, a) in enumerate(alleles)
                      if a['allele']['spdi']['inserted_sequence'] !=
                      a['allele']['spdi']['deleted_sequence']][0]
        var_spdi = allele_tup[0]['allele']['spdi']
        allele_idx = allele_tup[1]
        return {
            'ref_seq': var_spdi['deleted_sequence'],
            'alt_seq': var_spdi['inserted_sequence'],
            'position': var_spdi['position'],
            'allele_idx': allele_idx
        }

    def get_allele_annotation(self, rsnp_obj, allele_idx):
        var_allele_annotation = rsnp_obj['primary_snapshot_data'][
                                    'allele_annotations'][
                                    allele_idx]
        assembly_annot = var_allele_annotation['assembly_annotation']
        frequencies = var_allele_annotation['frequency']
        fs = [FrequencyStudy(
            project_name=freq['project_name'],
            allele_count=freq['allele_count'],
            total_count=freq['total_count']) for freq in frequencies]
        clinical_entries = [{
            'citation_csv': ','.join(map(lambda c: str(c), clin['citations'])),
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
