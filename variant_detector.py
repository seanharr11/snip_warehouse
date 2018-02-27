import json
import csv
import abc
from collections import namedtuple
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from snp_matcher.schema import (
    Gene, Allele, AlleleFrequency, AlleleClinicalDiseaseNames)
import sys


FrequencyStudy = namedtuple('FrequencyStudy', [
    'project_name',
    'allele_count',
    'total_count',
    'freq'])
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


class SNPVariantMatcher:
    def __init__(self):
        self._variant_output_filename = "./variants.csv"

        self.heterozygous_snp_cnt = 0
        self.homozygous_snp_cnt = 0
        self.unexpected_snp_cnt = 0
        self.rows_processed = 0

        self.matched_snps = {}
        self.variant_sources = []

        self.ref_variant_generator = open("refsnp-chr1.json", "r").readlines()

    def add_variant_source(self, variant_source):
        self.variant_sources.append(variant_source)

    def match(self):
        self.variant_source[0].load()
        self.matched_snps = {v.get_id(): v for v in self.variant_source[0]}
        for variant_source in self.variant_sources[1:]:
            variant_source.load()
            for variant in variant_source:
                if self.matched_snps.get(variant.get_id()):
                    self.matched_snaps[variant.get_id()].append(variant)

    def _print_status(self):
        if self.rows_processed % 10000 == 0:
            print(
                """Processed '{0}' Ref SNPs
                ---> Heterozygous Count: '{1}'
                ---> Homozygous Count: '{2}'
                ---> Untracked SNP Count: '{3}'
                """.format(str(self.rows_processed),
                           str(self.heterozygous_snp_cnt),
                           str(self.homozygous_snp_cnt),
                           str(self.unexpected_snp_cnt)))
        self.rows_processed += 1

    def detect_variants(self):
        field_names = ['refsnp_id', 'ref_pos', 'pos', 'distance', 'allele_cnt',
                       'total_cnt', 'allele_freq', 'variant_type',
                       'is_homozygous', 'is_untracked_variant', 'clinical',
                       'gene']
        with open(self.variant_output_filename, "w") as fp_w:
            variant_writer = csv.DictWriter(fp_w, fieldnames=field_names)
            variant_writer.writeheader()
            self.find_and_print_variants(variant_writer)

    def load_ref_snps(self):
        engine = create_engine("./data/snps.sql")
        for line in self.ref_variant_generator:
            self._check_and_log_status()
            rsnp_json = json.loads(line.decode('utf-8'))
            rsnp_placements = rsnp_json['primary_snapshot_data'][
                                    'placements_with_allele']
            refsnp_id = rsnp_json['refsnp_id']
            if not rsnp_placements:
                continue
            alleles = self.find_alleles_from_assembly(rsnp_placements)
            variant_allele = self.get_variant_allele(alleles)
            # TODO: Save this in a table
            session = sessionmaker()(bind=engine)
            allele_annotation = self.get_allele_annotation(
                                    rsnp_json,
                                    variant_allele['allele_idx'])
            genes = allele_annotation['genes']
            if len(genes) > 1:
                print(f"Multiple Genes: {genes}")
            gene = Gene(gene_id=genes[0]['id'],
                        locus=genes[0]['locus'],
                        name=genes[0]['name'])
            session.add(gene)
            allele = Allele(rsnp_id=refsnp_id,
                            ref_seq=variant_allele['ref_seq'],
                            alt_seq=variant_allele['alt_seq'],
                            position=variant_allele['position'],
                            gene_id=gene.id)
            session.add(allele)
            allele_freqs = [AlleleFrequency(
                project_name=af.project_name,
                allele_count=af.allele_count,
                total_count=af.total_count,
                allele_id=allele.id)
                    for af in allele_annotation['frequency_studies']]
            session.add(allele_freqs)
            allele_disease_name = AlleleClinicalDiseaseNames(
                allele_id=allele.id,
                disease_names=allele_annotation['disease_names'])
            session.add(allele_disease_name)
            session.commit()

    def find_alleles_from_assembly(rsnp_placements,
                                   assembly_name="GRCh38"):
        for rsnp_placement in rsnp_placements:
            annot = rsnp_placements.get('placement_annot')
            if not annot or not annot.get('seq_id_traits_by_assembly'):
                return
            assembly_info = annot['seq_id_traits_by_assembly']
            this_assembly_name = assembly_info.get("assembly_name") or None
            if this_assembly_name == assembly_name:
                alleles = rsnp_placement['alleles']
                return alleles

    def get_variant_allele(alleles):
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
        return {'frequency_studies': [FrequencyStudy(
                    project_name=freq['project_name'],
                    allele_count=freq['allele_count'],
                    total_count=freq['total_count'])
                        for freq in frequencies],
                'disease_names': "\n".join(
                    [clin['disease_names'] for clin in var_allele_annotation[
                        'clinical']]),
                # TODO: Multiple assembly annotations?
                'genes': assembly_annot[0]['genes']
                if assembly_annot else None,
                # TODO: Just store gene 'names' and 'locus'
                'seq_id': assembly_annot[0]['seq_id']}

    def load_snps(self):
        with open(self.snp_input_filename, "r") as fp_snps:
            dr_snps = csv.DictReader([l for l in fp_snps
                                      if not l.startswith('#')],
                                     delimiter='\t')
            for line in dr_snps:
                # rsid chromosome position genotype
                stripped_rsid = line['rsid'].replace("rs", "").replace("i", "")
                self.snps[stripped_rsid] = {
                    'genotype': line['genotype'],
                    'pos': int(line['position'])
                }

    if __name__ == "__main__":
        load_ref_snps(sys.argv[1])
        # my_snps = load_snps()
        # lookup_ref_snps(my_snps, ref_snp_filename)