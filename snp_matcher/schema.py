# coding: utf-8
from sqlalchemy import (
    BigInteger, Column, ForeignKey,
    Integer, Text, create_engine)
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class Gene(Base):
    __tablename__ = 'genes'

    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, index=True)
    tax_id = Column(Integer)
    locus = Column(Text, index=True)
    symbol = Column(Text)
    typ = Column(Text)
    description = Column(Text)


class GeneSnv(Base):
    __tablename__ = "gene_snvs"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, index=True)
    snv_id = Column(ForeignKey('genes.id'), index=True)
    locus = Column(Text, index=True)


class Snv(Base):
    __tablename__ = 'snvs'

    id = Column(Integer, primary_key=True)
    rsnp_id = Column(Text, index=True)
    ref_seq = Column(Text, index=True)
    alt_seq = Column(Text, index=True)
    position = Column(BigInteger)


class SnvFrequency(Base):
    __tablename__ = 'snv_frequencies'

    id = Column(Integer, primary_key=True)
    project_name = Column(Text)
    allele_count = Column(Integer)
    total_count = Column(Integer)
    snv_id = Column(ForeignKey('snvs.id'), index=True)


class SnvClinicalDiseaseName(Base):
    __tablename__ = 'snv_clinical_disease_names'

    id = Column(Integer, primary_key=True)
    snv_id = Column(ForeignKey('snvs.id'), index=True)
    disease_name_csv = Column(Text)
    clinical_significance_csv = Column(Text)
    citation_csv = Column(Text)


if __name__ == "__main__":
    engine = create_engine("sqlite:///../data/snps.sql")
    metadata.create_all(engine)
