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
    locus = Column(Text, index=True)
    name = Column(Text)


class Allele(Base):
    __tablename__ = 'alleles'

    id = Column(Integer, primary_key=True)
    rnsp_id = Column(Text, index=True)
    ref_seq = Column(Text)
    alt_seq = Column(Text)
    position = Column(BigInteger)
    gene_id = Column(ForeignKey('genes.id'), index=True)


class AlleleFrequency(Base):
    __tablename__ = 'allele_frequencies'

    id = Column(Integer, primary_key=True)
    project_name = Column(Text)
    allele_count = Column(Integer)
    total_count = Column(Integer)
    allele_id = Column(ForeignKey('alleles.id'), index=True)


class AlleleClinicalDiseaseNames(Base):
    __tablename__ = 'allele_clinical_disease_names'

    id = Column(Integer, primary_key=True)
    allele_id = Column(ForeignKey('alleles.id'), index=True)
    disease_names = Column(Text)


engine = create_engine("sqlite:///data/snps.sql")
metadata.create_all(engine)
