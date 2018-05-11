# coding: utf-8
import os
from sqlalchemy import (
     BigInteger, Integer, Text, String,
     ForeignKey, Column, create_engine)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


Base = declarative_base()
metadata = Base.metadata
engine = create_engine(os.environ["SNVS_DB_URL"])
metadata.bind = engine
smaker = sessionmaker(bind=engine)


class Gene(Base):
    __tablename__ = 'genes'

    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, index=True)
    tax_id = Column(Integer)
    locus = Column(Text, index=True)
    symbol = Column(Text)
    typ = Column(Text)
    description = Column(Text)


class GeneRefSnp(Base):
    __tablename__ = "gene_ref_snp_alleles"
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, index=True)
    ref_snp_allele_id = Column(ForeignKey('ref_snp_alleles.id'),
                               index=True)
    locus = Column(Text, index=True)


class RefSnpAllele(Base):
    __tablename__ = 'ref_snp_alleles'

    id = Column(Integer, primary_key=True)
    ref_snp_id = Column(Integer, index=True)
    ins_seq = Column(Text, index=True)
    del_seq = Column(Text, index=True)
    position = Column(BigInteger)
    chromosome = Column(String(2))


class RefSnpFrequencyStudy(Base):
    __tablename__ = 'ref_snp_allele_freq_studies'

    id = Column(Integer, primary_key=True)
    ref_snp_allele_id = Column(ForeignKey('ref_snp_alleles.id'),
                               index=True)
    name = Column(Text)
    allele_count = Column(Integer)
    total_count = Column(Integer)


class RefSnpClinicalDiseases(Base):
    __tablename__ = 'ref_snp_allele_clin_diseases'

    id = Column(Integer, primary_key=True)
    ref_snp_allele_id = Column(ForeignKey('ref_snp_alleles.id'),
                               index=True)
    disease_name_csv = Column(Text)
    clinical_significance_csv = Column(Text)
    citation_list = Column(ARRAY(Integer, dimensions=1))


class UserRefSnp(Base):
    __tablename__ = 'user_ref_snps'
    id = Column(Integer, primary_key=True)
    ref_snp_id = Column(Integer, index=True)
    user_id = Column(ForeignKey('users.id'), index=True)
    genotype = Column(String(2), index=True)
    sort = Column(Integer, index=True)


class User(Base):
    __tablename__ = 'users'

    id = Column(Integer, primary_key=True)
    email = Column(String)
    fname = Column(String)
    lname = Column(String)
    sex = Column(String(1))
    ethnicity = Column(String)


def init_db(database_name):
    os.system(f"echo 'DROP DATABASE {database_name};"
              f"CREATE DATABASE {database_name};'"
              " | psql -U SeanH postgres")
    metadata.create_all()
    session = smaker()
    session.add(
        User(fname="Sean", lname="Harrington", email="seanharr11@gmail.com"))

    session.commit()


if __name__ == "__main__":
    init_db()
