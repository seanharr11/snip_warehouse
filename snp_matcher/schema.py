# coding: utf-8
import asyncio
from gino import Gino

db = Gino()


class Gene(db.Model):
    __tablename__ = 'genes'

    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.Integer, index=True)
    tax_id = db.Column(db.Integer)
    locus = db.Column(db.Text, index=True)
    symbol = db.Column(db.Text)
    typ = db.Column(db.Text)
    description = db.Column(db.Text)


class GeneSnv(db.Model):
    __tablename__ = "gene_snvs"
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.Integer, index=True)
    snv_id = db.Column(db.ForeignKey('genes.id'), index=True)
    locus = db.Column(db.Text, index=True)


class Snv(db.Model):
    __tablename__ = 'snvs'

    id = db.Column(db.Integer, primary_key=True)
    rsnp_id = db.Column(db.Text, index=True)
    ref_seq = db.Column(db.Text, index=True)
    alt_seq = db.Column(db.Text, index=True)
    position = db.Column(db.BigInteger)


class SnvFrequency(db.Model):
    __tablename__ = 'snv_frequencies'

    id = db.Column(db.Integer, primary_key=True)
    project_name = db.Column(db.Text)
    allele_count = db.Column(db.Integer)
    total_count = db.Column(db.Integer)
    snv_id = db.Column(db.ForeignKey('snvs.id'), index=True)


class SnvClinicalDiseaseName(db.Model):
    __tablename__ = 'snv_clinical_disease_names'

    id = db.Column(db.Integer, primary_key=True)
    snv_id = db.Column(db.ForeignKey('snvs.id'), index=True)
    disease_name_csv = db.Column(db.Text)
    clinical_significance_csv = db.Column(db.Text)
    citation_csv = db.Column(db.Text)


class UserSnv(db.Model):
    __tablename__ = 'user_snvs'

    snv_id = db.Column(db.ForeignKey('snvs.id'), index=True)
    user_id = db.Column(db.ForeignKey('users.id'), index=True)
    genotype = db.Column(db.String(2), index=True)
    sort = db.Column(db.Integer, index=True)


class User(db.Model):
    __tablename__ = 'users'

    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String)
    fname = db.Column(db.String)
    lname = db.Column(db.String)
    sex = db.Column(db.String(1))
    ethnicity = db.Column(db.String)


async def main():
    await db.set_bind("postgresql://localhost/snvs")
    await db.gino.create_all()

if __name__ == "__main__":
    asyncio.get_event_loop().run_until_complete(main())
