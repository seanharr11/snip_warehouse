

            is_homo = False
            is_untracked_variant = False
            if variant_allele['ref'] not in matched_snp_genotype:
                is_homo = True
                if closest_allele_variant['alt'] * 2 != matched_genotype:
                    is_untracked_variant = True
                    self.unexpected_snp_cnt += 1
                else:
                    self.homozygous_snp_cnt += 1
            else:
                if closest_allele_variant['ref']*2 == matched_genotype:
                    # We do not have the less frequent allele!
                    continue
                else:
                    self.heterozygous_snp_cnt += 1
            variant_type = rsnp_json['primary_snapshot_data'][
                                    'variant_type']
            closest_allele_variant.update({
                'variant_type': variant_type,
                'refsnp_id': rsnp_json['refsnp_id'],
                'is_homo': is_homo,
                'is_untracked_variant': is_untracked_variant,
                'position': my_pos
            })
            closest_allele_variant.update(allele_annotation)
