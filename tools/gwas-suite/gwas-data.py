import numpy as np
import pandas as pd
import os

np.random.seed(42)
n_samples = 200
n_snps = 2000

# Sample IDs
samples = [f"S{i+1}" for i in range(n_samples)]

# Covariates
age = np.random.randint(20, 75, size=n_samples)
sex = np.random.choice([0,1], size=n_samples)  # 0 female, 1 male

# Genotypes
genotypes = np.zeros((n_samples, n_snps), dtype=np.int8)
maf_array = np.random.uniform(0.05, 0.5, size=n_snps)

for j in range(n_snps):
    p = maf_array[j]
    freqs = [(1-p)**2, 2*p*(1-p), p**2]
    genotypes[:, j] = np.random.choice([0,1,2], size=n_samples, p=freqs)

# True association signal at SNP index 1000
logit = -2 + 0.05*(age-50) + 0.4*sex + 0.6*genotypes[:,999]
prob = 1 / (1 + np.exp(-logit))
pheno = (np.random.rand(n_samples) < prob).astype(int)

# Save data
os.makedirs("gwas_data", exist_ok=True)

snp_ids = [f"rs{100000+j}" for j in range(n_snps)]
geno_df = pd.DataFrame(genotypes, index=samples, columns=snp_ids)
geno_df.to_csv("gwas_data/genotypes.tsv", sep="\t")

sample_df = pd.DataFrame({
    "sample_id": samples,
    "phenotype": pheno,
    "age": age,
    "sex": sex
})
sample_df.to_csv("gwas_data/phenotypes_covariates.tsv", sep="\t", index=False)

chroms = np.random.choice([str(i) for i in range(1,23)], size=n_snps)
positions = np.arange(100000, 100000 + 100*n_snps, 100)
snp_annot = pd.DataFrame({
    "snp_id": snp_ids,
    "chrom": chroms,
    "pos": positions,
    "maf": np.round(maf_array, 4)
})
snp_annot.to_csv("gwas_data/snp_annotation.tsv", sep="\t", index=False)

print("Generated data in ./gwas_data directory")
