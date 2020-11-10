# shmoof

Companion github to "Learning the heterogeneous hypermutation landscape of immunoglobulins from high-throughput repertoire data".

Context and position mutabilities of the model learned on out-of-frame lineages from 9 individuals of Ref. 7 can be found in `mutabilities_context.tsv` and `mutabilities_position.tsv`.

## Example of use

After cloning the repository:
```python
import pandas as pd
df_context = pd.read_csv("shmoof/mutabilities_context.tsv", sep="\t").set_index("Motif")
df_pos = pd.read_csv("shmoof/mutabilities_position.tsv", sep="\t").set_index("Position")

#VH1-69
germline = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCGGTGAAGGTCTCCTGCAAGTCTTCTGGAGGCACCTTCAGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGACTGGATGGGAGGGATCATCCCTATCTTTGGTACAGCAAACTACGCACAGAAGTTCCAGGGCAGAGTCACGATTACCGCGGACAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGGCACGGGAATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA" 
# one mutation, at position 225
mutant = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCGGTGAAGGTCTCCTGCAAGTCTTCTGGAGGCACCTTCAGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGACTGGATGGGAGGGATCATCCCTATCTTTGGTACAGCAAACTACGCACAGAAGTTCCAGGGCAGAGTCACGATTACCGCGGACAAATCTACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGGCACGGGAATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA"

pos, mut = next((k+1, m) for k, (g, m) in enumerate(zip(germline, mutant)) if g != m)

mutation_rate = df_context.loc[germline[pos-2:pos+3]]["Mutability"] * df_pos.loc[pos+1]["Mutability"]
probability_substitution = df_context.loc[germline[pos-2:pos+3]]["To " + mut]

print(f"Mutation rate: {mutation_rate:.2}")
print(f"Probability of finding this specific mutation (given a mutation happened) {probability_substitution:.2}")
```
