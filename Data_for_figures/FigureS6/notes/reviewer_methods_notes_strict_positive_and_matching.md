# Reviewer benchmark methods notes

## Strict positive experimental evidence

For the main reviewer-facing analysis, strict positive evidence was defined as external evidence in one of the following harmonized evidence-strength categories:

1. `strong_immunogenicity_and_presentation`
2. `strong_immunogenicity`
3. `presentation`

This definition intentionally excludes negative assay evidence, prediction-only concordance, and binding-only evidence. It is therefore focused on evidence that an epitope or overlapping peptide has either been experimentally presented, experimentally immunogenic, or both.

## Broad supportive evidence

The broader supportive endpoint additionally includes:

1. `prediction_concordance`
2. `binding_evidence`

This broader endpoint is useful as a sensitivity analysis, but it is less stringent because it includes evidence that may support binding or database prediction concordance without direct functional immunogenicity or presentation evidence.

## Exact peptide matching

Exact peptide matching requires the EGG candidate mutant peptide sequence to be identical to the peptide sequence in the external evidence database.

## Relaxed min-8 peptide matching

Relaxed min-8 matching records a match when one peptide is fully contained within the other and the contained sequence is at least 8 amino acids long. This allows a shorter experimentally observed ligand or T-cell epitope to match a longer EGG candidate peptide, or vice versa. This is biologically useful because experimentally presented MHC-I peptides are often 8-11 amino acids, while candidate tables can include longer class-II or source-window peptide sequences.

## HLA matching

Peptide-only matching was used as the main enrichment endpoint because the goal was to ask whether EGG prioritizes candidates related to experimentally supported peptide sequences. Peptide+HLA matching was also retained as a stricter secondary check, but it is expected to be lower because many public database records contain missing, broad, low-resolution, non-human, or non-canonical HLA annotations.
