#!/usr/bin/env python3
"""Flag transcription-factor genes from curated GO annotations.

A TF is defined as a gene whose product has *DNA-binding transcription factor
activity* (GO:0003700) -- i.e. binds DNA sequence-specifically to regulate
transcription. We take GO:0003700 and ALL its descendant terms (e.g. the
RNA-pol-II-specific children GO:0000981/0001227/0001228), then select every
C. elegans gene with an `enables` (aspect F) annotation to any of them.

Inputs : go-basic.obo (GO ontology), wb.gaf.gz (WormBase GO annotations)
Output : tf.tsv  (wbgene, is_tf)
"""
import gzip

ROOT = "GO:0003700"   # DNA-binding transcription factor activity

# --- parse GO is_a graph from the OBO ---
parents = {}            # term -> set(parent terms)
cur, obsolete = None, False
with open("go-basic.obo") as f:
    for line in f:
        line = line.rstrip("\n")
        if line == "[Term]":
            cur, obsolete = None, False
        elif line.startswith("id: GO:"):
            cur = line[4:]
            parents.setdefault(cur, set())
        elif line.startswith("is_obsolete: true"):
            obsolete = True
        elif line.startswith("is_a:") and cur and not obsolete:
            parents[cur].add(line.split()[1])

# children map, then BFS to collect all descendants of ROOT (+ ROOT itself)
children = {}
for c, ps in parents.items():
    for p in ps:
        children.setdefault(p, set()).add(c)
tf_terms, stack = {ROOT}, [ROOT]
while stack:
    for ch in children.get(stack.pop(), ()):
        if ch not in tf_terms:
            tf_terms.add(ch); stack.append(ch)
print(f"TF-activity GO terms (GO:0003700 + descendants): {len(tf_terms)}")

# --- select genes annotated (enables / aspect F, not NOT) to any TF term ---
tf_genes = set()
with gzip.open("wb.gaf.gz", "rt") as f:
    for line in f:
        if line.startswith("!"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        wbgene, qualifier, go_id, aspect = c[1], c[3], c[4], c[8]
        if aspect == "F" and "NOT" not in qualifier and go_id in tf_terms:
            tf_genes.add(wbgene)

with open("tf.tsv", "w") as f:
    f.write("wbgene\tis_tf\n")
    for wb in sorted(tf_genes):
        f.write(f"{wb}\t1\n")

print(f"C. elegans TF genes flagged: {len(tf_genes)}")
