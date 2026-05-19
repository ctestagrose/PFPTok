# PFPTok
A Dictionary-Compression Approach to Genomic Tokenization via Prefix-Free Parsing

## Installation
```
cd PFPTok
pip install -r requirements.txt
```

## Usage
PFPTok expects a list of sequences (strings).

Tokenizer Training:
```
from PFPTok.src.PFP_Tokenizer import TokenizerManager

def main():
    sequences = [
        ["ACGT" * 25, "TGCA" * 25],
        ["AAAA" * 25, "CCCC" * 25],
    ]

    tm = TokenizerManager()
    tok = tm.setup_tokenizer(sequences, w=6, d=117)

if __name__ == "__main__":
    main()
```

## Running Experiments
1. Curated MTB Gene Classification Experiments
2. Whole Genome Tokenization/Classfication Experiments
3. Ablation Experiments
4. DNALONGBENCH Experiments

### Curated MTB Gene Classification Experiments
Code for these experiments can be found in the "Curated_Genes_Experiment" directory. 
Data can be downloaded from "https://github.com/ctestagrose/LLMTB/tree/main/Data"

Running:


### Whole Genome Tokenization/Classification Experiment

### Abalation Experiments

### DNALONGBENCH Experimements
