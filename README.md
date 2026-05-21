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
   * The 10 example isolate files found in Sample_Data can be used to run the ablation script.
   * Update the two file paths:
     ```
         SEQUENCE_DIR="/path/to/train"
         TEST_SEQUENCE_DIR="/path/to/test"
     ```
   * Edit the configurations in the bash script (PFPTok below for example).
     ```
        # PFPTOK
        pfptok_quick)
            W_VALUES="5 10"
            D_VALUES="63 127 255"
            ;;
        pfptok_focused)
            W_VALUES="2000"
            D_VALUES="255 511 1021 4096"
            NUM_SEQUENCES=1000
            ;;
        pfptok_comprehensive)
            W_VALUES="20 100 1000 2000"
            D_VALUES="31 63 127 255 511 1021 4096"
            NUM_SEQUENCES=10
            ;;
     ```
   * Run the configuration for the tokenizer of choice
     ```
         bash run_ablation.sh <unigram|bpe|pfptok> <quick|focused|comprehensive>
     ```
   * Results will be printed to terminal and saved to the defined output location.
     
5. DNALONGBENCH Experiments

### Curated MTB Gene Classification Experiments
Code for these experiments can be found in the "Curated_Genes_Experiment" directory. 
Data can be downloaded from "https://github.com/ctestagrose/LLMTB/tree/main/Data"

Running:


### Whole Genome Tokenization/Classification Experiment

### Abalation Experiments

### DNALONGBENCH Experimements
