# PFPTok

A Dictionary-Compression Approach to Genomic Tokenization via Prefix-Free Parsing

## Overview

PFPTok is a genomic tokenizer that uses prefix-free parsing (PFP) to segment DNA sequences into variable-length phrases. Instead of learning a fixed vocabulary from subword statistics (as BPE and SentencePiece do), PFPTok uses a rolling hash to identify phrase boundaries that are deterministic and reproducible for any given `(w, d)` parameter pair. The resulting phrases form a dictionary that directly serves as the tokenizer vocabulary.

### How it works

1. A Karp-Rabin rolling hash slides a window of size `w` across the input sequence.
2. Positions where `hash(window) mod d == 0` are marked as trigger points.
3. The sequence is split at these triggers into overlapping phrases, then adjusted to produce non-overlapping tokens.
4. All observed phrases are collected into a vocabulary.
5. The vocabulary is wrapped into a HuggingFace `tokenizers.Tokenizer` (WordLevel model) for fast encoding.

The window size `w` and divisor `d` together control phrase granularity: smaller `d` values produce more frequent triggers (shorter phrases, larger vocabulary), while larger `d` values yield longer phrases and a more compact vocabulary.

## Installation

```bash
git clone https://github.com/ctestagrose/PFPTok.git
cd PFPTok
pip install tokenizers tqdm
```

> **Note:** The included `requirements.txt` reflects a full research environment. The core dependencies for PFPTok are just `tokenizers` and `tqdm`.

## Usage

### Training a tokenizer

`setup_tokenizer` expects a list of sequence groups, where each group is a list of DNA strings (e.g., genes per isolate, contigs per assembly).

```python
from pfptok import TokenizerManager

sequences = [
    ["ACGTACGTACGT" * 10, "TGCATGCATGCA" * 10],
    ["AAACCCTTTGGG" * 10, "GGGAAACCCTT" * 10],
]

tm = TokenizerManager()
tokenizer = tm.setup_tokenizer(sequences, w=6, d=117)

print(f"Vocabulary size: {tokenizer.get_vocab_size()}")
```

### Encoding sequences

After training, use `encode_sequences` to convert sequences to token ID lists. Each entry is a `[sequence]` singleton list.

```python
seqs_to_encode = [
    ["ACGTACGTACGTACGT"],
    ["TGCATGCATGCATGCA"],
]

token_ids, (unk_count, non_unk_count) = tm.encode_sequences(seqs_to_encode)
print(f"Token IDs: {token_ids}")
print(f"UNK rate: {unk_count / (unk_count + non_unk_count):.4f}")
```

### Saving and loading

```python
tm.save_tokenizer(tokenizer, "my_tokenizer.json")

tm2 = TokenizerManager()
loaded = tm2.load_tokenizer("my_tokenizer.json")
```

### Using the PFP parser directly

The `prefix_free_parse` function can be used standalone to inspect how a sequence gets segmented:

```python
from src.PFP_Tokenizer import prefix_free_parse

phrases = prefix_free_parse("ACGTACGTACGTACGTACGT", w=4, d=7)
print(phrases)  # List of non-overlapping phrase strings
```

## API Reference

### `TokenizerManager`

| Method | Description |
|--------|-------------|
| `setup_tokenizer(sequences, w, d)` | Train a PFP tokenizer on the provided sequences. Returns a `tokenizers.Tokenizer`. |
| `encode_sequences(sequences, seed=None)` | Encode sequences into token ID lists. Returns `(ids, (unk_count, non_unk_count))`. |
| `save_tokenizer(tokenizer, path)` | Save a trained tokenizer to a JSON file. |
| `load_tokenizer(path)` | Load a tokenizer from a JSON file. Returns a `tokenizers.Tokenizer`. |

### `prefix_free_parse(sequence, w, d, use_simple_hash=True)`

Parse a single sequence into non-overlapping phrases.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequence` | `str` | — | Input DNA sequence |
| `w` | `int` | `10` | Window size for the rolling hash |
| `d` | `int` | `127` | Hash divisor — triggers occur where `hash mod d == 0` |
| `use_simple_hash` | `bool` | `True` | Use Karp-Rabin rolling hash (`True`) or MD5 (`False`) |

### Parameters

| Parameter | Controls | Guidance                                                                                                                               |
|-----------|----------|----------------------------------------------------------------------------------------------------------------------------------------|
| `w` | Rolling hash window size | Larger values produce fewer, more context-dependent trigger points                                                                     |
| `d` | Hash modulus / trigger frequency | Larger values -> fewer triggers -> longer phrases -> smaller vocab. Smaller values -> more triggers -> shorter phrases -> larger vocab |

## Experiments

Experiment code for evaluating PFPTok against BPE and Unigram across MTB antibiotic resistance classification, DNALongBench benchmarks, and hyperparameter ablation studies is available in the [PFPTok-Experiments](https://github.com/ctestagrose/PFPTok-Experiments) repository.
