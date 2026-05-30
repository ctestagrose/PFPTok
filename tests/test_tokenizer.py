import pytest
from pfptok import TokenizerManager, prefix_free_parse


SEQ = "ACGTACGTACGTACGTACGT"
SEQUENCES = [
    ["ACGTACGTACGT" * 5, "TGCATGCATGCA" * 5],
    ["AAACCCTTTGGG" * 5, "GGGAAACCCTT" * 5],
]


def test_prefix_free_parse_returns_list():
    phrases = prefix_free_parse(SEQ, w=4, p=7)
    assert isinstance(phrases, list)
    assert len(phrases) > 0


def test_prefix_free_parse_reconstructs_sequence():
    phrases = prefix_free_parse(SEQ, w=4, p=7)
    assert "".join(phrases) == SEQ


def test_prefix_free_parse_md5():
    phrases = prefix_free_parse(SEQ, w=4, p=7, use_simple_hash=False)
    assert "".join(phrases) == SEQ


def test_prefix_free_parse_short_sequence():
    # sequence must be longer than w to produce phrases
    short = "ACGTACGT"
    phrases = prefix_free_parse(short, w=4, p=7)
    assert "".join(phrases) == short


def test_setup_tokenizer_returns_tokenizer():
    tm = TokenizerManager()
    tok = tm.setup_tokenizer(SEQUENCES, w=6, p=117)
    assert tok is not None
    assert tok.get_vocab_size() > 0


def test_vocab_contains_special_tokens():
    tm = TokenizerManager()
    tok = tm.setup_tokenizer(SEQUENCES, w=6, p=117)
    vocab = tok.get_vocab()
    for special in ["[UNK]", "[CLS]", "[SEP]", "[PAD]", "[MASK]"]:
        assert special in vocab


def test_encode_sequences():
    tm = TokenizerManager()
    tm.setup_tokenizer(SEQUENCES, w=6, p=117)
    seqs = [["ACGTACGTACGT"], ["TGCATGCATGCA"]]
    ids, (unk_count, non_unk_count) = tm.encode_sequences(seqs)
    assert isinstance(ids, list)
    assert len(ids) > 0
    assert unk_count + non_unk_count == len(ids)


def test_encode_unk_rate_low_for_seen_sequences():
    tm = TokenizerManager()
    tm.setup_tokenizer(SEQUENCES, w=6, p=117)
    seqs = [[s] for group in SEQUENCES for s in group]
    _, (unk_count, non_unk_count) = tm.encode_sequences(seqs)
    total = unk_count + non_unk_count
    assert total > 0
    assert unk_count / total == 0.0


def test_save_and_load_tokenizer(tmp_path):
    tm = TokenizerManager()
    tok = tm.setup_tokenizer(SEQUENCES, w=6, p=117)
    path = str(tmp_path / "tok.json")
    tm.save_tokenizer(tok, path)

    tm2 = TokenizerManager()
    loaded = tm2.load_tokenizer(path)
    assert loaded.get_vocab_size() == tok.get_vocab_size()
    assert tm2.w == 6
    assert tm2.p == 117


def test_load_restores_correct_encoding(tmp_path):
    tm = TokenizerManager()
    tok = tm.setup_tokenizer(SEQUENCES, w=6, p=117)
    path = str(tmp_path / "tok.json")
    tm.save_tokenizer(tok, path)

    tm2 = TokenizerManager()
    tm2.load_tokenizer(path)
    seqs = [["ACGTACGTACGT"]]
    ids1, _ = tm.encode_sequences(seqs)
    ids2, _ = tm2.encode_sequences(seqs)
    assert ids1 == ids2
