use super::*;

fn dummy_enst() -> EnsemblId {
    EnsemblId("ENST00000000000".to_string())
}

fn p(s: &str) -> Variant {
    parse_variant(s).expect("parse_variant failed")
}

fn r(s: &str, seq: &str) -> Variant {
    let entries = vec![(dummy_enst(), 0, p(s))];
    resolve_stop_variants(Some(&entries), seq)
        .unwrap()
        .into_iter()
        .next()
        .unwrap()
        .2
}

// ── non-stop variants (no sequence needed) ────────────────────────────────

#[test]
fn simple_substitution() {
    // "1A>1B" → begin=1, end=1, replaced="A", replacement="B"
    let v = p("1A>1B");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (1, 1, "A", "B"));
}

#[test]
fn one_to_multi() {
    // "1A>1BB" → begin=1, end=1, replaced="A", replacement="BB"
    let v = p("1A>1BB");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (1, 1, "A", "BB"));
}

#[test]
fn multi_to_one() {
    // "2AA>2B" → begin=2, end=3, replaced="AA", replacement="B"
    let v = p("2AA>2B");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (2, 3, "AA", "B"));
}

#[test]
fn multi_to_multi() {
    // "2AA>2BB" → begin=2, end=3, replaced="AA", replacement="BB"
    let v = p("2AA>2BB");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (2, 3, "AA", "BB"));
}

// ── stop variants (sequence = "ABCDE", len=5) ─────────────────────────────

#[test]
fn truncation_single_wt() {
    // "1A>1*" on "ABCDE": stop at position 1 deletes the whole protein →
    // begin=1, end=5, replaced="ABCDE" (the entire deleted range), replacement="" (Missing)
    let v = r("1A>1*", "ABCDE");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (1, 5, "ABCDE", ""));
}

#[test]
fn stop_extension_single_wt() {
    // "1A>1BB*" on "ABCDE" → begin=1, end=5, replaced="ABCDE", replacement="BB"
    let v = r("1A>1BB*", "ABCDE");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (1, 5, "ABCDE", "BB"));
}

#[test]
fn c_terminal_extension() {
    // "1*>1BB*" on "ABCDE" → begin=5, end=5, replaced="E", replacement="EBB"
    let v = r("1*>1BB*", "ABCDE");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (5, 5, "E", "EBB"));
}

// ── multi-char wt stop variants ───────────────────

#[test]
fn truncation_multi_wt() {
    // "1AB>1*" on "ABCDE" → begin=1, end=5, replaced="ABCDE", replacement=""
    // suffix starts after the raw "AB" (len=2): seq[2..]="CDE" appended to "AB"
    let v = r("1AB>1*", "ABCDE");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (1, 5, "ABCDE", ""));
}

#[test]
fn stop_extension_multi_wt() {
    // "1AB>1XY*" on "ABCDE" → begin=1, end=5, replaced="ABCDE", replacement="XY"
    // suffix starts after the raw "AB" (len=2): seq[2..]="CDE" appended to "AB"
    let v = r("1AB>1XY*", "ABCDE");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (1, 5, "ABCDE", "XY"));
}

// ── mid-sequence stop ─────────────────────────────────────────────────────

#[test]
fn truncation_mid_sequence() {
    // "3C>3*" on "ABCDE": stop at position 3 truncates the protein to "AB",
    // i.e. deletes the tail "CDE" (positions 3-5) → begin=3, end=5,
    // replaced="CDE" (the deleted range, not the whole sequence), replacement=""
    let v = r("3C>3*", "ABCDE");
    assert_eq!((v.begin, v.end, v.replaced.as_str(), v.replacement.as_str()), (3, 5, "CDE", ""));
}
