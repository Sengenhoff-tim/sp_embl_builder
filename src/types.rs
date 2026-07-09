// ============================================================================
// Core domain types shared across modules
// ============================================================================

/// Represents a single amino acid change: residues `begin..=end` (`replaced`) are
/// swapped for `replacement`. An empty `replacement` means the range is deleted
/// ("Missing"). Only the affected range is stored, not the full before/after
/// sequence. `begin`/`end` are in canonical coordinates, unless `isoform_ref` is
/// set — then they're in that isoform's own coordinates (see below).
/// `isoform_ref`, when set, is the isoform accession (e.g. "P04637-2") this
/// variant is specific to. It's written into the FT position field as
/// "P04637-2:45" rather than a bare "45" — ProtGraph's biopython-based FT parser
/// reads the "accession:position" form as a remote cross-reference, setting
/// `location.ref`, which ProtGraph then uses to route the feature onto that
/// isoform's own vertex chain (tagged `isoform_accession`/`isoform_position` by
/// its VAR_SEQ reconstruction) instead of the canonical chain. A bare canonical
/// position instead lands on every chain — canonical's and any isoform's — that
/// still carries that same position, which is how canonical-only variants reach
/// isoforms without any remapping on our end.
/// Equality and Hash deliberately ignore `id` and `isoform_ref` — both are
/// presentation-only, not part of the variant's identity.
/// `end = begin + len(replaced) - 1` (single-residue has begin == end).
#[derive(Debug, Clone)]
pub(crate) struct Variant {
    pub(crate) id: String,
    pub(crate) begin: usize,
    pub(crate) end: usize,
    pub(crate) replaced: String,
    pub(crate) replacement: String,
    pub(crate) isoform_ref: Option<UniprotId>,
}

impl PartialEq for Variant {
    fn eq(&self, other: &Self) -> bool {
        self.begin == other.begin
            && self.end == other.end
            && self.replaced == other.replaced
            && self.replacement == other.replacement
    }
}

impl Eq for Variant {}

impl std::hash::Hash for Variant {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.begin.hash(state);
        self.end.hash(state);
        self.replaced.hash(state);
        self.replacement.hash(state);
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct EnsemblId(pub(crate) String);
impl EnsemblId {
    pub(crate) fn as_str(&self) -> &str { &self.0 }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct UniprotId(pub(crate) String);
impl UniprotId {
    pub(crate) fn as_str(&self) -> &str { &self.0 }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct Sequence(pub(crate) String);
impl Sequence {
    pub(crate) fn as_str(&self) -> &str { &self.0 }
}

/// A single UniProt isoform: its full accession (e.g. P31946-2) and reconstructed sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Isoform(pub(crate) UniprotId, pub(crate) Sequence);
