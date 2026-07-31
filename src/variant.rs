use crate::types::UniProtIsoId;

/// Represents a single amino acid change: residues `begin..=end` (`replaced`) are
/// swapped for `replacement`. An empty `replacement` means the range is deleted
/// ("Missing"). Only the affected range is stored, not the full before/after
/// sequence. `begin`/`end` are in canonical coordinates, unless a variant is
/// isoforms specific. Isoform specific variants are written as "ISO_ID:POS" in the FT line,
/// e.g., "P04637-2:45" rather than a bare "45".
/// Equality and Hash deliberately ignores `id`, because sample specific varaints may 
/// match a known variant from UniProt or EBI.
/// `end = begin + len(replaced) - 1` (single-residue has begin == end).
#[derive(Debug, Clone)]
pub(crate) struct Variant {
    pub(crate) isoform: Option<UniProtIsoId>,
    pub(crate) id: String,
    pub(crate) begin: usize,
    pub(crate) end: usize,
    pub(crate) aa_ref: String,
    pub(crate) aa_new: String,
}

//skips id field
impl PartialEq for Variant {
    fn eq(&self, other: &Self) -> bool {
        self.isoform == other.isoform
            && self.begin == other.begin
            && self.end == other.end
            && self.aa_ref == other.aa_ref
            && self.aa_new == other.aa_new
    }
}

impl Eq for Variant {}

//skips id field
impl std::hash::Hash for Variant {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.isoform.hash(state);
        self.begin.hash(state);
        self.end.hash(state);
        self.aa_ref.hash(state);
        self.aa_new.hash(state);
    }
}

impl Variant {
    /// Merge two identical variants by joining their IDs with `|`.
    pub(crate) fn merge(&mut self, other: &Variant) {
        if self.id.is_empty() {
            self.id = other.id.clone();
        } else {
            self.id.push('|');
            self.id.push_str(&other.id);
        }
    }
}

impl Variant {
    pub(crate) fn resolve_stop(&mut self, seq: &str) -> anyhow::Result<()> {
        if self.aa_ref == "*" {
            let last = seq
                .chars()
                .last()
                .ok_or_else(|| anyhow::anyhow!("cannot resolve stop: sequence is empty"))?
                .to_string();
            let mt = self.aa_new.trim_end_matches('*');
            self.begin = seq.len();
            self.end = seq.len();
            self.aa_ref = last.clone();
            self.aa_new = format!("{}{}", last, mt);
        } else if self.aa_new.ends_with('*') || self.aa_new == "del" {
            let suffix_start = self
                .begin
                .checked_add(self.aa_ref.len())
                .and_then(|n| n.checked_sub(1))
                .ok_or_else(|| anyhow::anyhow!("invalid begin/aa_ref length combination"))?;
            let suffix = seq.get(suffix_start..).ok_or_else(|| {
                anyhow::anyhow!(
                    "suffix_start {} out of bounds for sequence of length {}",
                    suffix_start,
                    seq.len()
                )
            })?;
            self.end = seq.len();
            self.aa_ref = format!("{}{}", self.aa_ref, suffix);
            self.aa_new = self.aa_new.trim_end_matches('*').to_string();
        }
        Ok(())
    }
}

pub(crate) trait JoinVariants {
    fn join(&mut self, variants_b: Vec<Variant>);
}

impl JoinVariants for Vec<Variant> {
    fn join(&mut self, variants_b: Vec<Variant>) {
        for variant_b in variants_b {
            if let Some(variant_a) = self.iter_mut().find(|v| *v == &variant_b) {
                variant_a.merge(&variant_b);
            } else {
                self.push(variant_b);
            }
        }
    }
}

