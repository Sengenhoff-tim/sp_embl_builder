    use crate::types::UniProtIsoId;


    #[derive(Debug, Clone)]
    pub(crate) struct Variant {
        pub(crate) isoform: Option<UniProtIsoId>,
        pub(crate) id: String,
        pub(crate) begin: usize,
        pub(crate) end: usize,
        pub(crate) aa_ref: Option<String>, 
        pub(crate) aa_new: Option<String>, // None encodes missing sequence
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
                self.id.push_str(&other.id);
            }
        }
    }

    fn is_valid_aa(s: &str) -> bool {
        s.chars().all(|c| matches!(c, 'G' | 'A' | 'S' | 'P' | 'V' | 'T' | 'C' | 'L' | 'I' | 'N' | 'D' | 'Q' | 'K' | 'E' | 'M' | 'H' | 'F' | 'U' | 'R' | 'Y' | 'W' | 'O' | 'J' | 'X' | 'Z' | 'B'))
    }


    impl Variant {
        pub(crate) fn is_valid(&self) -> bool{
            if let Some(aa_ref) = &self.aa_ref {
                if !is_valid_aa(&aa_ref) {
                    return false;
                }
            }
            if let Some(aa_new) = &self.aa_new {
                if !is_valid_aa(&aa_new) {
                    return false;
                }
            }
            true
        }
    }

impl Variant {
    pub(crate) fn resolve_stop(&mut self, seq: &str) -> anyhow::Result<()> {
        if seq.is_empty() {
            return Err(anyhow::anyhow!(
                "Isoform sequence is empty. Variant at position [{}] is omitted",
                self.begin
            ));
        }

        if let Some(aa_ref) = self.aa_ref.clone() && let Some(aa_new) = self.aa_new.clone() {
            // "del" sequence: alternative missing encoding in EBI. Positions are assumed to stay correct
            if aa_new == "del" {
                self.aa_new = None;
                return Ok(());
            }

            let ref_star_idx = aa_ref.find('*');
            let new_star_idx = aa_new.find('*');

            if ref_star_idx.is_none() && new_star_idx.is_none() {
                return Ok(());
            }

            tracing::warn!(
                "Variant {} contains a stop codon marker [*] (aa_ref='{}', aa_new='{}'); \
                 best-effort truncating at first [*] to fit target encoding",
                self.id, aa_ref, aa_new
            );

            // aa_new, truncated at its first '*' if present (position-independent).
            let new_trunc = match new_star_idx {
                Some(idx) => aa_new[..idx].to_string(),
                None => aa_new.clone(),
            };

            if ref_star_idx == Some(0) {
                // aa_ref is only "*": no real reference residue exists at this
                // position (downstream can't append), so synthesize an anchor
                // from the isoform's last residue and fold the insertion into
                // aa_new: e.g. seq="AAA", ref="*", new="BB*" -> ref="A", new="ABB",
                // begin=end=len(seq), i.e. "replace the last A with ABB".
                let last = seq.chars().last().unwrap();
                self.aa_ref = Some(last.to_string());
                self.aa_new = Some(format!("{}{}", last, new_trunc));
                self.begin = seq.len();
                self.end = seq.len();
            } else {
                self.aa_new = Some(new_trunc);

                if let Some(idx) = ref_star_idx {
                    // Residues precede the '*': stop is lost further downstream,
                    // so the replaced span truly runs to the end of the isoform.
                    self.aa_ref = Some(aa_ref[..idx].to_string());
                    self.end = seq.len();
                } else if new_star_idx.is_some() {
                    // No '*' in aa_ref, but aa_new gains a stop: everything
                    // downstream of this position in the isoform is truncated.
                    self.end = seq.len();
                }
            }
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
