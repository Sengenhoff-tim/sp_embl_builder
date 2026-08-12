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
                self.id.push_str("_");
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

        match (self.aa_ref.clone(), self.aa_new.clone()) {
            // Neither field present: nothing to resolve.
            (None, None) => Ok(()),


            
            // Only aa_new present, aa_ref missing: same — nothing to resolve.
            (None, Some(_)) => Err(anyhow::anyhow!(
                "Invalid variant {}: ref seq must be set",
                self.id
            )),

            // Only aa_ref present, aa_new missing: -> Map to missing
            (Some(aa_ref), None) => self.resolve_star(seq, &aa_ref, ""),
            
            // Both present: the only shape resolve_stop actually handles.
            (Some(aa_ref), Some(aa_new)) => self.resolve_star(seq, &aa_ref, &aa_new),
        }
    }

    fn resolve_star(&mut self, seq: &str, aa_ref: &str, aa_new: &str) -> anyhow::Result<()> {
        // "del" sequence: alternative missing encoding in EBI. Positions are assumed to stay correct
        if aa_new == "del" {
            self.aa_new = None;
            return Ok(());
        }

        let ref_star = aa_ref.find('*');
        let new_star = aa_new.find('*');

        if ref_star.is_none() && new_star.is_none() {
            return Ok(());
        }

        // Expected/common shape: both fields end in a single trailing '*'.
        // Anything else (mid-string '*', lone '*' ref, asymmetric '*') is
        // unusual enough to call out explicitly.
        let is_common_case = aa_ref.ends_with('*') && aa_new.ends_with('*');
        if !is_common_case {
            tracing::info!(
                "Variant {} contains a stop codon marker [*] (aa_ref='{}', aa_new='{}'); \
                 best-effort truncating at first [*] to fit target encoding",
                self.id, aa_ref, aa_new
            );
        }

        self.apply_stop_truncation(seq, aa_ref, aa_new, ref_star, new_star);

        // After resolution, if the replacement collapsed to nothing, encode
        // that as "missing sequence" on both fields
        if self.aa_new.as_deref() == Some("") {
            self.aa_ref = None;
            self.aa_new = None;
        }

        Ok(())
    }

    /// Given the positions of the first '*' in `aa_ref`/`aa_new` (if any),
    /// truncate the fields and adjust begin/end to fit the target encoding.
    /// Precondition: at least one of `ref_star`, `new_star` is `Some`.
    fn apply_stop_truncation(
        &mut self,
        seq: &str,
        aa_ref: &str,
        aa_new: &str,
        ref_star: Option<usize>,
        new_star: Option<usize>,
    ) {
        let new_trunc = match new_star {
            Some(idx) => &aa_new[..idx],
            None => aa_new,
        };

        match (ref_star, new_star) {
            // aa_ref is only "*": no real reference residue exists at this
            // position (downstream can't append), so synthesize an anchor
            // from the isoform's last residue and fold the insertion into
            // aa_new: e.g. seq="AAA", ref="*", new="BB*" -> ref="A", new="ABB",
            // begin=end=len(seq), i.e. "replace the last A with ABB".
            (Some(0), _) => {
                let last = seq.chars().last().unwrap();
                self.aa_ref = Some(last.to_string());
                self.aa_new = Some(format!("{last}{new_trunc}"));
                self.begin = seq.len();
                self.end = seq.len();
            }

            // Residues precede the '*' in aa_ref (idx > 0 guaranteed by match
            // order): stop is lost/gained further downstream either way, so
            // the replaced span runs to the true end of the isoform.
            (Some(idx), _) => {
                self.aa_ref = Some(aa_ref[..idx].to_string());
                self.aa_new = Some(new_trunc.to_string());
                self.end = seq.len();
            }

            // No '*' in aa_ref, but aa_new gains a stop: everything downstream
            // of this position in the isoform is truncated.
            (None, Some(_)) => {
                self.aa_new = Some(new_trunc.to_string());
                self.end = seq.len();
            }

            (None, None) => unreachable!("caller guarantees at least one of ref_star/new_star is Some"),
        }
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
