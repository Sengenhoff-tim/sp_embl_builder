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
            self.id.push('|');
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
        if seq.len() < 1 {
            return Err(anyhow::anyhow!("Isoform sequence is empty. Variant at position [{}] is ommited", self.begin));
        }
        if let Some(aa_ref) = self.aa_ref.clone() && let Some(aa_new) = self.aa_new.clone(){
            // "del" sequence: alternative missing encoding in EBI. Positions are assumed to stay correct
            if aa_new == "del" {
                self.aa_new = None;

                return Ok(());
            }

            // resolve stop lost: set aa_ref to last(seq), aa_new' to last(seq) + aa_new
            if aa_ref.ends_with("*") {
                if aa_ref.len() <= 1 {
                    return Err(anyhow::anyhow!("Malformed variant: reference sequence ends with [*] and is longer than one character"));
                }
                self.end = seq.len();
                // safe beacause seq.len() > 1
                let last = seq.chars().last().unwrap();

                self.aa_ref = Some(last.to_string());
                self.aa_new = Some(format!("{}{}", last, aa_new.trim_end_matches("*")));

            }

            // resolve stop gain: strip * and, set self.end to len(seq)
            if aa_new.ends_with('*') {
                self.end = seq.len();
                // Missing sequence: set ref and new to None
                if aa_new.len() <= 1 {
                    self.aa_ref = None;
                    self.aa_new = None;
                }
                else {
                    if let Some(suffix_end) = self
                        .begin
                        .checked_add(aa_ref.len())
                        .and_then(|n| n.checked_sub(1)) 
                        {
                        if let Some(seq_suffix) = seq.get(self.begin..suffix_end){
                            self.aa_ref = Some(format!("{}{}", aa_ref, seq_suffix));
                            self.aa_new = Some(aa_new.trim_end_matches("*").to_string())
                        } else {
                            return Err(
                                anyhow::anyhow!(
                                    "Variant {} claims residues {}..{} (aa_ref='{}', aa_new='{}').",
                                    self.id, self.begin, self.begin + aa_ref.len().saturating_sub(1),
                                    aa_ref, aa_new
                                )
                            )
                        }
                    } else {
                        return Err(
                            anyhow::anyhow!(
                                "Calculation for (begin={}, aa_ref='{}', len={}): \
                                computing begin + len(aa_ref) - 1 overflowed",
                                self.begin, aa_ref, aa_ref.len()
                            )
                        )
                    }   
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
