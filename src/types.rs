// ============================================================================
// Core domain types shared across modules
// ============================================================================

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct EnsemblId(pub(crate) String);
impl EnsemblId {
    pub(crate) fn as_str(&self) -> &str { &self.0 }
}

impl std::str::FromStr for EnsemblId {
    type Err = std::convert::Infallible;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(s.to_string()))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) enum UniProtId {
    Id(UniProtCanonId),
    Iso(UniProtIsoId),
}

impl UniProtId {
    pub(crate) fn as_str(&self) -> &str {
        match self {
            Self::Id(id) => id.as_str(),
            Self::Iso(id) => id.as_str(),
        }
    }
}

impl std::str::FromStr for UniProtId {
    type Err = std::convert::Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.contains('-') {
            Ok(Self::Iso(UniProtIsoId(s.to_string())))
        } else {
            Ok(Self::Id(UniProtCanonId(s.to_string())))
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Deserialize)]
pub(crate) struct UniProtCanonId(pub(crate) String);

impl UniProtCanonId {
    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::str::FromStr for UniProtCanonId {
    type Err = std::convert::Infallible;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(s.to_string()))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Deserialize)]
pub(crate) struct UniProtIsoId(pub(crate) String);

impl UniProtIsoId {
    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::str::FromStr for UniProtIsoId {
    type Err = std::convert::Infallible;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(s.to_string()))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct Sequence(pub(crate) String);
impl Sequence {
    pub(crate) fn as_str(&self) -> &str { &self.0 }
}

impl std::str::FromStr for Sequence {
    type Err = std::convert::Infallible;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(s.to_string()))
    }
}

/// A single UniProt isoform: its full accession (e.g. P31946-2) and reconstructed sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Isoform(pub(crate) UniProtIsoId, pub(crate) Sequence);
