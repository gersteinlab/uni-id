use cosmwasm_std::{Addr, Binary};
use schemars::JsonSchema;
use secret_toolkit::storage::{AppendStore, Item, Keymap};
use serde::{Deserialize, Serialize};

/// Contract configuration.
pub static CONFIG: Item<Config> = Item::new(b"config");

/// Symmetric key protecting the MP layer.
pub static DATA_KEY: Item<Binary> = Item::new(b"datakey");

/// Location of the encrypted MP layer.
pub static PAYLOAD: Item<Payload> = Item::new(b"payload");

/// Authorized recipients, keyed by address.
pub static WHITELIST: Keymap<Addr, Grant> = Keymap::new(b"whitelist");

/// Append-only record of key releases.
pub static ACCESS_LOG: AppendStore<AccessRecord> = AppendStore::new(b"accesslog");

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
pub struct Config {
    /// Sole administrator of this contract.
    pub owner: Addr,
    /// LP Uni-ID this contract governs, as a lowercase hex SHA-256 digest.
    pub uni_id: String,
    /// True once the data key has been generated.
    pub key_initialized: bool,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
pub struct Payload {
    /// Resolvable location of the ciphertext.
    pub uri: String,
    /// Lowercase hex SHA-256 digest of the ciphertext at `uri`.
    pub content_hash: String,
    /// Block height at which this pointer was last written.
    pub updated_at: u64,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
pub struct Grant {
    pub granted_at: u64,
    /// Optional label, readable only by the owner.
    pub note: Option<String>,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
pub struct AccessRecord {
    pub recipient: Addr,
    pub height: u64,
    pub time: u64,
}
