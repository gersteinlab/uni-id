use cosmwasm_std::{Addr, Binary};
use schemars::JsonSchema;
use secret_toolkit::permit::Permit;
use serde::{Deserialize, Serialize};

use crate::state::{AccessRecord, Payload};

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
pub struct InstantiateMsg {
    /// Lowercase hex SHA-256 LP Uni-ID of the data subject.
    pub uni_id: String,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ExecuteMsg {
    /// Generate the data key and return it to the owner. Owner only, once.
    InitKey {},

    /// Set the location and ciphertext digest of the MP layer. Owner only.
    SetPayload { uri: String, content_hash: String },

    /// Add a recipient to the whitelist. Owner only.
    Grant {
        recipient: String,
        note: Option<String>,
    },

    /// Remove a recipient from the whitelist. Owner only.
    Revoke { recipient: String },

    /// Return the data key and payload location. Whitelisted recipients only.
    RequestAccess {},
}

/// Execute responses, returned in the response `data` field.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ExecuteAnswer {
    InitKey { key: Binary },
    SetPayload { status: String },
    Grant { status: String },
    Revoke { status: String },
    RequestAccess { key: Binary, payload: Payload },
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum QueryMsg {
    Config {},
    WithPermit {
        permit: Permit,
        query: QueryWithPermit,
    },
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum QueryWithPermit {
    /// Owner or a whitelisted recipient.
    Payload {},
    /// Owner only.
    Whitelist {},
    /// Owner only.
    AccessLog { page: u32, page_size: u32 },
    /// Any signer, about their own address.
    IsAuthorized {},
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum QueryAnswer {
    Config {
        owner: Addr,
        uni_id: String,
        key_initialized: bool,
        payload_set: bool,
    },
    Payload {
        payload: Payload,
    },
    Whitelist {
        entries: Vec<WhitelistEntry>,
    },
    AccessLog {
        records: Vec<AccessRecord>,
        total: u32,
    },
    IsAuthorized {
        authorized: bool,
    },
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, JsonSchema)]
pub struct WhitelistEntry {
    pub recipient: Addr,
    pub granted_at: u64,
    pub note: Option<String>,
}
