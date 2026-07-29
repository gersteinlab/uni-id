use cosmwasm_std::StdError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ContractError {
    #[error("{0}")]
    Std(#[from] StdError),

    #[error("unauthorized: this operation is restricted to the data subject")]
    NotOwner {},

    #[error("unauthorized: sender is not on the whitelist")]
    NotWhitelisted {},

    #[error("the data key has already been generated and cannot be regenerated")]
    KeyAlreadyInitialized {},

    #[error("the data key has not been generated yet")]
    KeyNotInitialized {},

    #[error("no payload pointer has been set")]
    PayloadNotSet {},

    #[error("consensus randomness unavailable in this block")]
    NoRandomness {},

    #[error("uni_id must be a 64-character lowercase hex SHA-256 digest")]
    InvalidUniId {},

    #[error("content_hash must be a 64-character lowercase hex SHA-256 digest")]
    InvalidContentHash {},

    #[error("uri must be non-empty and at most 512 characters")]
    InvalidUri {},

    #[error("the data subject cannot be added to their own whitelist")]
    SelfGrant {},
}
