use cosmwasm_std::{
    entry_point, to_binary, Addr, Binary, Deps, DepsMut, Env, MessageInfo, Response, StdError,
    StdResult,
};
use secret_toolkit::permit::{validate, Permit};
use sha2::{Digest, Sha256};

use crate::error::ContractError;
use crate::msg::{
    ExecuteAnswer, ExecuteMsg, InstantiateMsg, QueryAnswer, QueryMsg, QueryWithPermit,
    WhitelistEntry,
};
use crate::state::{
    AccessRecord, Config, Grant, Payload, ACCESS_LOG, CONFIG, DATA_KEY, PAYLOAD, WHITELIST,
};

/// Storage prefix used by the permit module to track revoked permits.
const PREFIX_REVOKED_PERMITS: &str = "revoked_permits";

/// Domain separation tag for key derivation.
const KDF_DOMAIN: &[u8] = b"uniid-mp-access/aes256-data-key/v1";

const MAX_URI_LEN: usize = 512;

// ---------------------------------------------------------------------------
// Instantiate
// ---------------------------------------------------------------------------

#[entry_point]
pub fn instantiate(
    deps: DepsMut,
    _env: Env,
    info: MessageInfo,
    msg: InstantiateMsg,
) -> Result<Response, ContractError> {
    let uni_id = msg.uni_id.to_lowercase();
    if !is_hex_sha256(&uni_id) {
        return Err(ContractError::InvalidUniId {});
    }

    CONFIG.save(
        deps.storage,
        &Config {
            owner: info.sender.clone(),
            uni_id,
            key_initialized: false,
        },
    )?;

    Ok(Response::new().add_attribute("action", "instantiate"))
}

// ---------------------------------------------------------------------------
// Execute
// ---------------------------------------------------------------------------

#[entry_point]
pub fn execute(
    deps: DepsMut,
    env: Env,
    info: MessageInfo,
    msg: ExecuteMsg,
) -> Result<Response, ContractError> {
    match msg {
        ExecuteMsg::InitKey {} => try_init_key(deps, env, info),
        ExecuteMsg::SetPayload { uri, content_hash } => {
            try_set_payload(deps, env, info, uri, content_hash)
        }
        ExecuteMsg::Grant { recipient, note } => try_grant(deps, env, info, recipient, note),
        ExecuteMsg::Revoke { recipient } => try_revoke(deps, info, recipient),
        ExecuteMsg::RequestAccess {} => try_request_access(deps, env, info),
    }
}

/// Generate a 256-bit data key from consensus randomness and return it to the
/// owner in the encrypted response.
fn try_init_key(deps: DepsMut, env: Env, info: MessageInfo) -> Result<Response, ContractError> {
    let mut config = CONFIG.load(deps.storage)?;
    if info.sender != config.owner {
        return Err(ContractError::NotOwner {});
    }
    if config.key_initialized {
        return Err(ContractError::KeyAlreadyInitialized {});
    }

    let seed = env
        .block
        .random
        .ok_or(ContractError::NoRandomness {})?;

    let mut hasher = Sha256::new();
    hasher.update(KDF_DOMAIN);
    hasher.update(seed.as_slice());
    hasher.update(env.contract.address.as_bytes());
    hasher.update(config.uni_id.as_bytes());
    let key = Binary::from(hasher.finalize().as_slice());

    DATA_KEY.save(deps.storage, &key)?;
    config.key_initialized = true;
    CONFIG.save(deps.storage, &config)?;

    Ok(Response::new()
        .add_attribute("action", "init_key")
        .set_data(to_binary(&ExecuteAnswer::InitKey { key })?))
}

/// Record where the encrypted MP layer lives and what it hashes to.
fn try_set_payload(
    deps: DepsMut,
    env: Env,
    info: MessageInfo,
    uri: String,
    content_hash: String,
) -> Result<Response, ContractError> {
    let config = CONFIG.load(deps.storage)?;
    if info.sender != config.owner {
        return Err(ContractError::NotOwner {});
    }
    if uri.is_empty() || uri.len() > MAX_URI_LEN {
        return Err(ContractError::InvalidUri {});
    }
    let content_hash = content_hash.to_lowercase();
    if !is_hex_sha256(&content_hash) {
        return Err(ContractError::InvalidContentHash {});
    }

    PAYLOAD.save(
        deps.storage,
        &Payload {
            uri,
            content_hash,
            updated_at: env.block.height,
        },
    )?;

    Ok(Response::new()
        .add_attribute("action", "set_payload")
        .set_data(to_binary(&ExecuteAnswer::SetPayload {
            status: "ok".to_string(),
        })?))
}

/// Add a recipient to the whitelist.
fn try_grant(
    deps: DepsMut,
    env: Env,
    info: MessageInfo,
    recipient: String,
    note: Option<String>,
) -> Result<Response, ContractError> {
    let config = CONFIG.load(deps.storage)?;
    if info.sender != config.owner {
        return Err(ContractError::NotOwner {});
    }

    let recipient = deps.api.addr_validate(&recipient)?;
    if recipient == config.owner {
        return Err(ContractError::SelfGrant {});
    }

    WHITELIST.insert(
        deps.storage,
        &recipient,
        &Grant {
            granted_at: env.block.height,
            note,
        },
    )?;

    Ok(Response::new()
        .add_attribute("action", "grant")
        .set_data(to_binary(&ExecuteAnswer::Grant {
            status: "ok".to_string(),
        })?))
}

/// Remove a recipient from the whitelist. Blocks future releases only.
fn try_revoke(
    deps: DepsMut,
    info: MessageInfo,
    recipient: String,
) -> Result<Response, ContractError> {
    let config = CONFIG.load(deps.storage)?;
    if info.sender != config.owner {
        return Err(ContractError::NotOwner {});
    }

    let recipient = deps.api.addr_validate(&recipient)?;
    WHITELIST.remove(deps.storage, &recipient)?;

    Ok(Response::new()
        .add_attribute("action", "revoke")
        .set_data(to_binary(&ExecuteAnswer::Revoke {
            status: "ok".to_string(),
        })?))
}

/// Release the data key and payload location to a whitelisted recipient.
fn try_request_access(
    deps: DepsMut,
    env: Env,
    info: MessageInfo,
) -> Result<Response, ContractError> {
    let config = CONFIG.load(deps.storage)?;

    if WHITELIST.get(deps.storage, &info.sender).is_none() {
        return Err(ContractError::NotWhitelisted {});
    }
    if !config.key_initialized {
        return Err(ContractError::KeyNotInitialized {});
    }

    let key = DATA_KEY.load(deps.storage)?;
    let payload = PAYLOAD
        .may_load(deps.storage)?
        .ok_or(ContractError::PayloadNotSet {})?;

    ACCESS_LOG.push(
        deps.storage,
        &AccessRecord {
            recipient: info.sender.clone(),
            height: env.block.height,
            time: env.block.time.seconds(),
        },
    )?;

    Ok(Response::new()
        .add_attribute("action", "request_access")
        .set_data(to_binary(&ExecuteAnswer::RequestAccess { key, payload })?))
}

// ---------------------------------------------------------------------------
// Query
// ---------------------------------------------------------------------------

#[entry_point]
pub fn query(deps: Deps, env: Env, msg: QueryMsg) -> StdResult<Binary> {
    match msg {
        QueryMsg::Config {} => query_config(deps),
        QueryMsg::WithPermit { permit, query } => permit_query(deps, env, permit, query),
    }
}

/// Public metadata. Exposes no key material or whitelist entries.
fn query_config(deps: Deps) -> StdResult<Binary> {
    let config = CONFIG.load(deps.storage)?;
    to_binary(&QueryAnswer::Config {
        owner: config.owner,
        uni_id: config.uni_id,
        key_initialized: config.key_initialized,
        payload_set: PAYLOAD.may_load(deps.storage)?.is_some(),
    })
}

fn permit_query(
    deps: Deps,
    env: Env,
    permit: Permit,
    query: QueryWithPermit,
) -> StdResult<Binary> {
    let signer = validate(
        deps,
        PREFIX_REVOKED_PERMITS,
        &permit,
        env.contract.address.to_string(),
        None,
    )?;
    let signer = deps.api.addr_validate(&signer)?;
    let config = CONFIG.load(deps.storage)?;
    let is_owner = signer == config.owner;

    match query {
        QueryWithPermit::Payload {} => {
            if !is_owner && WHITELIST.get(deps.storage, &signer).is_none() {
                return Err(StdError::generic_err("not authorized to view the payload"));
            }
            let payload = PAYLOAD
                .may_load(deps.storage)?
                .ok_or_else(|| StdError::generic_err("no payload pointer has been set"))?;
            to_binary(&QueryAnswer::Payload { payload })
        }

        QueryWithPermit::Whitelist {} => {
            require_owner(is_owner)?;
            let mut entries: Vec<WhitelistEntry> = Vec::new();
            for item in WHITELIST.iter(deps.storage)? {
                let (recipient, grant): (Addr, Grant) = item?;
                entries.push(WhitelistEntry {
                    recipient,
                    granted_at: grant.granted_at,
                    note: grant.note,
                });
            }
            to_binary(&QueryAnswer::Whitelist { entries })
        }

        QueryWithPermit::AccessLog { page, page_size } => {
            require_owner(is_owner)?;
            let total = ACCESS_LOG.get_len(deps.storage)?;
            let records = ACCESS_LOG
                .paging(deps.storage, page, page_size)
                .unwrap_or_default();
            to_binary(&QueryAnswer::AccessLog { records, total })
        }

        QueryWithPermit::IsAuthorized {} => to_binary(&QueryAnswer::IsAuthorized {
            authorized: is_owner || WHITELIST.get(deps.storage, &signer).is_some(),
        }),
    }
}

fn require_owner(is_owner: bool) -> StdResult<()> {
    if is_owner {
        Ok(())
    } else {
        Err(StdError::generic_err(
            "unauthorized: this query is restricted to the data subject",
        ))
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// True for a 64-character lowercase hex string.
fn is_hex_sha256(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use cosmwasm_std::testing::{mock_dependencies, mock_env, mock_info};
    use cosmwasm_std::{from_binary, OwnedDeps};

    const UNI_ID: &str = "3f786850e387550fdab836ed7e6dc881de23001b0000000000000000000000ab";
    const CIPHER_HASH: &str = "aa61ce8b5e0b4b1e1c2f5d3e4a5b6c7d8e9f000112233445566778899aabbccd";
    const URI: &str = "s3://uniid-mp/NA12878.mp.enc";

    fn env_with_randomness() -> Env {
        let mut env = mock_env();
        env.block.random = Some(Binary::from(&[7u8; 32][..]));
        env
    }

    fn setup() -> OwnedDeps<
        cosmwasm_std::MemoryStorage,
        cosmwasm_std::testing::MockApi,
        cosmwasm_std::testing::MockQuerier,
    > {
        let mut deps = mock_dependencies();
        instantiate(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            InstantiateMsg {
                uni_id: UNI_ID.to_string(),
            },
        )
        .unwrap();
        deps
    }

    fn init_key(deps: &mut OwnedDeps<
        cosmwasm_std::MemoryStorage,
        cosmwasm_std::testing::MockApi,
        cosmwasm_std::testing::MockQuerier,
    >) -> Binary {
        let res = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::InitKey {},
        )
        .unwrap();
        match from_binary(&res.data.unwrap()).unwrap() {
            ExecuteAnswer::InitKey { key } => key,
            _ => panic!("wrong answer variant"),
        }
    }

    fn set_payload(deps: &mut OwnedDeps<
        cosmwasm_std::MemoryStorage,
        cosmwasm_std::testing::MockApi,
        cosmwasm_std::testing::MockQuerier,
    >) {
        execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::SetPayload {
                uri: URI.to_string(),
                content_hash: CIPHER_HASH.to_string(),
            },
        )
        .unwrap();
    }

    fn grant(
        deps: &mut OwnedDeps<
            cosmwasm_std::MemoryStorage,
            cosmwasm_std::testing::MockApi,
            cosmwasm_std::testing::MockQuerier,
        >,
        who: &str,
    ) {
        execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::Grant {
                recipient: who.to_string(),
                note: Some("IRB-2026-0142".to_string()),
            },
        )
        .unwrap();
    }

    #[test]
    fn instantiate_rejects_malformed_uni_id() {
        let mut deps = mock_dependencies();
        let err = instantiate(
            deps.as_mut(),
            mock_env(),
            mock_info("alice", &[]),
            InstantiateMsg {
                uni_id: "not-a-digest".to_string(),
            },
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::InvalidUniId {}));
    }

    #[test]
    fn key_is_256_bits_and_returned_to_owner() {
        let mut deps = setup();
        let key = init_key(&mut deps);
        assert_eq!(key.as_slice().len(), 32);
    }

    #[test]
    fn key_cannot_be_regenerated() {
        let mut deps = setup();
        init_key(&mut deps);
        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::InitKey {},
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::KeyAlreadyInitialized {}));
    }

    #[test]
    fn non_owner_cannot_init_key_or_grant() {
        let mut deps = setup();
        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("mallory", &[]),
            ExecuteMsg::InitKey {},
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::NotOwner {}));

        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("mallory", &[]),
            ExecuteMsg::Grant {
                recipient: "mallory".to_string(),
                note: None,
            },
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::NotOwner {}));
    }

    #[test]
    fn whitelisted_recipient_receives_the_same_key_the_owner_got() {
        let mut deps = setup();
        let owner_key = init_key(&mut deps);
        set_payload(&mut deps);
        grant(&mut deps, "carl");

        let res = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("carl", &[]),
            ExecuteMsg::RequestAccess {},
        )
        .unwrap();

        match from_binary(&res.data.unwrap()).unwrap() {
            ExecuteAnswer::RequestAccess { key, payload } => {
                assert_eq!(key, owner_key);
                assert_eq!(payload.uri, URI);
                assert_eq!(payload.content_hash, CIPHER_HASH);
            }
            _ => panic!("wrong answer variant"),
        }
    }

    #[test]
    fn non_whitelisted_request_is_rejected() {
        let mut deps = setup();
        init_key(&mut deps);
        set_payload(&mut deps);

        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("mallory", &[]),
            ExecuteMsg::RequestAccess {},
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::NotWhitelisted {}));
    }

    #[test]
    fn revoked_recipient_cannot_request_again() {
        let mut deps = setup();
        init_key(&mut deps);
        set_payload(&mut deps);
        grant(&mut deps, "carl");

        execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("carl", &[]),
            ExecuteMsg::RequestAccess {},
        )
        .unwrap();

        execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::Revoke {
                recipient: "carl".to_string(),
            },
        )
        .unwrap();

        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("carl", &[]),
            ExecuteMsg::RequestAccess {},
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::NotWhitelisted {}));
    }

    #[test]
    fn request_before_payload_is_set_fails() {
        let mut deps = setup();
        init_key(&mut deps);
        grant(&mut deps, "carl");

        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("carl", &[]),
            ExecuteMsg::RequestAccess {},
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::PayloadNotSet {}));
    }

    #[test]
    fn access_log_records_each_release() {
        let mut deps = setup();
        init_key(&mut deps);
        set_payload(&mut deps);
        grant(&mut deps, "carl");
        grant(&mut deps, "dana");

        for who in ["carl", "dana", "carl"] {
            execute(
                deps.as_mut(),
                env_with_randomness(),
                mock_info(who, &[]),
                ExecuteMsg::RequestAccess {},
            )
            .unwrap();
        }

        assert_eq!(ACCESS_LOG.get_len(&deps.storage).unwrap(), 3);
    }

    #[test]
    fn owner_cannot_whitelist_themselves() {
        let mut deps = setup();
        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::Grant {
                recipient: "alice".to_string(),
                note: None,
            },
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::SelfGrant {}));
    }

    #[test]
    fn set_payload_rejects_bad_hash_and_empty_uri() {
        let mut deps = setup();
        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::SetPayload {
                uri: URI.to_string(),
                content_hash: "deadbeef".to_string(),
            },
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::InvalidContentHash {}));

        let err = execute(
            deps.as_mut(),
            env_with_randomness(),
            mock_info("alice", &[]),
            ExecuteMsg::SetPayload {
                uri: String::new(),
                content_hash: CIPHER_HASH.to_string(),
            },
        )
        .unwrap_err();
        assert!(matches!(err, ContractError::InvalidUri {}));
    }

    #[test]
    fn public_config_query_exposes_no_private_state() {
        let mut deps = setup();
        init_key(&mut deps);
        grant(&mut deps, "carl");

        let res = query(deps.as_ref(), mock_env(), QueryMsg::Config {}).unwrap();
        match from_binary(&res).unwrap() {
            QueryAnswer::Config {
                owner,
                uni_id,
                key_initialized,
                payload_set,
            } => {
                assert_eq!(owner, Addr::unchecked("alice"));
                assert_eq!(uni_id, UNI_ID);
                assert!(key_initialized);
                assert!(!payload_set);
            }
            _ => panic!("wrong answer variant"),
        }
    }
}
