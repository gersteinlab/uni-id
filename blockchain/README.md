# Blockchain access control for the MP layer

A Secret Network smart contract for controlling access to the most-private (MP) layer of a partitioned genomic file.

One contract instance governs one individual's MP layer. It generates and holds the symmetric key that encrypts that layer, maintains a private whitelist of authorized recipients, and stores a pointer to the encrypted data. Contract state is sealed inside a trusted execution environment, so the key and the whitelist are not readable by node operators. The encrypted data itself is stored off-chain.

## Prerequisites

- Rust 1.75 or later, with the WASM target:
  ```sh
  rustup target add wasm32-unknown-unknown
  ```
- [Docker](https://docs.docker.com/get-docker/), for producing an optimized build
- [`secretcli`](https://docs.scrt.network/secret-network-documentation/development/tools-and-libraries/secret-cli), configured against a node

## Install

```sh
git clone https://github.com/gersteinlab/uni-id.git
cd uni-id/blockchain
cargo build
```

## Test

```sh
cargo test
```

## Build for deployment

From `uni-id/blockchain`:

```sh
docker run --rm -v "$PWD":/contract \
  --mount type=volume,source=uniid_cache,target=/contract/target \
  --mount type=volume,source=registry_cache,target=/usr/local/cargo/registry \
  ghcr.io/scrtlabs/secret-contract-optimizer
```

This produces `contract.wasm.gz` in this directory.

## Deploy

Set your chain and key. Endpoints change, so check the [Secret Network docs](https://docs.scrt.network/secret-network-documentation/development/resources-api-contract-addresses/connecting-to-the-network) for the current testnet.

```sh
export CHAIN_ID=pulsar-3
export NODE=https://rpc.pulsar.scrttestnet.com
export KEY=alice
```

Upload the code:

```sh
secretcli tx compute store contract.wasm.gz \
  --from $KEY --gas 5000000 --chain-id $CHAIN_ID --node $NODE -y
```

Find the resulting code ID:

```sh
secretcli query compute list-code --node $NODE
```

Instantiate a contract for one individual, identified by their 64-character hex Uni-ID:

```sh
export CODE_ID=<code id from above>

secretcli tx compute instantiate $CODE_ID \
  '{"uni_id":"3f786850e387550fdab836ed7e6dc881de23001b0000000000000000000000ab"}' \
  --from $KEY --label uniid-mp-alice --gas 200000 \
  --chain-id $CHAIN_ID --node $NODE -y
```

Get the contract address:

```sh
secretcli query compute list-contract-by-code $CODE_ID --node $NODE
export CONTRACT=<contract address>
```

## Usage

The account that instantiates the contract is its owner and is the only account that can generate the key, set the payload pointer, or modify the whitelist.

Responses to execute messages are encrypted to the sender. After broadcasting a transaction, read the decrypted response with `secretcli query compute tx <TXHASH>`.

### 1. Generate the data key (owner)

```sh
secretcli tx compute execute $CONTRACT '{"init_key":{}}' \
  --from $KEY --gas 200000 --chain-id $CHAIN_ID --node $NODE -y
```

Read the returned key:

```sh
secretcli query compute tx <TXHASH> --node $NODE
```

The response contains a base64-encoded 256-bit key. Decode it and use it to encrypt the MP layer locally, for example with `openssl enc -aes-256-gcm`. This can only be run once per contract.

### 2. Encrypt and upload the MP layer

Done outside the contract. Encrypt the MP layer produced by the partitioning pipeline with the key from step 1, upload the ciphertext to wherever you host it, and compute its SHA-256:

```sh
sha256sum NA12878.mp.enc
```

### 3. Record where the data lives (owner)

```sh
secretcli tx compute execute $CONTRACT \
  '{"set_payload":{"uri":"s3://uniid-mp/NA12878.mp.enc","content_hash":"<sha256 hex>"}}' \
  --from $KEY --gas 200000 --chain-id $CHAIN_ID --node $NODE -y
```

### 4. Authorize a recipient (owner)

```sh
secretcli tx compute execute $CONTRACT \
  '{"grant":{"recipient":"secret1carl...","note":"IRB-2026-0142"}}' \
  --from $KEY --gas 200000 --chain-id $CHAIN_ID --node $NODE -y
```

`note` is optional and visible only to the owner.

To remove a recipient:

```sh
secretcli tx compute execute $CONTRACT \
  '{"revoke":{"recipient":"secret1carl..."}}' \
  --from $KEY --gas 200000 --chain-id $CHAIN_ID --node $NODE -y
```

Revoking blocks future requests. It does not affect a recipient who has already retrieved the key.

### 5. Request access (recipient)

```sh
secretcli tx compute execute $CONTRACT '{"request_access":{}}' \
  --from carl --gas 200000 --chain-id $CHAIN_ID --node $NODE -y

secretcli query compute tx <TXHASH> --node $NODE
```

The response contains the data key and the payload URI and hash. Fetch the ciphertext from the URI, verify its SHA-256 against `content_hash`, and decrypt with the key.

Requests from accounts not on the whitelist are rejected.

## Queries

Public. Returns the owner, Uni-ID, and whether the key and payload have been set:

```sh
secretcli query compute query $CONTRACT '{"config":{}}' --node $NODE
```

All other queries require a signed permit, which proves control of an address without a transaction. Generate one with [secret.js](https://github.com/scrtlabs/secret.js) and submit it as:

```json
{
  "with_permit": {
    "permit": { "params": { ... }, "signature": { ... } },
    "query": { "whitelist": {} }
  }
}
```

| Query | Who can run it |
|---|---|
| `payload` | Owner or a whitelisted recipient |
| `whitelist` | Owner |
| `access_log` | Owner. Takes `page` and `page_size` |
| `is_authorized` | Any signer, about their own address |

## Messages

| Execute | Caller |
|---|---|
| `init_key` | Owner. Once only |
| `set_payload` | Owner |
| `grant` | Owner |
| `revoke` | Owner |
| `request_access` | Whitelisted recipients |
