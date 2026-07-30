pragma circom 2.1.6;

// Merkle tree templates over Poseidon. No main component; included by the proof
// circuits. Requires circomlib on the include path (-l node_modules/circomlib/circuits).

include "poseidon.circom";
include "switcher.circom";

// leaf = Poseidon(chrom, pos, refCode, altCode)
template VariantLeaf() {
    signal input chrom;
    signal input pos;
    signal input refCode;
    signal input altCode;
    signal output leaf;

    component h = Poseidon(4);
    h.inputs[0] <== chrom;
    h.inputs[1] <== pos;
    h.inputs[2] <== refCode;
    h.inputs[3] <== altCode;
    leaf <== h.out;
}

// One level of a Merkle climb. index is the current node's side (0 left, 1 right).
// parent = Poseidon(left, right) with children ordered by index.
template MerkleLevel() {
    signal input curHash;
    signal input sibling;
    signal input index;
    signal output parent;

    index * (index - 1) === 0;      // index is a bit

    component sw = Switcher();
    sw.sel <== index;
    sw.L   <== curHash;
    sw.R   <== sibling;

    component h = Poseidon(2);
    h.inputs[0] <== sw.outL;
    h.inputs[1] <== sw.outR;
    parent <== h.out;
}

// Membership: recompute the root from a leaf and its authentication path, and
// constrain it to equal the public root. No valid path exists for a non-member,
// so a proof can only be produced for a leaf in the tree.
template Membership(depth) {
    signal input chrom;
    signal input pos;
    signal input refCode;
    signal input altCode;
    signal input pathElements[depth];
    signal input pathIndices[depth];
    signal input root;

    component lh = VariantLeaf();
    lh.chrom   <== chrom;
    lh.pos     <== pos;
    lh.refCode <== refCode;
    lh.altCode <== altCode;

    signal cur[depth + 1];
    cur[0] <== lh.leaf;

    component levels[depth];
    for (var i = 0; i < depth; i++) {
        levels[i] = MerkleLevel();
        levels[i].curHash <== cur[i];
        levels[i].sibling <== pathElements[i];
        levels[i].index   <== pathIndices[i];
        cur[i + 1] <== levels[i].parent;
    }

    root === cur[depth];
}
