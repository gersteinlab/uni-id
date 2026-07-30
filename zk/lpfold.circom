pragma circom 2.1.6;

// Fold step for the all-variants proof, in Nova-Scotia's interface.
//
//   step_in[2]  = [acc, root]        folding state
//   step_out[2] = [Poseidon(acc, leaf), root]
//   external_inputs[4 + 2*depth]     one variant's fields and Merkle path:
//     [0..3]              chrom, pos, refCode, altCode
//     [4 .. 3+depth]      pathElements
//     [4+depth .. end]    pathIndices
//
// Each step checks the variant is in the allowed-set tree (root) and folds its
// leaf into the accumulator. root is carried through unchanged.

include "poseidon.circom";
include "mtree.circom";

template LPFoldStep(depth) {
    signal input step_in[2];
    signal output step_out[2];

    signal input chrom;
    signal input pos;
    signal input refCode;
    signal input altCode;
    signal input pathElements[depth];
    signal input pathIndices[depth];

    signal acc  <== step_in[0];
    signal root <== step_in[1];

    component lh = Poseidon(4);
    lh.inputs[0] <== chrom;
    lh.inputs[1] <== pos;
    lh.inputs[2] <== refCode;
    lh.inputs[3] <== altCode;
    signal leaf <== lh.out;

    signal cur[depth + 1];
    cur[0] <== leaf;
    component levels[depth];
    for (var i = 0; i < depth; i++) {
        levels[i] = MerkleLevel();
        levels[i].curHash <== cur[i];
        levels[i].sibling <== pathElements[i];
        levels[i].index   <== pathIndices[i];
        cur[i + 1] <== levels[i].parent;
    }
    root === cur[depth];

    component a = Poseidon(2);
    a.inputs[0] <== acc;
    a.inputs[1] <== leaf;

    step_out[0] <== a.out;
    step_out[1] <== root;
}

component main { public [step_in] } = LPFoldStep(23);
