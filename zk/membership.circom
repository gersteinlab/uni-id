pragma circom 2.1.6;

// Single-variant membership proof. Public input: root. Private: the variant and
// its authentication path. depth must match the tree (23 for a ~6.3M-leaf set).

include "mtree.circom";

component main { public [root] } = Membership(23);
