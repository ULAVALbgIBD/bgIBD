
function new_pairs = reconcile_pairs(pairs, posterior, range)

new_pairs = pairs;
consistency = ibd2_consistent(new_pairs);
new_pairs = ibd2_reconcile(new_pairs, posterior);
new_pairs = ibd1_reconcile(new_pairs, posterior, range);
consistency = ibd2_consistent(new_pairs);


end