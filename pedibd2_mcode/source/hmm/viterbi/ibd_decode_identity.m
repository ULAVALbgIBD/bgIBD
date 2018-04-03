
function [output_dec_state output_del] = ibd_decode_identity(identity, ob)

ct_total = 0;
t = cputime;

oblist = ob + 1;
[output_dec_state, output_del] = get_viterbi_locus_log(identity, oblist);

ct_total = ct_total + cputime - t

end




















