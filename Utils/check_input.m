function check_input(M)

chi = M.nf - M.ne + M.nv;
assert(chi == 2, 'Algorithm only works for genus 0 surfaces.');
assert(isempty(M.idx_bound), 'Algorithm does not work for surfaces with boundary.');