#include "common.h"


arma::vec to_sigmoid(const arma::vec &outbound_vec, const arma::vec &lower_bound, const arma::vec &upper_bound) {
    arma::vec inbound_vec(outbound_vec.n_elem);
    for (arma::uword i = 0; i < outbound_vec.n_elem; i++) {
        if (outbound_vec(i) > constants::UPP_EXP) {
            inbound_vec(i) = upper_bound(i);
        } else if (outbound_vec(i) < constants::LOW_EXP) {
            inbound_vec(i) = lower_bound(i);
        } else {
            inbound_vec(i) = lower_bound(i) + (upper_bound(i) - lower_bound(i)) * (exp(outbound_vec(i)) / (1 + exp(outbound_vec(i))));
        }
    }

    return inbound_vec;
}

arma::vec from_sigmoid(const arma::vec &inbound_vec, const arma::vec &lower_bound, const arma::vec &upper_bound) {
    arma::vec outbound_vec(inbound_vec.n_elem);
    for (arma::uword i = 0; i < outbound_vec.n_elem; i++) {
        if (inbound_vec(i) < constants::LN_TOL + lower_bound(i)) {
            outbound_vec(i) = constants::LOW_EXP;
        } else if (upper_bound(i) < inbound_vec(i) +  constants::LN_TOL) {
            outbound_vec(i) = constants::UPP_EXP;
        } else {
            outbound_vec(i) = - log((upper_bound(i) - lower_bound(i)) / (inbound_vec(i) - lower_bound(i)) - 1);
        }
    }

    return outbound_vec;
}


arma::mat unique_rows(const arma::mat& m) {

    arma::uvec ulmt = arma::zeros<arma::uvec>(m.n_rows);

    for (arma::uword i = 0; i < m.n_rows; i++) {
        for (arma::uword j = i + 1; j < m.n_rows; j++) {
            if (arma::approx_equal(m.row(i), m.row(j), "absdiff", constants::TOL)) {
                ulmt(j) = 1;
                break;
            }
        }
    }

    return m.rows(find(ulmt == 0));
}


int compute_batch_size(int multi_start, int n_factor) {
    return int((128 * 4 * 3^2) / (multi_start * n_factor^2));
}