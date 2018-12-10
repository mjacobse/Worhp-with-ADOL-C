/* Copyright (C) 2018 Luis Lüttgens - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GPL license.
 *
 * You should have received a copy of the GPL license with
 * this file. If not, please write to: luis.luett@googlemail.com
 */

#ifndef CORE_WORHP_ADOLC_INTERFACE_HPP_
#define CORE_WORHP_ADOLC_INTERFACE_HPP_

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <worhp/worhp.h>

#include <cassert>
#include <vector>
#include <algorithm>

#include "adolc_symbols.hpp"
#include "../user/user_interface.hpp"

namespace worhpAD {

void DF_pattern(Workspace* wsp) {
    for (int i = 0; i < wsp->DF.nnz; ++i) {
            wsp->DF.row[i] = adolc::cind_f[i] + 1;
    }
}

void DG_pattern(Workspace* wsp) {
    for (int i = 0; i < wsp->DG.nnz; ++i) {
        wsp->DG.row[i] = adolc::rind_g[i] + 1;
        wsp->DG.col[i] = adolc::cind_g[i] + 1;
    }
    SortWorhpMatrix(&wsp->DG);
}

void HM_pattern(Workspace* wsp) {
    for (int i = 0; i < adolc::nnz_L; ++i) {
        // transpose
        wsp->HM.row[i] = adolc::cind_L[i] + 1;
        wsp->HM.col[i] = adolc::rind_L[i] + 1;
    }

    std::set<int> missingDiagonalElems {};

    for (size_t i = 0; i < user::opt_n; ++i) {
        missingDiagonalElems.insert(i);
    }

    for (int i = 0; i < adolc::nnz_L; ++i) {
        if (adolc::rind_L[i] == adolc::cind_L[i]) {
            missingDiagonalElems.erase(adolc::rind_L[i]);
        }
    }

    int i = 0;
    for (const auto idx : missingDiagonalElems) {
        wsp->HM.row[adolc::nnz_L + i] = (idx + 1);
        wsp->HM.col[adolc::nnz_L + i] = (idx + 1);
        ++i;
    }

    SortWorhpMatrix(&wsp->HM);
}

void DF(Workspace* wsp, OptVar* opt) {
    double grad_f[user::opt_n];

    gradient(adolc::tag_f, user::opt_n, opt->X, grad_f);

    std::transform(grad_f, grad_f + wsp->DF.nnz, wsp->DF.val,
        [wsp](const double &a) -> double {
            return a * wsp -> ScaleObj;
        });
}

void DG(Workspace* wsp, OptVar* opt) {
    sparse_jac(adolc::tag_g, user::opt_m, user::opt_n, adolc::reuse_pattern, opt->X,
              &wsp->DG.nnz, &adolc::rind_g, &adolc::cind_g, &adolc::jacval, adolc::options_g);

    for (int i = 0; i < wsp->DG.nnz; ++i) {
        assert(wsp->DG.row[i] == adolc::rind_g[wsp->DG.perm[i] - 1] + 1);
        assert(wsp->DG.col[i] == adolc::cind_g[wsp->DG.perm[i] - 1] + 1);
        wsp->DG.val[i] = adolc::jacval[wsp->DG.perm[i] - 1];
    }
}

void HM(Workspace* wsp, OptVar* opt) {
    sparse_hess(adolc::tag_L, user::opt_n, adolc::reuse_pattern, opt->X, &adolc::nnz_L,
               &adolc::rind_L, &adolc::cind_L, &adolc::hessval, adolc::options_L);

    for (int i = 0; i < wsp->HM.nnz; ++i) {
        if (wsp->HM.perm[i] > adolc::nnz_L) {
            // zero diagonal entry, only to obey WORHP's structure
            wsp->HM.val[i] = 0.0;
            continue;
        }
        assert(wsp->HM.row[i] == adolc::cind_L[wsp->HM.perm[i] - 1] + 1);
        assert(wsp->HM.col[i] == adolc::rind_L[wsp->HM.perm[i] - 1] + 1);
        wsp->HM.val[i] = adolc::hessval[wsp->HM.perm[i] - 1];
    }
}

void UserF(Workspace* wsp, OptVar* opt) {
    user::eval_obj(opt->X, opt->F);
    opt->F *= wsp->ScaleObj;
}

void UserG(OptVar* opt) {
    user::eval_constraints(opt->X, opt->G);
}

}  // namespace worhpAD

void UserF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt) {
    worhpAD::UserF(wsp, opt);
}

void UserG(OptVar *opt, Workspace *wsp, Params *par, Control *cnt) {
    worhpAD::UserG(opt);
}

void UserDF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt) {
    worhpAD::DF(wsp, opt);
}

void UserDG(OptVar *opt, Workspace *wsp, Params *par, Control *cnt) {
    worhpAD::DG(wsp, opt);
}

void UserHM(OptVar *opt, Workspace *wsp, Params *par, Control *cnt) {
    worhpAD::HM(wsp, opt);
}

#endif  // CORE_WORHP_ADOLC_INTERFACE_HPP_
