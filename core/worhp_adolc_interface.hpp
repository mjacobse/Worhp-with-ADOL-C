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

#include <tuple>
#include <vector>
#include <algorithm>
#include <set>

#include "adolc_symbols.hpp"
#include "data_types.hpp"
#include "../user/user_interface.hpp"

namespace worhpAD {

/* intermediate function hiding the lambda in the actual source code
 * Intend is to sort UserDG according to WORHP's desire
 */
template<typename T = MatrixEntry>
auto sortDG(T a = T{}, T b = T{}) {
    return [](auto a, auto b) -> bool {
        // 1) the column index of a is smaller than b's
        // 2) the column indices are equal and the rows ar compared
        // 3) the column index of b is smaller than a's

        if (a.getCol().to_int() <  b.getCol().to_int()) {
            return true;
        } else if (a.getCol().to_int() == b.getCol().to_int()) {
            return a.getRow().to_int() <  b.getRow().to_int();
        } else {
            return false;
        }
    };
}

template<typename T = MatrixEntry>
auto sortHM(size_t sizeHM, T a = T{}, T b = T{}) {
    return [sizeHM](auto a, auto b)->bool {
        // 1) a is in the lower triangle and b on the main diagonal
        // 2) a is on the main diagonal and b in the lower triangle
        // 3) both element are on the main diagonal
        // 4) no  diagonal elements, column-major ordering is performed

        if (a.getRow().to_int() != a.getCol().to_int() &&
            b.getRow().to_int() == b.getCol().to_int()) {
            return true;
        } else if (a.getRow().to_int() == a.getCol().to_int() &&
                   b.getRow().to_int() != b.getCol().to_int()) {
            return false;
        } else if (a.getRow().to_int() == a.getCol().to_int() &&
                   b.getRow().to_int() == b.getCol().to_int()) {
            return a.getRow().to_int() < b.getRow().to_int();
        } else {
        return a.getRow().to_int() + a.getCol().to_int() * sizeHM <
               b.getRow().to_int() + b.getCol().to_int() * sizeHM;
        }
    };
}

void DF_pattern(Workspace* wsp) {
    for (int i = 0; i < wsp->DF.nnz; ++i) {
            wsp->DF.row[i] = adolc::cind_f[i] + 1;
    }
}

void DG_pattern(Workspace* wsp) {
            /*
         * This piece of code is intended to reorder the arrays rind_g and cind_g. ADOL-C
         * returns them in row major order WORHP needs them in column order. The three components
         * of the sparsity pattern rowIdx (rind_g), colIdx (cind_g) and value at this position (jacval).
         * They are gathered in tupels these tupel are then ordered in column-major order.
         * Finally the sparsity pattern of DG is filled with with values.
         */

        std::vector<MatrixLocation> sparseDG;
        for (int i = 0; i < wsp->DG.nnz; ++i) {
            sparseDG.emplace_back(Row(adolc::rind_g[i]), Col(adolc::cind_g[i]));
        }

        std::sort(sparseDG.begin(), sparseDG.end(), worhpAD::sortDG<MatrixLocation>());

        for (int i = 0; i < wsp->DG.nnz; ++i) {
            wsp->DG.row[i] = sparseDG[i].getRow().to_int() +1;
            wsp->DG.col[i] = sparseDG[i].getCol().to_int() +1;
        }
}

void HM_pattern(Workspace* wsp) {
    /*
     * This piece of code is intended to reorder the arrays rind_L and cind_L. ADOL-C
     * returns them in row-major order and the upper triangular matrix.
     * WORHP needs them in a special order. Each diagonal element is bigger than any
     * non-diagonal element. Internally both partitions are sorted in column-major order.
     * 
     * NOTE: WORHP requires a full diagonal independent of this actual sparsity structure,
     * this is considered in this block of code. After the elements are sorted they written to
     * The wsp.HM pattern with tranposed indices.
     */ 

    std::vector<MatrixLocation> sparseHM;
    std::set<int> missingDiagonalElems {};

    for (size_t i = 0; i < user::opt_n; ++i) {
        missingDiagonalElems.insert(i);
    }

    for (int i = 0; i < adolc::nnz_L; ++i) {
        sparseHM.emplace_back(Row(adolc::rind_L[i]), Col(adolc::cind_L[i]));

        if (adolc::rind_L[i] == adolc::cind_L[i]) {
            missingDiagonalElems.erase(adolc::rind_L[i]);
        }
    }

    for (const auto idx : missingDiagonalElems) {
        sparseHM.emplace_back(Row(idx), Col(idx));
    }

    std::sort(sparseHM.begin(), sparseHM.end(), worhpAD::sortHM<MatrixLocation>(user::opt_n));

    for (size_t i = 0; i < sparseHM.size(); ++i) {
        wsp->HM.row[i] = sparseHM[i].getCol().to_int() + 1;
        wsp->HM.col[i] = sparseHM[i].getRow().to_int() + 1;
    }
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

    std::vector<MatrixEntry> sparseDG;
    for (int i = 0; i < wsp->DG.nnz; ++i) {
        sparseDG.emplace_back(Row(adolc::rind_g[i]),
                              Col(adolc::cind_g[i]),
                              Value(adolc::jacval[i]));
    }

    std::sort(sparseDG.begin(), sparseDG.end(), worhpAD::sortDG<>());

    for (int i = 0; i < wsp->DG.nnz; ++i) {
        wsp->DG.val[i] = sparseDG[i].getVal().to_double();
    }
}

void HM(Workspace* wsp, OptVar* opt) {
    sparse_hess(adolc::tag_L, user::opt_n, adolc::reuse_pattern, opt->X, &adolc::nnz_L,
               &adolc::rind_L, &adolc::cind_L, &adolc::hessval, adolc::options_L);

    std::vector<MatrixEntry>  sparseHM;
    std::set<int> missingDiagonalElems {};

    for (size_t i = 0; i < user::opt_n; ++i) {
        missingDiagonalElems.insert(i);
    }

    for (int i = 0; i < adolc::nnz_L; ++i) {
        sparseHM.emplace_back(Row(adolc::rind_L[i]),
                              Col(adolc::cind_L[i]),
                              Value(adolc::hessval[i]));

        if (adolc::rind_L[i] == adolc::cind_L[i]) {
            missingDiagonalElems.erase(adolc::rind_L[i]);
        }
    }

    for (const auto idx : missingDiagonalElems) {
        sparseHM.emplace_back(Row(idx), Col(idx), Value(0));
    }

    std::sort(sparseHM.begin(), sparseHM.end(), worhpAD::sortHM<>(opt->n));

    for (size_t i =0; i < sparseHM.size(); ++i) {
        wsp->HM.val[i] = sparseHM[i].getVal().to_double();
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
