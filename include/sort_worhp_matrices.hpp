/* Copyright (C) 2018 Luis Lüttgens - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the EPL or GPL license.
 *
 * You should have received a copy of the XYZ license with
 * this file. If not, please write to: luis.luett@googlemail.com
 */

#ifndef INCLUDE_SORT_WORHP_MATRICES_HPP_
    #define INCLUDE_SORT_WORHP_MATRICES_HPP_

    #include <tuple>
    #include "data_types.hpp"

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

    // typedef of MatrixEntry with weak types
    typedef std::tuple<int, int, double> MatrixEntry_weak;

#endif  // INCLUDE_SORT_WORHP_MATRICES_HPP_
