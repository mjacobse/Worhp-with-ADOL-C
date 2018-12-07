/* Copyright (C) 2018 Luis Lüttgens - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the EPL or GPL license.
 *
 * You should have received a copy of the XYZ license with
 * this file. If not, please write to: luis.luett@googlemail.com
 */

#ifndef CORE_DATA_TYPES_HPP_
#define CORE_DATA_TYPES_HPP_

// strong type for row index
class Row {
 public:
      Row() {}
      explicit Row(int value) : value_(value) {}
      int to_int() const { return value_; }

 private:
      int value_;
};

// strong type for column index
class Col {
 public:
      Col() {}
      explicit Col(int value) : value_(value) {}
      int to_int() const { return value_; }

 private:
      int value_;
};

// strong type for matrix value
class Value {
 public:
      Value() {}
      explicit Value(double value) : value_(value) {}
      double to_double() const { return value_; }

 private:
      double value_;
};

class MatrixLocation {
 public:
      MatrixLocation() {}
      explicit MatrixLocation(Row row, Col col)
         :
         row_(row),
         col_(col)
      {}

      Row getRow()   const {return row_;}
      Col getCol()   const {return col_;}

 private:
      Row row_;
      Col col_;
};

class MatrixEntry {
 public:
      MatrixEntry() {}
      explicit MatrixEntry(Row row, Col col, Value val)
         :
         row_(row),
         col_(col),
         val_(val)
      {}

        Row getRow()   const {return row_;}
        Col getCol()   const {return col_;}
        Value getVal() const {return val_;}

 private:
      Row row_;
      Col col_;
      Value val_;
};

#endif  // CORE_DATA_TYPES_HPP_
