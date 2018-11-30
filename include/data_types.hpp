#ifndef DATA_TYPES_HPP
    #define DATA_TYPES_HPP
    
    // strong type for row index
    class Row
    {
        public:
            explicit Row(){};
            explicit Row(int value) : value_(value) {}
            int to_int() const { return value_; }
        private:
            int value_;
    };

    // strong type for column index
    class Col
    {
        public:
            explicit Col(){};
            explicit Col(int value) : value_(value) {}
            int to_int() const { return value_; }
        private:
            int value_;
    };

    // strong type for matrix value
    class Value
    {
        public:
            explicit Value(){};
            explicit Value(double value) : value_(value) {}
            double to_double() const { return value_; }
        private:
            double value_;
    };

    class MatrixEntry_strong
    {
        public:
            explicit MatrixEntry_strong(){};
            explicit MatrixEntry_strong(Row row, Col col,Value val)
            :
                row_(row),
                col_(col),
                val_(val)
            {};

            Row getRow()   const {return row_;}
            Col getCol()   const {return col_;}
            Value getVal() const {return val_;}

        private:
            Row row_;
            Col col_;
            Value val_;
    };

#endif