//////////////////////////////////////////////
// Your Name:  Katherine Shanahan
// CS Email:   shanahan@cs.wisc.edu
// CS Login:   shanahan
// 
// CS 368, Spring 2018
// Matrix.hpp
//
//////////////////////////////////////////////

#ifndef LECTURE9_MATRIX_HPP
#define LECTURE9_MATRIX_HPP

#include <iostream>
#include <vector>
#include <stdexcept> // includes runtime_error, which is a subclass of exception

/**
 * @brief This class derives from runtime_error
 *        It is thrown in the case of Non-postive matrix dimensions
 */
class NonPositiveDimensionException : public std::runtime_error {
public:
    // constructor that calls the base class constructor in its initializer list
    NonPositiveDimensionException() : std::runtime_error("dimensions must be positive")
    {// no other code in the constructor
    }
};

/**
 * @brief This class derives from runtime_error
 *        It is thrown in the case of any dimension mismatch
 */
class DimensionMismatchException : public std::runtime_error {
public:
    DimensionMismatchException() : std::runtime_error("dimensions do not match")
    {// no other code in the constructor
    }
};

class RowIndexOutOfBoundsException : public std::runtime_error {
public:
    RowIndexOutOfBoundsException() : std::runtime_error("row index out of bounds")
    {
    }
};

///////////////////////////////////////////////////////
// DECLARATION OF THE MATRIX TEMPLATED CLASS
// YOU WILL NEED TO ADD MORE PROTOTYPES HERE
///////////////////////////////////////////////////////
template<typename T>
class Matrix {
private:
    int rows;
    int cols;
    std::vector<std::vector<T>> data;

public:
    Matrix();
    Matrix(int r, int c);
    void print() const;

    // we need operator[] to return a modifiable L-value, that is, a non-const reference
    // so that we can store a value in a cell of a Matrix
    std::vector<T> & operator[](const int index); 

    // we need a second operator[] in order to make our operator+ below work
    // the operator+ function below needs a const Matrix as a parameter
    // thus this second opertor[] function returns a const reference
    const std::vector<T> & operator[](const int index) const; 

    // operator+ has a const reference parameter, promises not to modify this object
    // and returns a const value
    const Matrix<T> operator+(const Matrix<T> &rhs) const;
    
    // Getter function for Matrix row dimension
    int getRows();
 
    // Getter function for Matrix column dimension
    int getCols();
    
    // Subtracts matrices
    const Matrix<T> operator-(const Matrix<T> &rhs) const;
    
    // Checks equality of matrices
    bool operator==(const Matrix<T> &rhs);
    
    // Checks inequality of matrices
    bool operator!=(const Matrix<T> &rhs);
 
    // Multiplies matrices
    const Matrix<T> operator*(const Matrix<T> &rhs) const;
    
    // Compound adds matrices
    Matrix<T> operator+=(const Matrix<T> &rhs);
    
    // Compound subtracts matrices
    Matrix<T> operator-=(const Matrix<T> &rhs);

    // Compound multiples matrices
    Matrix<T> operator*=(const Matrix<T> &rhs);


    /**
    * @brief overrides the << operator for Matrix<T>
    *        because this is a friend function, and because we have a templated class
    *        we need to define the function inside our class declaration
    * @param os the stream that we are using << with
    * @param obj  the Matrix<T> we are trying to insert into the stream
    * @return a non-const reference, which allows us to chain << operators
    */
    friend std::ostream& operator<<(std::ostream& os, const Matrix<T> &obj) {
        for (auto rowIt = obj.data.begin(); rowIt != obj.data.end(); ++rowIt) {
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                if (colIt != rowIt->end() - 1) {
                    os << *colIt << " ";
                } else {
                    os << *colIt;
                }
            }
            os << std::endl;
        }
        return os;
    }
   
    /*
    * @brief Allows multiplication of a matrix with a scalar value.
    *        Each position is multiplied by the scalar and
    *        added to a new Matrix.
    *
    * @return Matrix<T> with new values of multiplied Matrix
    */
    friend Matrix<T> operator*(const T &scalar, Matrix<T> &rhs) {
        int rows = rhs.rows;
        int cols = rhs.cols;
        Matrix<T> result(rows,cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = scalar * rhs[i][j];
            }
        }
        return result;
    }

    /*
    * @brief Allows multiplication of a matrix with a scalar value.
    *        Along with the other operator* function, multiplication
    *        can take place with a scalar * Matrix or a Matrix * scalar
    *
    * @return Matrix<T> which is an updated Matrix of values.
    *
    */
    friend Matrix<T> operator*(Matrix<T> &rhs, const T &scalar) {
        int rows = rhs.rows;
        int cols = rhs.cols;
        Matrix<T> result(rows,cols);
        result = scalar * rhs;
        return result;
    }
    
    /*
    * @brief operator*= works with the multiplication functions
    *        to be able to compound multiply a Matrix with a scalar
    *
    * @return Matrix<T> which holds the updated values of the Matrix
    *
    */
    friend Matrix<T> operator*=(Matrix<T> &rhs, const T &scalar) {
        rhs = scalar * rhs;
        return rhs;
    }
};

///////////////////////////////////////////////////////
// IMPLEMENTATION OF THE MATRIX TEMPLATED CLASS
// YOU WILL NEED TO ADD MORE FUNCTION DEFINITIONS HERE
///////////////////////////////////////////////////////


/**
 * @brief default Constructor for Matrix<T>
 *        calls the other constructor through an initializer list
 */

template<typename T>
Matrix<T>::Matrix() : Matrix(1,1)
{
}


/**
 * @brief Constructor for Matrix<T>
 *        Matrices store their data in a single vector called data
 *        We assume the matrix will use row-major ordering
 * @param r the number of rows (non-negative)
 * @param c the number of cols (non-negative)
 */
template<typename T>
Matrix<T>::Matrix(int r, int c) {
    if (r <= 0 || c <= 0) {
        throw NonPositiveDimensionException();
    }
    rows = r;
    cols = c;
    data.resize(rows);
    for (int r = 0; r < rows; ++r) {
        data[r].resize(cols);
    }
}

/**
 * @brief prints out to the terminal the elements of this matrix
 *        we assume row-major ordering of the elements
 */
template<typename T>
void Matrix<T>::print() const {
    for (auto rowit = data.begin(); rowit != data.end(); ++rowit) {
        std::vector<T> rowData = *rowit;
        for (auto colit = rowData.begin(); colit != rowData.end(); ++colit) {
            std::cout  << *colit << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * @brief accesses the row index of this Matrix
 *        we assume row-major ordering and a start index of 0
 *        this non-const function returns a non-const reference
 * @param index     the index which corresponds to a vector (row) in the matrix
 * @return a non-const vector reference that is suitable for a Left-value
 */
template<typename T>
std::vector<T> & Matrix<T>::operator[](const int index) {
    if (index < 0 || index >= getRows()) {
       throw RowIndexOutOfBoundsException();
    }
    return data[index];
}

/**
 * @brief accesses the row index of this Matrix
 *        we assume row-major ordering and a start index of 0
 * @param       the index which corresponds to a vector (row) in the matrix
 * @return      a const vector referene that is suitable for a Right-value
 */
template<typename T>
const std::vector<T> & Matrix<T>::operator[](const int index) const {
    if (index < 0 || index >= (*this).rows) {
        throw RowIndexOutOfBoundsException();
    }
    return data[index];
}

/**
 * @brief adds the Matrix<T> on the right side of the + operator to the matrix on the left
 *        is called on the Matrix on the left
 * @param rhs a const reference to the Matrix on the right of the + operator
 * @return a const Matrix that represents the sum
 */
template<typename T>
const Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const {
    Matrix<T> lhs = *this;
    // Checks that the rows and cols have the same dimensions
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        throw DimensionMismatchException();
    }
    int rows = rhs.rows;
    int cols = rhs.cols;
    // Creates a new Matrix to hold the new values
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = lhs[i][j] + rhs[i][j];
        }
    }
    return result;
}

/*
* @brief Gets and returns the number of rows in this Matrix
*
* @return an int as the value of rows
*/
template<typename T>
int Matrix<T>::getRows(){
    return rows;
}

/*
* @brief Gets and return the number of columns in this Matrix
*
* @return and int which is the number of columns
*/
template<typename T>
int Matrix<T>::getCols(){
    return cols;
}

/*
* @brief Subtracts one Matrix from another. Takes the second Matrix
*        as a parameter and operater is called on this Matrix. 
*        Method does not change the Matrices being subtracted.
*
* @return Returns a Matrix which holds the new subtracted matrix.
*/
template<typename T>
const Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const {
    Matrix<T> lhs = *this;
    // Check that they have the same number of rows and cols
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        throw DimensionMismatchException();
    }
    int rows = rhs.rows;
    int cols = rhs.cols;
    // Create new Matrix to hold subtracted values
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = lhs[i][j] - rhs[i][j];
        }
    }
    return result;
}

/*
* @brief Compares two matrices to see if they are the same.
*        If rows or cols numbers don't match, automatically 
*        not equal.
*
* @return bool which is true if they are the same rows/cols
*         and values. Else return false
*/
template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &rhs) {
    Matrix<T> lhs = *this;
    // Check that they have the same number rows and cols
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        return false;
    }
    int rows = rhs.rows;
    int cols = rhs.cols;
    // Check all values to make sure they are the same
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (lhs[i][j] != rhs[i][j]) {
                return false;
            }
        }
    }
    return true;
}

/*
* @brief Determines if matrices are not equal. Return true if they
*        aren't equal. Return false if they are equal. Uses 
*        operator== to determine equality
*
* @return bool that determines equality
*/
template<typename T>
bool Matrix<T>::operator!=(const Matrix<T> &rhs) {
    Matrix<T> lhs = *this;
    // Uses operator== function
    // Returns false if they are equal
    // Returns true if they aren't equal
    if (lhs == rhs){
        return false;
    } else {
        return true;
    }
}

/*
* @brief Multiplication of two matrices. New Matrix will have
*        number of cols of param Matrix and number of rows
*        of called upon Matrix. 
*
* @return Matrix which holds the values of the new Matrix
*/
template<typename T>
const Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const {
    Matrix<T> lhs = *this;
    // Check that cols of this matrix and rows of
    // param matrix are the same.
    if (lhs.cols != rhs.rows) {
         throw DimensionMismatchException();
    }
    int rows = lhs.rows;
    int cols = rhs.cols;
    int sum = 0;
    // Create new matrix that will hold values
    Matrix<T> result(rows, cols);
    // Initialize values with 0 for addition
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = 0;
        }
    }
    // Iterate through rows of this matrix, cols of this
    // matrix and cols of param matrix. 
    // Multiply and add at each position
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < lhs.cols; j++) {
            for (int k = 0; k < cols; k++) {
                result[i][k] += lhs[i][j] * rhs [j][k];
            }
        }
    }
    return result;
}

/*
* @brief Operater+= adds two matrices. Uses operator+ 
*        function for main addition. Directly changes
*        this Matrix.
*
* @return Matrix<T> is the updated Matrix that called
*         this method.
*/
template<typename T>
Matrix<T> Matrix<T>::operator+=(const Matrix<T> &rhs) {
    *this = (*this + rhs);
    return *this;
}

/*
* @brief Operator-= subtracts two matrices. Uses operator-
*        function for main subtraction. Directly changes 
*        this Matrix.
*
* @return Matrix<T> is the updated Matrix that called
*         this method.
*/
template<typename T>
Matrix<T> Matrix<T>::operator-=(const Matrix<T> &rhs) {
    *this = *this - rhs;
    return *this;
}

/*
* @brief Operator*= multiplies two matrices. Uses operator*
*        function for main multiplication. Directly changes
*        this Matrix.
* 
* @return Matrix<T> is the update Matrix that called this
*         method.
*/
template<typename T>
Matrix<T> Matrix<T>::operator*=(const Matrix<T> &rhs) {
    *this = *this * rhs;
    return *this;
}

#endif //LECTURE9_MATRIX_HPP

